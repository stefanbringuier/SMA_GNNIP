import sys

import numpy as np
import paths
from ase.db import connect
from ase.io import Trajectory, read, write
from ase.md import MDLogger
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.units import GPa, bar, fs, kg, m
from Calculators import get_ase_calculator


def enforce_uppertriangular_cell(cell, threshold=1.0e-10):
    zero_mask = np.abs(cell) < threshold
    cell[zero_mask] = 0

    # Calculate the upper triangular matrix
    new_cell = np.triu(cell)
    return new_cell


def write_exyz(traj_filepath):
    traj = read(traj_filepath, index=":", format="traj")
    filename = "".join([str(traj_filepath).split(".")[0], ".extxyz"])
    write(filename, traj, format="extxyz")
    return None


class CustomMDLogger(MDLogger):
    """
    `cumulative_time` is in fs
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cumulative_time = 0.0 / (1000 * fs)
        self.density = (0.0 * 1.0e-3 * kg) / ((1.0e-2 * m) ** 3)

    def __call__(self):
        # Obtain potential and kinetic energy, and temperature
        epot = self.atoms.get_potential_energy()
        ekin = self.atoms.get_kinetic_energy()
        temp = self.atoms.get_temperature()
        global_natoms = self.atoms.get_global_number_of_atoms()

        # Normalize per atom if necessary
        if self.peratom:
            epot /= global_natoms
            ekin /= global_natoms

        # Update cumulative time based on the condition
        if self.dyn is not None:
            current_time = self.dyn.get_time() / (1000 * fs)
            try:
                factor = (self.cumulative_time + current_time) % self.cumulative_time
            except ZeroDivisionError:
                factor = 0.0e0

            if not factor != 0.0e0:
                display_time = self.cumulative_time + current_time
            else:
                display_time = self.cumulative_time

        dat = (display_time, epot + ekin, epot, ekin, temp)

        # Include stress in the log if requested
        if self.stress:
            stress_data = tuple(self.atoms.get_stress(include_ideal_gas=True) / GPa)
            dat += stress_data

        # self.add_density()

        # dat += (self.density)

        # Write data to the logfile
        self.logfile.write(self.fmt % dat)
        self.logfile.flush()

    def update_time(self, time):
        self.cumulative_time = time / (1000 * fs)

    def add_density(self):
        grams = 1.0e-3 * kg
        cm = 1.0e-2 * m
        mass = self.atoms.get_masses().sum() * grams
        volume = self.atoms.cell.volume * cm**3
        if volume > 0:  # Check to avoid division by zero
            self.density = mass / volume
        else:
            self.density = 0.0  # or handle this case as needed
        # Add to formatter
        self.hdr += " Density [g/cm^3]"
        self.fmt += " %10.6f\n"


def simulate_dynamics(
    dbname,
    model_name,
    formula="NiTi",
    structure_name="B2",
    super_cell=(5, 5, 5),
    initial_temp=425.00,
    final_temp=125.00,
    temp_step=10.00,
    timesteps=10_000,
    dt=1.0,
    out_interval=5000,
):
    """
    Simulate the NPT ensemble starting from a starting temperature to a final temperature.

    NOTES:
    1. The ramping down temperature doesnt seem to work as with LAMMPS, this probably needs to be
    done in a quasi-transient fashion.

    """
    db = connect(dbname)
    entry = db.get(formula=formula, structure_name=structure_name, model_name=model)
    structure = entry.toatoms()
    cell = enforce_uppertriangular_cell(structure.get_cell())
    structure.set_cell(cell, scale_atoms=True)
    structure = structure.repeat(super_cell)

    # Get calculator
    calculator = get_ase_calculator(model=model_name)
    structure.set_calculator(calculator)

    # Step 1: Equilibrate
    # Assign distribution of momenta to target temperature
    MaxwellBoltzmannDistribution(structure, temperature_K=initial_temp)
    equil_dynamics = NPT(
        structure,
        dt * fs,
        temperature_K=initial_temp,
        externalstress=1.013 * bar,
        ttime=100 * fs,
        pfactor=1000 * fs,
        mask=(1, 1, 1),
    )

    equil_dynamics.run(10)

    # Logging
    log_path = (
        paths.data
        / model_name.upper()
        / f"{formula}_{structure_name}_{model_name}_MD-NPT.log"
    )

    # Logger, tie to equil_dynamics but not used
    logger = CustomMDLogger(
        equil_dynamics, structure, log_path, header=True, stress=True, mode="w"
    )

    # Step 2: Phase transformation
    # Doing a quasi-equilibrium dynamics, hold 10 ps at NPT in 10 K decrements
    # Trajectory
    traj_path = (
        paths.data
        / model_name.upper()
        / f"{formula}_{structure_name}_{model_name}_MD-npt.traj"
    )
    trajectory = Trajectory(traj_path, "w", structure)

    cumalitive_time = 0.0e0 * fs
    for t in np.arange(initial_temp, final_temp, -temp_step):
        # Mask should only allow for x, y, z, and xz componenets to change.
        dynamics = NPT(
            structure,
            dt * fs,
            temperature_K=t,
            externalstress=1.013 * bar,
            ttime=100 * fs,
            pfactor=1000 * fs,
            mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]]),
        )

        dynamics.attach(logger, interval=out_interval)
        dynamics.attach(trajectory, interval=out_interval)
        dynamics.run(timesteps)

        logger.update_time(timesteps * dt * fs)
        cumalitive_time += timesteps * dt * fs

    trajectory.close()

    write_exyz(traj_path)


if __name__ == "__main__":
    dbname = paths.data / sys.argv[1]
    model = sys.argv[2]
    formula = sys.argv[3]
    structure = sys.argv[4]
    simulate_dynamics(dbname, model, formula=formula, structure_name=structure)
