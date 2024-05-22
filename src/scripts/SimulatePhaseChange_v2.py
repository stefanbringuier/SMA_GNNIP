import csv
import sys

import numpy as np
import paths

# from ase.parallel import world
# from ase.md import MDLogger
from ase.db import connect
from ase.io import Trajectory, read, write
from ase.md.npt import NPT
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.neighborlist import NeighborList
from ase.units import GPa, bar, fs, kg, m
from Calculators import get_ase_calculator

GRAMS = kg * (1 / 1000)
CM = m * (1 / 100)


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


def write_to_csv(file_path, line="", mode="a"):
    """
    Writes data to a CSV file.

    Parameters:
    file_path (str): The path of the file to write to.
    data (list): The data to write. If it's a list of lists, multiple rows will be written.
    mode (str): The file opening mode ('w' for write, 'a' for append).
    """
    with open(file_path, mode, newline="") as file:
        writer = csv.writer(file)
        if mode == "w":
            csv_headers = [
                "Cumulative Time (ps)",
                "Total Energy (eV)",
                "Potential Energy per Atom (eV/atom)",
                "Temperature (K)",
                "Volume (Å³)",
                "Density (g/cm³)",
                "Lx (Å)",
                "Ly (Å)",
                "Lz (Å)",
                "α (deg)",
                "β (deg)",
                "γ (deg)",
                "σ_xx (GPa)",
                "σ_yy (GPa)",
                "σ_zz (GPa)",
                "σ_yz (GPa)",
                "σ_xz (GPa)",
                "σ_xy (GPa)",
                "NiTi_Order",
            ]
            writer.writerow(csv_headers)
        else:
            writer.writerow(line)


def get_orderparam_niti(
    structure,
    cutoff=2.725,
    d0B19=2.460232,
    d1B19=2.646524,
    d0B2=2.611067,
    d1B2=2.611067,
):
    """
    Calculate the order parameter for NiTi alloys based on the method described by Mutter et al.

    This function computes the order parameter for each atom in a provided ASE Structure object. The order parameter
    calculation is specific to NiTi alloys and is based on the distances to the nearest neighbors, taking into account
    periodic boundary conditions. The heuristic constants used in the calculation can be adjusted.

    Parameters:
    structure (ASE Structure): The Structure object representing the atomic structure.
    cutoff (float, optional): The cutoff distance for neighbor search. Default is 2.725 Ångströms.
    d0B19 (float, optional): Heuristic constant for NiTi, characteristic distance. Default is 2.460232 Ångströms.
    d1B19 (float, optional): Heuristic constant for NiTi, characteristic distance. Default is 2.646524 Ångströms.
    d0B2 (float, optional): Heuristic constant for NiTi, characteristic distance. Default is 2.611067 Ångströms.
    d1B2 (float, optional): Heuristic constant for NiTi, characteristic distance. Default is 2.611067 Ångströms.

    Returns:
    numpy.ndarray: An array of order parameters for each atom in the Structure object.

    Raises:
    ValueError: If the number of neighbors for any atom is less than 8, as required for the calculation.

    Notes:
    - The function is specific to NiTi alloys and may not be applicable to other materials.
    - The default heuristic constants are based on the study by Mutter et al. and should be adjusted if different
      parameters are required.
    - Ensure the Structure object contains correct periodic boundary conditions for accurate calculations.
    """
    orderp = np.zeros(len(structure))

    # Define the cutoff radius for neighbor search
    cutoffs = [cutoff] * len(structure)

    # Build neighbor list
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(structure)

    # Calculate order parameter
    for i in range(len(structure)):
        indices, offsets = nl.get_neighbors(i)
        itype = structure[i].symbol

        dist = []
        for j, offset in zip(indices, offsets):
            jtype = structure[j].symbol
            if itype == jtype:  # Only Ni-Ti bonds
                continue

            # Consider periodic boundary conditions
            rij = (
                structure.positions[j]
                + np.dot(offset, structure.get_cell())
                - structure.positions[i]
            )
            dij = np.linalg.norm(rij)
            dist.append(dij)

        # Sort distances and check the number of neighbors
        dist.sort()
        if len(dist) < 8:
            raise ValueError("Not enough neighbors for atom {}".format(i))

        d0 = np.sum(dist[:6]) / 6.0
        d1 = np.sum(dist[6:8]) / 2.0
        chi = (d0 * (d1B2 + d1B19) - d1 * (d0B2 + d0B19)) / (d0B2 * (d0B19 - d1B19))
        orderp[i] = chi

    return orderp


def simulate_dynamics(
    dbname,
    model_name,
    chemsys="NiTi",
    structure_name="B2",
    super_cell=(3, 3, 3),
    initial_temp=550.00,
    final_temp=125.00,
    temp_step=10.00,
    timesteps=25_000,
    dt=1.0,
    out_interval=5000,
    pfactor=100 * fs,
    ramp_direction=None,  # New kwarg to determine ramp direction
):
    db = connect(dbname)
    entry = db.get(
        chemsys=chemsys, structure_name=structure_name, model_name=model_name
    )
    structure = entry.toatoms()
    cell = enforce_uppertriangular_cell(structure.get_cell())
    structure.set_cell(cell, scale_atoms=True)
    structure = structure.repeat(super_cell)
    natoms = structure.get_global_number_of_atoms()
    calculator = get_ase_calculator(model=model_name)
    structure.set_calculator(calculator)
    if chemsys == "NiTi":
        structure.new_array("order_parameter", np.zeros(natoms))

    # Determine ramp direction based on initial and final temperatures
    if ramp_direction is None:
        ramp_direction = "RampUp" if final_temp > initial_temp else "RampDown"

    # CSV file setup
    csv_file_path = (
        paths.data
        / model_name.upper()
        / f"{chemsys}_{structure_name}_{model_name}_MD-NPT_{ramp_direction}.csv"
    )

    # Write headers to CSV file
    write_to_csv(csv_file_path, mode="w")

    # Equilibration step
    MaxwellBoltzmannDistribution(structure, temperature_K=initial_temp)
    equil_dynamics = NPT(
        structure,
        dt * fs,
        temperature_K=initial_temp,
        externalstress=1.013 * bar,
        ttime=5 * fs,
        pfactor=pfactor,
        mask=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    )
    equil_dynamics.run(10_000)

    # Trajectory setup
    traj_path = (
        paths.data
        / model_name.upper()
        / f"{chemsys}_{structure_name}_{model_name}_MD-NPT_{ramp_direction}.traj"
    )
    trajectory = Trajectory(traj_path, "w", structure)

    # Temperature stepping
    temperature_range = np.arange(
        initial_temp,
        final_temp,
        temp_step if ramp_direction == "RampUp" else -temp_step,
    )
    cumulative_time = 0.0e0 * fs
    for t in temperature_range:
        dynamics = NPT(
            structure,
            dt * fs,
            temperature_K=t,
            externalstress=1.013 * bar,
            ttime=5 * fs,
            pfactor=pfactor,
            mask=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        )
        dynamics.attach(trajectory, interval=out_interval)
        dynamics.run(timesteps)

        # Update cumulative time and other properties
        cumulative_time += (
            timesteps * dt * fs / (1000 * fs)
        )  # Convert to ps, formatted to 2 decimal places
        etotal = round(structure.get_total_energy(), 4)
        pote_atom = round(structure.get_potential_energy() / natoms, 4)
        volume = round(structure.get_volume(), 2)  # Volume rounded to 2 decimal places
        lx, ly, lz = map(lambda x: round(x, 2), structure.get_cell().lengths())
        alpha, beta, gamma = map(lambda x: round(x, 2), structure.get_cell().angles())
        temperature = round(structure.get_temperature(), 2)
        density = round((structure.get_masses().sum() / volume) * (CM**3 / GRAMS), 3)
        stress = structure.get_stress(include_ideal_gas=True) / GPa
        stress_formatted = [round(s, 3) for s in stress]

        properties = (
            cumulative_time,
            etotal,
            pote_atom,
            temperature,
            volume,
            density,
            lx,
            ly,
            lz,
            alpha,
            beta,
            gamma,
        ) + tuple(stress_formatted)

        if chemsys == "NiTi":
            orderparam = get_orderparam_niti(structure)
            structure.order_parameter = orderparam
            mean_op = np.average(orderparam)
            properties += (mean_op,)
        else:
            properties = +(0,)

        write_to_csv(csv_file_path, line=properties)

    trajectory.close()
    write_exyz(traj_path)


if __name__ == "__main__":
    dbname = paths.data / sys.argv[1]
    model = sys.argv[2]
    chemsys = sys.argv[3]
    structure = sys.argv[4]
    # pfactor = sys.argv[5]
    simulate_dynamics(dbname, model, chemsys=chemsys, structure_name=structure)
