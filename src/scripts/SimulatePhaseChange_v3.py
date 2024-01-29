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
from scipy.spatial import ConvexHull, Voronoi

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
    write(
        filename,
        traj,
        format="extxyz",
    )
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
                "Atomic Volume (Å³)",
            ]
            writer.writerow(csv_headers)
        else:
            writer.writerow(line)


def custom_md_logger(dynamics, structure, file_path, interval=10):
    """
    Custom logger for MD simulations. Writes simulation data to a CSV file.

    Parameters:
    dynamics (ASE dynamics object): The dynamics object from which to extract simulation data.
    structure (ASE Atoms object): The Atoms object representing the simulated system.
    file_path (str): Path to the CSV file for logging.
    interval (int): The number of steps between logging data.
    """

    def write_data():
        if dynamics.nsteps % interval == 0:
            natoms = structure.get_global_number_of_atoms()
            time = dynamics.get_time() / (1000 * fs)  # Convert to ps
            etotal = structure.get_total_energy()
            pote_atom = structure.get_potential_energy() / natoms
            temp = structure.get_temperature()
            volume = structure.get_volume()
            density = (structure.get_masses().sum() / volume) * (CM**3 / GRAMS)
            cell_lengths = structure.get_cell().lengths()
            cell_angles = structure.get_cell().angles()
            stress = structure.get_stress(include_ideal_gas=True) / GPa

            structure.atomic_volume = get_atomic_volumes(structure)
            avg_atomic_vol = np.average(structure.atomic_volume)

            # Calculate NiTi order parameter if applicable
            if "NiTi" in chemsys:
                orderparam = get_orderparam_niti(structure)
                structure.order_parameter = orderparam
                mean_op = np.average(orderparam)
            else:
                mean_op = 0

            line = (
                [time, etotal, pote_atom, temp, volume, density]
                + list(cell_lengths)
                + list(cell_angles)
                + list(stress)
                + [mean_op, avg_atomic_vol]
            )
            write_to_csv(file_path, line=line, mode="a")

    return write_data


def get_atomic_volumes(structure, cutoff=2.825):
    cell = structure.get_cell()
    positions = structure.get_positions()

    # Build neighbor list
    cutoffs = [cutoff] * len(structure)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(structure)

    # Extend positions to include periodic images
    extended_positions = positions.copy()
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                if not (i == 0 and j == 0 and k == 0):
                    displacement = i * cell[0] + j * cell[1] + k * cell[2]
                    extended_positions = np.vstack(
                        [extended_positions, positions + displacement]
                    )

    # Perform Voronoi tessellation
    voronoi = Voronoi(extended_positions)

    volumes = []
    for atom_index in range(len(structure)):
        # Identify vertices of the Voronoi cell for the atom
        region_index = voronoi.point_region[atom_index]
        vertices = voronoi.vertices[voronoi.regions[region_index]]

        # Calculate volume of the Voronoi cell
        if len(vertices) > 0 and not -1 in vertices:
            volumes.append(ConvexHull(vertices).volume)
        else:
            volumes.append(
                np.nan
            )  # Can't calculate volume for cells with missing vertices

    return volumes


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
    super_cell=(8, 8, 8),
    initial_temp=500.00,
    final_temp=150.00,
    ramp_rate=0.001,  # K/fs
    timesteps=250_000,
    dt=1.0,
    out_interval=1_000,
    pfactor=1000 * fs,
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
    structure.new_array("atomic_volume", np.zeros(natoms))

    # CSV file setup
    csv_file_path = (
        paths.data
        / model_name.upper()
        / f"{chemsys}_{structure_name}_{model_name}_MD-NPT_DynamicTemp.csv"
    )

    # Write headers to CSV file
    write_to_csv(csv_file_path, mode="w")

    # Step 1: Equilibrate to starting temperature
    MaxwellBoltzmannDistribution(structure, temperature_K=initial_temp)
    equil_dynamics = NPT(
        structure,
        dt * fs,
        temperature_K=initial_temp,
        externalstress=1.013 * bar,
        ttime=75 * fs,
        pfactor=pfactor,
        mask=np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    )
    equil_dynamics.run(25_000)

    # Corrected Dynamic temperature ramping setup
    current_set_temp = initial_temp
    total_temp_change = final_temp - initial_temp
    total_simulation_time = timesteps * dt  # Total simulation time in fs
    actual_ramp_rate = total_temp_change / total_simulation_time  # K/fs
    temp_change_per_step = actual_ramp_rate * dt  # Temperature change per timestep

    # Initialize MD simulation with NPT ensemble
    dynamics = NPT(
        structure,
        dt * fs,
        temperature_K=current_set_temp,
        externalstress=1.013 * bar,
        ttime=75 * fs,
        pfactor=pfactor,
        mask=np.array([[1, 0, 1], [0, 1, 0], [1, 0, 1]]),
    )

    # Function to update temperature at each step
    def update_temperature():
        nonlocal current_set_temp
        current_set_temp += temp_change_per_step  # Incremental temperature change
        dynamics.set_temperature(temperature_K=current_set_temp)

    # Attach the temperature update function to the dynamics
    temp_update_interval = max(1, int(dt / (actual_ramp_rate / ramp_rate)))
    dynamics.attach(update_temperature, interval=temp_update_interval)

    # Attach Trajectory setup
    traj_path = (
        paths.data
        / model_name.upper()
        / f"{chemsys}_{structure_name}_{model_name}_MD-NPT_DynamicTemp.traj"
    )
    trajectory = Trajectory(traj_path, "w", structure)

    dynamics.attach(trajectory, interval=out_interval)

    # MDLogger for logging (can be customized as needed)
    logger_func = custom_md_logger(
        dynamics, structure, csv_file_path, interval=out_interval
    )
    dynamics.attach(logger_func)

    # Run the dynamics
    dynamics.run(timesteps)

    # Post-simulation tasks
    trajectory.close()

    write_exyz(traj_path)


if __name__ == "__main__":
    dbname = paths.data / sys.argv[1]
    model = sys.argv[2]
    chemsys = sys.argv[3]
    structure = sys.argv[4]
    # pfactor = sys.argv[5]
    simulate_dynamics(
        dbname, model, chemsys=chemsys, structure_name=structure, timesteps=10
    )
