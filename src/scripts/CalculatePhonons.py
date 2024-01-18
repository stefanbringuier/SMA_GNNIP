import paths
from ase.db import connect
from ase.phonons import Phonons
from ase.spacegroup import crystal, get_spacegroup
from ase.spectrum.band_structure import BandStructure
from Config import get_phonon_config
from PlotPhonons import plot_default_phonons
from Structures import get_bandpath

EV_to_THz = 241.799050402293e0


def get_phonons(
    structure,
    potential,
    bandpath,
    supercell=(8, 8, 8),
    doskpts=(30, 30, 30),
    dospts=150,
    displacement=0.03,
    center_refcell=True,
    lorentz_width=4.5e-4,
):
    """
    Calculate phonon properties of a given structure.

    Args:
        structure (ase.Atoms): Atomic structure for phonon calculations.
        potential (tuple): Tuple containing potential model name and calculator object.
        bandpath (ase.dft.kpoints.BandPath): Object defining the band path for phonon dispersion.
        supercell (tuple of int, optional): Supercell dimensions. Defaults to (8, 8, 8).
        doskpts (tuple of int, optional): k-point mesh for density of states. Defaults to (30, 30, 30).
        dospts (int, optional): Number of points for density of states. Defaults to 150.
        displacement (float, optional): Atom displacement in Angstrom for phonon calculation. Defaults to 0.03.
        center_refcell (bool, optional): Whether to center the reference cell. Defaults to True.
        lorentz_width (float, optional): Lorentzian broadening width. Defaults to 4.5e-4.

    Returns:
        tuple: Tuple containing band structure, density of states, and dynamical matrix (ndarray)
    """

    model = potential[0].upper()
    structure_name = structure.info["structure_name"]
    phfolder = str(paths.data / model / structure_name)

    phonons = Phonons(
        structure,
        name=phfolder + "_PhononCalcFiles",
        calc=potential[1],
        supercell=supercell,
        delta=displacement,
        center_refcell=center_refcell,
    )

    phonons.run()
    phonons.read(acoustic=True)
    phonons.clean()

    dyn_mat = phonons.D_N

    # Calculate the band structure for provide path
    # NOTE: No longer need modes here, will use dyn_mat if needed
    omegas = phonons.band_structure(
        bandpath.kpts, modes=False, born=False, verbose=False
    )
    bs = BandStructure(bandpath, energies=omegas[None])
    dos = phonons.get_dos(kpts=doskpts).sample_grid(npts=dospts, width=lorentz_width)

    return bs, dos, dyn_mat


def calculate_phonons(
    dbname,
    chemsys,
    structure_name,
    potential,
    strain=[
        -0.02,
        -0.01,
        0.00,
        0.01,
        0.02,
    ],
    plotting=False,
):
    """
    Calculate phonon properties under various strain conditions.

    Args:
        dbname (str): Name of the database for storing results.
        chemsys (str): The chemical system, e.g. NiTi, PtTi, NiAl3
        structure_name (str): Name of the structure to analyze.
        potential (tuple): Tuple containing potential model name and calculator object.
        strain (list of float, optional): List of fractional strain values to apply. Defaults to pre-defined list.
        plotting (bool, optional): If True, plots the phonon band structure and DOS. Defaults to False.

    Returns:
        None: The function updates the database with calculated data but does not return anything.

    Notes:
       - The database entry for the phonons also contains the real-space dynamical matrix which
         can be used to plot along different band paths or get modes. An example on how to do this
       >>> from ase.db import connect
       >>> from ase.phonons import Phonons
       >>> db = connect("Results.json")
       >>> entry = db.get(model_name="Zhong",structure_name="B2")
       >>> ph = Phonons(entry.toatoms())
       >>> D_N = entry['strain_phonons']['0.00']['dyn_mat']
       >>> ph.D_N = D_N
       >>> ph.get_band_structure(NEW_PATH_YOU_WANT)
    """
    db = connect(dbname, append=True)
    potname = potential[0]
    calculator = potential[1]

    supercell, displacement, bspts = get_phonon_config(structure_name, potential[0])

    # Extract relaxed-structure from database
    entry = db.get(
        chemsys=chemsys, model_name=potname, structure_name=structure_name, relaxed=True
    )
    spg = entry.spacegroup
    cell_relaxed = entry.cell

    structure = crystal(entry.toatoms(), spacegroup=spg)
    # NOTE: Structure cell must be same as cell_relaxed
    structure.set_calculator(calculator)
    structure.info["structure_name"] = entry.structure_name
    bandpath, directions = get_bandpath(structure)

    strain_phonons = {}
    for e in strain:
        mod_structure = structure.copy()
        mod_structure.info["epsilon"] = e
        mod_structure.set_calculator(calculator)
        a, b, c, alpha, beta, gamma = cell_relaxed.cellpar()
        # NOTE: Isotropic strain tensor
        a_s, b_s, c_s = [x * (1.0 + e) for x in (a, b, c)]
        mod_structure.set_cell([a_s, b_s, c_s, alpha, beta, gamma], scale_atoms=True)

        # NOTE: This won't catch all abuse. If structure is monoclinic angle
        # might not be conserved yet still be same spacegroup.
        new_spg = get_spacegroup(mod_structure, symprec=1e-06)
        if new_spg.no != spg:
            ValueError(
                f"Symmetry was broken! Determined strained crystal spacegroup is {new_spg}."
            )

        stress = mod_structure.get_stress()
        volume = mod_structure.get_cell().volume
        bandstructure, dos, dyn_mat = get_phonons(
            mod_structure,
            potential,
            bandpath,
            supercell=supercell,
            displacement=displacement,
        )

        if plotting:
            # Plot raw output bandstructure
            plot_default_phonons(
                bandstructure,
                dos,
                {
                    "chemsys": chemsys,
                    "vol": volume,
                    "strain": e * 100.0,
                    "structure": structure.info["structure_name"],
                    "potname": potname,
                },
                of=str(paths.data / potential[0].upper()),
            )

        eigenvalues = bandstructure.energies * EV_to_THz
        ev_list_by_band = [
            eigenvalues[:, :, i].tolist() for i in range(eigenvalues.shape[2])
        ]

        dos_weights_list = dos.get_weights().tolist()
        dos_energies = dos.get_energies() * EV_to_THz
        dos_energies_list = dos_energies.tolist()

        # Create a dictionary for band structure data
        bandstructure_data = {
            "labels": bandstructure.get_labels()[2],
            "label_locs": bandstructure.get_labels()[1].tolist(),
            "q_vec_labels": directions,
            "q_points": bandstructure.get_labels()[0].tolist(),
            "band_index": {
                str(i): ev_list_by_band[i][0] for i in range(len(ev_list_by_band))
            },
        }

        dos_data = {"weights": dos_weights_list, "energies": dos_energies_list}

        frmt_e = f"{e*100.0:.2f}"
        strain_phonons[frmt_e] = {
            "bandstructure": bandstructure_data,
            "dos": dos_data,
            "dyn_mat": dyn_mat,
            "volume": volume,
            "stress": stress,
            "supercell": supercell,
        }

    db.update(entry.id, data={"strain_phonons": strain_phonons})

    return None
