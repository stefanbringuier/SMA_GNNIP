
import os
import numpy as np
import matplotlib.pyplot as plt
import paths

from ase.phonons import Phonons
from ase.db.core import Database, connect
from ase.spacegroup import get_spacegroup
from ase.spacegroup import crystal
from ase.io.trajectory import Trajectory
from ase.io import write
#from ase.constraints import StrainFilter, UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry

from Structures import *
from Configs import *

EV_to_THz = 241.799050402293e0


def PlotBandstructure(bs, dos, info, ymaxlim=9.0, of="."):
    """
    Plots a single phonon band structure and DOS.
    """
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])

    xticks = bs.get_labels()[1].tolist()
    xlabels = bs.get_labels()[2]

    for xt in xticks:
        ax.axvline(x=xt, color="lightgray")
    ax.axhline(y=0.0, color="black")

    q_points = bs.get_labels()[0]
    energies_THz = EV_to_THz * bs.energies
    for n in range(len(energies_THz[0, 1, :])):
        omega = energies_THz[0, :, n]
        ax.plot(q_points, omega, "k-", lw=2)

    ax.set_title(
        f"Volume: {info['vol']:.3f}Å^3\
        Strain: {info['strain']:.2f}% \
        Potential: {info['potname']}",
        fontsize=10,
    )

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    ax.set_xlim((0.0, q_points[-1]))
    ax.set_ylim((-2.0, ymaxlim))
    ax.set_ylabel("Frequency ($\mathrm{THz}$)", fontsize=22)

    dosax.fill_between(
        dos.get_weights(),
        EV_to_THz * dos.get_energies(),
        y2=0,
        color="grey",
        edgecolor="k",
        lw=1,
    )

    ymaxlim = np.max(energies_THz) * 1.05
    dosax.set_ylim((-2.0, ymaxlim))
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS", fontsize=18)
    fig.savefig(
        os.path.join(
            of,
            f"{info['structure']}_{info['potname']}_{info['strain']:.3f}_Phonons.png",
        )
    )
    fig.clear()

    return None


def WriteModeWisual(
    fileprefix,
    phonons,
    structure,
    point="M",
    branches=[0, 1],
    repeat=(8, 8, 8),
    amp_scale=25.7e-3,
    mode_branches=[
        ("TA1", "mode.0.traj"),
        ("TA2", "mode.1.traj"),
    ],
    image_rotation="-36x,26.5y,-25z",
    pts=100,
    of=".",
):
    """
    Creates trajectory files and then gif images of modes at a given symmetry point.

    NOT TESTED with showyourwork flow
    """

    # Write modes for specific q-vector to trajectory files:
    bandpath = structure.cell.bandpath(npoints=pts)
    p = bandpath.special_points[point]
    phonons.write_modes(
        [l / 2 for l in p],
        branches=branches,
        repeat=repeat,
        kT=amp_scale,
        center=True,
    )

    for branch in mode_branches:
        trajfile = phonons.name + "." + branch[1]
        with Trajectory(trajfile, "r") as traj:
            write(
                os.path.join(of, f"{fileprefix}.{point}Point.Mode.{branch[0]}.gif"),
                traj,
                interval=50,
                rotation=image_rotation,
            )
            write(
                os.path.join(of, f"{fileprefix}.{point}Point.Mode.{branch[0]}.extxyz"),
                traj,
            )

    # clean-up write modes
    for branch in mode_branches:
        os.remove(phonons.name + "." + branch[1])


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
    displacement is in Angstrom and is atom movement.
    """

    phfolder = (
        str(paths.data / potential[0].upper() / structure.info["structure_name"])
        + "-phonon-calc-files"
    )
    
    phonons = Phonons(
        structure,
        name=phfolder,
        calc=potential[1],
        supercell=supercell,
        delta=displacement,
        center_refcell=center_refcell,
    )

    phonons.run()
    phonons.read(acoustic=True)
    phonons.clean()

    bs = phonons.get_band_structure(bandpath, verbose=False)
    dos = phonons.get_dos(kpts=doskpts).sample_grid(npts=dospts, width=lorentz_width)

    return bs, dos, phonons


def calculate_phonons(
    dbname,
    structure_name,
    potential,
    strain=[
        -0.11000,
        -0.09375,
        -0.07750,
        -0.06125,
        -0.04500,
        -0.02875,
        -0.01250,
        0.00375,
        0.02000,
    ],
):
    db = connect(dbname, append=True)
    potname = potential[0]
    calculator = potential[1]

    supercell, displacement, bspts = GetPhononCalcConfig(structure_name, potential[0])

    entry = db.get(model_name=potname, structure_name=structure_name)
    spg = entry.spacegroup
    cell_relaxed = entry.cell

    structure = crystal(entry.toatoms(), spacegroup=spg)
    # NOTE: Structure cell must be same as cell_relaxed
    structure.set_calculator(calculator)
    structure.info["structure_name"] = entry.structure_name
    bandpath,directions = get_bandpath(structure)

    # Ideally we want to enfoce these filters. We cannot break symmetry
    # and strain tensor must comply. But this doesn't work.
    # structure.set_constraint(FixSymmetry(structure, symprec=1.0e-4))
    # sf = StrainFilter(structure)

    strain_phonons = {}
    for e in strain:
        mod_structure = structure.copy()
        mod_structure.set_calculator(calculator)
        a,b,c,alpha,beta,gamma = cell_relaxed.cellpar()
        # NOTE: Isotropic strain
        a_s,b_s,c_s = [ x*(1.0 + e) for x in (a,b,c) ]
        mod_structure.set_cell([a_s,b_s,c_s,alpha,beta,gamma],scale_atoms=True)

        # NOTE: This won't catch all abuse. If structure is monoclinic angle
        # might not be conserved yet still be same spacegroup.
        new_spg = get_spacegroup(mod_structure,symprec=1e-06)
        if new_spg.no != spg:
            ValueError(
                f"Symmetry was broken! Determined strained crystal spacegroup is {new_spg}."
            )

        stress = mod_structure.get_stress()
        volume = mod_structure.get_cell().volume
        bandstructure, dos, ph = get_phonons(
            mod_structure,
            potential,
            bandpath,
            supercell=supercell,
            displacement=displacement,
        )

        # Plot raw output bandstructure
        PlotBandstructure(
            bandstructure,
            dos,
            {
                "vol": volume,
                "strain": e * 100.0,
                "structure": structure.info["structure_name"],
                "potname": potname,
            },
            of=str(paths.data / potential[0].upper() ),
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
            "volume": volume,
            "stress": stress,
        }

    db.update(entry.id, data={"strain_phonons": strain_phonons})

    # write_mode_visual(f"{potname}_a={a:.3f}_NiTi-B2", ph, B2,of=of)
    return None
