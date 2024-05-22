from ase.io import write
from ase.io.trajectory import Trajectory


def visualize_phonon_mode(
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
