import io
import sys
from math import cos, degrees, radians, sin, sqrt

import matplotlib.pyplot as plt
import paths
from ase.visualize.plot import plot_atoms
from ovito.io.ase import ase_to_ovito
from ovito.modifiers import (
    DeleteSelectedModifier,
    ExpressionSelectionModifier,
    ReplicateModifier,
)
from ovito.pipeline import Pipeline, StaticSource
from ovito.qt_compat import QtCore
from ovito.vis import (
    CoordinateTripodOverlay,
    OpenGLRenderer,
    OSPRayRenderer,
    TachyonRenderer,
    TextLabelOverlay,
    Viewport,
)
from PIL import Image, ImageDraw, ImageFont
from PySide6.QtCore import QBuffer, QIODevice
from Structures import get_structure


def ovito_visualize_and_save(structures, filename, renderer="OpenGL"):
    rendered_images = []

    for structure in structures:
        a, b, c = structure.cell.lengths()
        alpha, beta, gamma = [radians(a) for a in structure.cell.angles()]
        # Convert the ASE structure to an OVITO DataCollection
        data_collection = ase_to_ovito(structure)

        # Create a pipeline and set the data collection as its source
        pipeline = Pipeline(source=StaticSource(data=data_collection))

        # Replicate the structure
        n = 3
        pipeline.modifiers.append(
            ReplicateModifier(num_x=n, num_y=n - 1, num_z=n, adjust_box=False)
        )

        # Select and delete atoms outside the original box
        # expression = (
        #     "Position.X > {:.4f} || Position.Y > {:.4f} || Position.Z > {:.4f}".format(
        #         pipeline.source.data.cell[0][0],
        #         pipeline.source.data.cell[1][1],
        #         pipeline.source.data.cell[2][2],
        #     )
        # )
        # pipeline.modifiers.append(ExpressionSelectionModifier(expression=expression))
        # pipeline.modifiers.append(DeleteSelectedModifier())

        pipeline.compute()

        vis_element = pipeline.compute().particles.vis
        vis_element.scaling = 0.50
        vis_cell = pipeline.source.data.cell.vis
        vis_cell.line_width = 0.05

        # NOTES: We make a coordinate system, and change the vectors
        # to align with the unit cell
        tripod = CoordinateTripodOverlay()
        tripod.size = 0.07
        tripod.offset_x = xof = 0.120
        tripod.offset_y = yof = 0.040
        tripod.axis1_color = [0, 0, 0]
        tripod.axis1_label = f"a = {a:0.2f}"
        tripod.axis2_color = [0, 0, 0]
        tripod.axis2_dir = [cos(gamma), -5.0 * sin(gamma), 0.0]
        tripod.axis2_label = f"b = {b:0.2f}"
        tripod.axis3_color = [0, 0, 0]
        tripod.axis3_label = f"c = {c:0.2f}"
        asy = (cos(alpha) - cos(beta) * cos(gamma)) / sin(gamma)
        asz = 1 / sin(gamma) * sqrt(1 - cos(gamma) ** 2 - asy**2)
        tripod.axis3_dir = [
            cos(beta),
            asy,
            asz,
        ]
        beta_label = TextLabelOverlay(
            text=f"<a>&beta;<a/>={degrees(beta):0.2f}",
            alignment=QtCore.Qt.AlignmentFlag.AlignLeft
            | QtCore.Qt.AlignmentFlag.AlignBottom,
            offset_x=xof + 0.04,
            offset_y=yof + 0.04,
            font_size=0.04,
            text_color=(0, 0, 0),
        )
        tripod.outline_enabled = True
        tripod.outline_color = [0, 0, 0]
        tripod.font_size = 0.6
        tripod.line_width = 0.1
        tripod.outline_enabled = False
        tripod.outline_color = [0, 0, 0]
        # tripod.font = 'Arial,10,-1,5,50,0,0,0,0,0'  # Example font specification

        # Set up the viewport for rendering
        viewport = Viewport()
        viewport.type = Viewport.Type.Ortho
        viewport.camera_pos = (1.21445, -0.0105608, 1.31315)
        viewport.camera_dir = (-0.0990148, 0.990148, -0.0990148)
        viewport.camera_up = (-0.00985234, 0.0985234, 0.995086)
        viewport.fov = n * 3.81262
        viewport.overlays.append(tripod)
        viewport.overlays.append(beta_label)

        pipeline.add_to_scene()

        if renderer == "OpenGL":
            renderer = OpenGLRenderer()
        # Set up the TachyonRenderer
        else:
            renderer = TachyonRenderer()
            renderer.ambient_occlusion = False  # Modify as needed
            renderer.direct_light = True  # Modify as needed
            renderer.direct_light_intensity = 1.0
            renderer.shadows = False

        of = str(
            paths.figures
            / f"Rendered_{structure.info['chemsys']}_{structure.info['structure_name']}.png"
        )
        qimage = viewport.render_image(
            filename=of, renderer=renderer, crop=True, alpha=True
        )

        # Convert QImage to PIL Image
        buffer = QBuffer()
        buffer.open(QIODevice.WriteOnly)
        qimage.save(buffer, "PNG")
        pil_image = Image.open(io.BytesIO(buffer.data()))
        rendered_images.append(pil_image)

        pipeline.remove_from_scene()

    # Combine rendered images into a single row
    total_width = sum(img.width for img in rendered_images)
    max_height = max(img.height for img in rendered_images)
    combined_image = Image.new("RGBA", (total_width, max_height))

    x_offset = 0
    for img in rendered_images:
        combined_image.paste(img, (x_offset, 0))
        x_offset += img.width

    combined_image.save(filename)

    return None


def ase_visualize_and_save(
    structures,
    filename,
    rotation_view="10x,5y,0z",
    scale=1,
    radii=0.5,
    supercell=(1, 1, 1),
):
    # Number of structures
    num_structures = len(structures)

    # Define figure size (Letter width in inches, height adjusted)
    letter_width = 8.5  # Standard Letter width in inches
    height = letter_width / 3  # Adjust height as needed

    # Create a figure with a single row and multiple columns
    fig, axes = plt.subplots(1, num_structures, figsize=(letter_width, height))

    if num_structures == 1:
        axes = [axes]  # Make axes iterable if only one plot

    # Loop over each structure and subplot
    for ax, structure in zip(axes, structures):
        plot_atoms(
            structure.repeat(supercell),
            ax,
            radii=radii,
            rotation=rotation_view,
            scale=scale,
        )

        # Turn off axis
        ax.set_axis_off()

    # Add subplot labels
    for i, ax in enumerate(axes):
        ax.text(
            -0.05,
            0.95,
            f"({chr(97+i)})",
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            horizontalalignment="left",
        )

    plt.tight_layout()
    plt.savefig(filename, dpi=600)


if __name__ == "__main__":
    chemsys = sys.argv[1]
    structure_list = sys.argv[2:]
    structures = [get_structure(chemsys, s) for s in structure_list]
    filename = paths.figures / f"{chemsys}_VisualizedStructures.png"
    ovito_visualize_and_save(structures, str(filename))
