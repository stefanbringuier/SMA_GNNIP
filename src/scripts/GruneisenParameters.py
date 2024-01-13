import paths
import sys
import numpy as np
import numpy.polynomial.polynomial as poly

# TODO: Move to ConfigPlots and use for all.
import matplotlib.pyplot as plt
import matplotlib.style as style
# Tufte-inspired style
#style.use('seaborn-whitegrid')
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['figure.dpi'] = 300


from ase.db import connect

from TableConfig import ORDER

# Function to apply the signed log transformation
def log_transform(val):
    return np.sign(val) * np.log10(abs(val) + 1)

def find_index(data, label):
    """
    locates the index for the q-point value that is closest to the point of interest.

    Arguments:
        data (dict): Phonon data
        label (str): Symmetry point label

    Returns:
        int: index of the closest matched occurance.
    """
    if label not in data["bandstructure"]["labels"]:
        return "Label not found"

    label_index = data["bandstructure"]["labels"].index(label)
    label_loc = data["bandstructure"]["label_locs"][label_index]
    q_points = data["bandstructure"]["q_points"]
    closest_index = min(
        range(len(q_points)), key=lambda i: abs(q_points[i] - label_loc)
    )

    return closest_index


def gruneisen(phonons, n=3):
    """Calculates the derivative of the frequency with respect to volume for all q-points.

    This function computes the derivative for the Gruneisen parameter. It fits a polynomial of degree `n`
    to the frequency-vs-volume data for each mode at each wave vector and then evaluates the derivative
    of this polynomial at the equilibrium volume.

    Args:
        phonons (dict): A dictionary with keys being the strain values (as strings, e.g., "0.00") and values
            being dictionaries containing 'volume' and 'bandstructure' data.
        n (int, optional): The degree of the polynomial to fit to the frequency-vs-volume data. Defaults to 3.

    Returns:
        np.ndarray: A 2D array where each element [i, j] is the derivative of the frequency of the i-th mode
        at the j-th wave vector with respect to the volume, evaluated at the equilibrium volume.

    """

    volumes = []
    ef_list = []
    for k in phonons.keys():
        volumes.append(phonons[k]["volume"])
        bands = []
        for band in phonons[k]["bandstructure"]["band_index"].values():
            bands.append(band)
        ef_list.append(bands)
    ef = np.array(ef_list)

    v0 = phonons["0.00"]["volume"]
    domega_dv0 = np.zeros((ef.shape[1], ef.shape[2]))
    for i in range(ef.shape[1]):
        for j in range(ef.shape[2]):
            # Fit a polynomial to the i-th mode
            coefs = poly.polyfit(volumes, ef[:, i, j], n)
            # Differentiate the polynomial
            derivative = poly.polyder(coefs)
            # Evaluate the derivative at the equilibrium volume
            domega_dv0[i, j] = poly.polyval(v0, derivative)

    return domega_dv0


def mode_gruneisen(phonons, q, n=3):
    """
    Calculate the specific mode Gruneisen (see function `gruneisen` for details).
    """

    volumes = []
    ef_list = []
    for k in phonons.keys():
        volumes.append(phonons[k]["volume"])
        mode_index = find_index(phonons[k], q)
        mode_bands = []
        for band in phonons[k]["bandstructure"]["band_index"].values():
            mode_bands.append(band[mode_index])
        ef_list.append(mode_bands)
    ef = np.array(ef_list)

    v0 = phonons["0.00"]["volume"]
    domega_dv0 = np.zeros(ef.shape[1])
    for i in range(ef.shape[1]):
        coefs = poly.polyfit(volumes, ef[:, i], n)
        derivative = poly.polyder(coefs)
        # Evaluate the derivative at the equilibrium volume
        domega_dv0[i] = poly.polyval(v0, derivative)

    return domega_dv0


def calculate_mode_gruneisen(dbentry, q="M"):
    """Calculate M-mode Gruneisen parameter.

    The Gruneisen parameter \(\gamma\) is given by:

    \[
    \gamma = -\frac{V}{\omega}\left(\frac{\partial \omega}{\partial V}\right)_{T, P}
    \]

    where \( V \) is the volume, \( \omega \) is the phonon frequency, and the derivative is taken at constant temperature and pressure.

    Arguments:
        dbentry (ASE.db.row): The specific database entry.
        q (str): The symmetry point to calculate at

    Returns:
        dict: Specific phonon data for ASE database entry.
        float: The mode Grunesien parameter

    """

    data = dbentry.data
    phonons = data["strain_phonons"]
    v0 = phonons["0.00"]["volume"]
    mode_index = find_index(phonons["0.00"], q)
    omega0 = [
        b[mode_index] for b in phonons["0.00"]["bandstructure"]["band_index"].values()
    ]
    gamma = -v0 / np.array(omega0) * mode_gruneisen(phonons, q)
    return omega0, gamma


def calculate_gruneisen(dbname, model_name, structure_name):
    """Calculate the Gruneisen parameters.

    The Gruneisen parameter \(\gamma\) is given by:

    \[
    \gamma = -\frac{V}{\omega}\left(\frac{\partial \omega}{\partial V}\right)_{T, P}
    \]

    where \( V \) is the volume, \( \omega \) is the phonon frequency, and the derivative is taken at constant temperature and pressure.

    Arguments:
        dbname (str): The path to the ASE database
        model_name (str): The model to calculate the Gruneisen parameters for, e.g., Zhong.
        structure_name (str): The structure name, e.g., B2.

    Returns:
        dict: Specific phonon data from ASE database
        dict: Gruneisen parameters.

    """

    with connect(dbname) as db:
        entry = next(
            db.select(structure_name=structure_name, model_name=model_name), None
        )
        if entry is None:
            return None
        data = entry.data
        phonons = data["strain_phonons"]
        v0 = phonons["0.00"]["volume"]

        # Initialize dictionary to store Gruneisen parameters for all q-points
        gruneisen_all_q = {
            "q_points": phonons["0.00"]["bandstructure"]["q_points"],
            "band_index": {},
        }

        # Iterate over all q-points
        for q_index in range(len(gruneisen_all_q["q_points"])):
            omega0 = [
                b[q_index]
                for b in phonons["0.00"]["bandstructure"]["band_index"].values()
            ]

            gamma = np.expand_dims(-v0 / np.array(omega0), axis=1) * gruneisen(phonons)

            for band_index, g_value in enumerate(gamma, start=1):
                gruneisen_all_q["band_index"].setdefault(band_index, []).append(
                    g_value[q_index]
                )

        return phonons, gruneisen_all_q


def generate_mode_gruneisen_table(
    dbname, chemsys, output_filename, structure_name="B2", q="M"
):
    """
    Create a LaTeX table of the M point Gruneisen parameter.
    """
    db = connect(dbname)
    unique_structure_types = set(
        entry.structure_name for entry in db.select(chemsys=chemsys)
    )

    with open(output_filename, "w") as file:
        file.write("\\begin{table}[h]\n")
        file.write("\\centering\n")

        # Begin tabular environment
        file.write("\\begin{tabular}{|c|c|}\n")
        file.write("\\hline\n")
        file.write('Structure & \\textbf{M}-Mode Gr\\"{u}neisen Parameters \\\\\n')
        file.write("\\hline\n")

        for name in unique_structure_types:
            if name != structure_name:
                continue
            structure_name = name.replace("_", "-")
            file.write(f"{structure_name} & ")
            file.write("\\begin{tabular}{c c c c c c c}\n")  # Inner table header

            # NOTE: The labels are based on low to high eigenvalues which is not
            # Guarenteed to be correct! The proper way is to identify what the
            # cross product of the mode eigenvector  w.r.t the q-vector. This needs
            # to be done in the CalculatePhonons.py
            file.write("Model & TA$_1$ & TA$_2$ & LA & TO$_1$ & TO$_2$ & LO \\\\\n")
            file.write("\\hline\n")

            entries = db.select(chemsys=chemsys, structure_name=name)

            # Filter and sort based on models defined in ORDER
            sentries = sorted(
                [
                    entry
                    for entry in entries
                    if entry.model_name in ORDER[chemsys]["model"]
                ],
                key=lambda entry: ORDER[chemsys]["model"][entry.model_name],
            )

            for entry in sentries:
                model_name = entry.get("model_name", "NA")
                _, gamma = calculate_mode_gruneisen(entry)

                file.write(
                    f"{model_name} & {gamma[0]:.2f} & {gamma[1]:.2f} & {gamma[2]:.2f} & {gamma[3]:.2f} & {gamma[4]:.2f} & {gamma[5]:.2f} \\\\\n"
                )

            file.write("\\hline\n")
            file.write("\\end{tabular}\\\\\n")
            file.write("\\hline\n")
        file.write("\\end{tabular}\n")
        file.write(
            "\\caption{Mode Gruneisen parameters for different NiTi structures.}\n"
        )
        file.write("\\label{tab:mode_gruneisen_niti}\n")
        file.write("\\end{table}\n")


def generate_mode_gruneisen_plot(
    dbname, chemsys, output_filename, structure_name="B2", q="M"
):
    """
    Plot for given q-point the gruneisen param vs. eigenfrequency
    """
    db = connect(dbname)
    unique_models = db.select(chemsys=chemsys, structure_name=structure_name)

    # Filter and sort based on models defined in ORDER
    sentries = sorted(
        [
            entry
            for entry in unique_models
            if entry.model_name in ORDER[chemsys]["model"]
        ],
        key=lambda entry: ORDER[chemsys]["model"][entry.model_name],
    )
    # Figure size for one-column width (inches)
    plt.figure(figsize=(3.5, 2.5))

    for entry in sentries:
        model_name = entry.get("model_name", "NA")
        omega0, gammas = calculate_mode_gruneisen(entry)
        plt.plot(omega0, log_transform(gammas), label=f"{model_name}", marker='o', markersize=4, linestyle='-')

    plt.xlabel("Eigenfrequency [THz]")
    m = r'\text{sign}(\gamma) \cdot \ln(|\gamma|)$'
    plt.ylabel(r"Gruneisen Parameter [{q}-Mode,{m}]")
    plt.legend()
    plt.tight_layout()  # Adjusts the plot to fit into the figure area.
    plt.savefig(output_filename, format='png')  # Save as PDF for high quality
    plt.close()


def plot_gruneisen(gruneisen_data, output_filename):
    if gruneisen_data is None:
        print("No data to plot.")
        return

    q_points = gruneisen_data["q_points"]
    num_bands = len(gruneisen_data["band_index"])

    plt.figure()
    for band_index in range(1, num_bands + 1):
        gamma_values = gruneisen_data["band_index"][band_index]
        plt.plot(q_points, gamma_values, label=f"Band {band_index}")
    plt.ylim((-50, 50))
    plt.xlabel("Q-point")
    plt.ylabel("Gruneisen Parameter")
    plt.title("Mode Gruneisen Parameters for Each Band")
    plt.legend()
    plt.savefig(output_filename)


def plot_gruneisen_frequency(phonons, gruneisen_data, output_filename):
    if gruneisen_data is None:
        UserWarning("No data to plot.")

    q_points = gruneisen_data["q_points"]
    num_bands = len(gruneisen_data["band_index"])

    plt.figure(figsize=(10, 6))
    for band_index in range(1, num_bands + 1):
        gamma_values = gruneisen_data["band_index"][band_index]
        frequencies = phonons["0.00"]["bandstructure"]["band_index"][band_index - 1]

        # Scatter plot for each band
        plt.scatter(frequencies, gamma_values, label=f"Band {band_index}")

    plt.ylim((-50, 50))
    # plt.xlabel("Frequency")
    # plt.ylabel("Gruneisen Parameter")
    plt.title("Gruneisen Parameters vs. Frequency")
    plt.legend()
    plt.savefig(output_filename)


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    structure_name = sys.argv[3]
    q = sys.argv[4]
    output_file = paths.output / f"Table_{chemsys}_{q}_ModeGruneisen.tex"
    generate_mode_gruneisen_table(
        paths.data / dbname, chemsys, output_file, structure_name=structure_name, q=q
    )
    generate_mode_gruneisen_plot(
        paths.data / dbname, chemsys, paths.figures / f"Plot_{chemsys}_{q}_ModeGruneisen.png", structure_name=structure_name, q=q
    )
