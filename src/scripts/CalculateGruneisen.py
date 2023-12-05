import paths
import numpy as np
import numpy.polynomial.polynomial as poly
from ase.db import connect


def FindIndex(data, label):
    """
        locates the index for the q-point value that is closest to the
    point of interest.
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


def ModeGruneisen(phonons, q, n=3):
    """
    Calculate the specific mode Gruneisen parameter via a polynomial fit.
    """

    volumes = []
    ef_list = []
    for k in phonons.keys():
        volumes.append(phonons[k]["volume"])
        mode_index = FindIndex(phonons[k], q)
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


def Gruneisen(phonons, n=3):
    """
    Calculate the Gruneisen parameters via a polynomial fit.
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


# def CalculateModeGruneisen(dbname,model_name,structure_name,q="M"):
#    with connect(dbname) as db:
#        entry = next(db.select(structure_name=structure_name, model_name=model_name),None)
#        data = entry.data
#        phonons = data["strain_phonons"]
#        v0 = phonons["0.00"]["volume"]
#        mode_index = FindIndex(phonons["0.00"],q)
#        omega0 = [ b[mode_index] for b in phonons["0.00"]["bandstructure"]["band_index"].values()]
#        gamma = -v0/np.array(omega0) * ModeGruneisen(phonons,q)
#    return phonons,gamma


def CalculateModeGruneisen(dbentry, q="M"):
    data = dbentry.data
    phonons = data["strain_phonons"]
    v0 = phonons["0.00"]["volume"]
    mode_index = FindIndex(phonons["0.00"], q)
    omega0 = [
        b[mode_index] for b in phonons["0.00"]["bandstructure"]["band_index"].values()
    ]
    gamma = -v0 / np.array(omega0) * ModeGruneisen(phonons, q)
    return phonons, gamma


def CalculateGruneisen(dbname, model_name, structure_name):
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
            gamma = np.expand_dims(-v0 / np.array(omega0), axis=1) * Gruneisen(phonons)
            for band_index, g_value in enumerate(gamma, start=1):
                gruneisen_all_q["band_index"].setdefault(band_index, []).append(
                    g_value[q_index]
                )

        return phonons, gruneisen_all_q
