
import matplotlib.pyplot as plt


def PlotGruneisen(gruneisen_data):
    if gruneisen_data is None:
        print("No data to plot.")
        return

    q_points = gruneisen_data['q_points']
    num_bands = len(gruneisen_data['band_index'])

    plt.figure(figsize=(10, 6))
    for band_index in range(1, num_bands + 1):
        gamma_values = gruneisen_data['band_index'][band_index]
        plt.plot(q_points, gamma_values, label=f'Band {band_index}')
    plt.ylim((-50,50))
    plt.xlabel('Q-point')
    plt.ylabel('Gruneisen Parameter')
    plt.title('Mode Gruneisen Parameters for Each Band')
    plt.legend()
    plt.grid(True)
    plt.savefig("/tmp/GruneisenAll.png")

def PlotGruneisenFrequency(ph,gruneisen_data):
    if gruneisen_data is None:
        print("No data to plot.")
        return

    q_points = gruneisen_data['q_points']
    num_bands = len(gruneisen_data['band_index'])
    
    plt.figure(figsize=(10, 6))
    for band_index in range(1, num_bands + 1):
        gamma_values = gruneisen_data['band_index'][band_index]
        frequencies = ph["0.00"]["bandstructure"]["band_index"][band_index-1]

        # Scatter plot for each band
        plt.scatter(frequencies, gamma_values, label=f'Band {band_index}')

    plt.ylim((-50, 50))
    plt.xlabel('Frequency')
    plt.ylabel('Gruneisen Parameter')
    plt.title('Gruneisen Parameters vs. Frequency')
    plt.legend()
    plt.grid(True)
    plt.savefig("/tmp/Gruneisen_vs_Frequency.png")

