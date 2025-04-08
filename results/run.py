
def run_analysis(csv_path1, csv_path2, out_path):
    import sys
    import os
    sys.path.append(os.path.abspath('..'))

    from gwpy.timeseries import TimeSeries
    import h5py
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd 
    from matplotlib.animation import FuncAnimation
    from IPython.display import HTML
    from lib.fromACC2 import MassPoint, GWStrainCalculator
    from matplotlib.patches import Circle

    data2 = pd.read_csv(csv_path1, header=None, names=["Time", "mass1_x", "mass1_y"])
    data3 = pd.read_csv(csv_path2, header=None, names=["Time", "mass2_x", "mass2_y"])
    data = pd.merge(data2, data3, on="Time", how="left")

    trajectory = []
    for _, row in data.iterrows():
        mass1_point = MassPoint(mass=25, x=row['mass1_x'], y=row['mass1_y'])
        mass2_point = MassPoint(mass=30, x=row['mass2_x'], y=row['mass2_y'])
        trajectory.append([mass1_point, mass2_point])

    calculator = GWStrainCalculator()
    r_observer = 1e20
    theta = 0
    phi = 0
    dt = data['Time'].iloc[1] - data['Time'].iloc[0]

    H_plus = []
    H_cross = []
    for i in range(1, len(trajectory) - 1):
        masses_t1 = trajectory[i - 1]
        masses_t2 = trajectory[i]
        masses_t3 = trajectory[i + 1]
        h_plus, h_cross = calculator.calculate_strain_components(
            masses_t1, masses_t2, masses_t3, dt, r_observer, theta, phi
        )
        H_plus.append(h_plus)
        H_cross.append(h_cross)

    fig, ax = plt.subplots(figsize=(12, 12))
    fig.patch.set_facecolor('lightgrey')
    ax.set_facecolor('lightgrey')
    ax.axis('off')
    ax.set_aspect('equal', adjustable='datalim')

    time = data['Time'][:-2]
    margin_factor = 1.1
    x_min, x_max = min(data['mass1_x'].min(), data['mass2_x'].min()), max(data['mass1_x'].max(), data['mass2_x'].max())
    y_min, y_max = min(data['mass1_y'].min(), data['mass2_y'].min()), max(data['mass1_y'].max(), data['mass2_y'].max())
    x_margin = (x_max - x_min) * (margin_factor - 1)
    y_margin = (y_max - y_min) * (margin_factor - 1)

    back = plt.Rectangle((x_min - x_margin, y_min - y_margin), 
                         (x_max - x_min) + 2 * x_margin, 
                         (y_max - y_min) + 2 * y_margin, 
                         color="lightgrey")

    star_1 = Circle((data['mass1_x'][0], data['mass1_y'][0]), 5, color='blue')
    star_2 = Circle((data['mass2_x'][0], data['mass2_y'][0]), 5, color='red')
    trajectory1, = ax.plot([], [], linestyle='dotted', color='blue', linewidth=1)
    trajectory2, = ax.plot([], [], linestyle='dotted', color='red', linewidth=1)

    ax.add_patch(back)
    ax.add_patch(star_1)
    ax.add_patch(star_2)
    trajectory1.set_data([], [])
    trajectory2.set_data([], [])

    inset_ax = fig.add_axes([0.70, 0.025, 0.25, 0.25], facecolor='white')
    inset_ax.set_xlim(0, max(time))
    inset_ax.set_ylim(
        min(np.min(H_plus), np.min(H_cross)) * 1.1, 
        max(np.max(H_plus), np.max(H_cross)) * 1.1
    )
    inset_ax.tick_params(colors='black')
    inset_ax.spines[:].set_color('black')
    inset_ax.xaxis.label.set_color('black')
    inset_ax.yaxis.label.set_color('black')
    inset_ax.set_xlabel('Time (s)')
    inset_ax.set_ylabel('Strain Amplitude')
    inset_ax.grid(True, color='black', linestyle=':', linewidth=0.5)

    waveform_plus, = inset_ax.plot([], [], linestyle='-', color='blue', label='$h_+$')
    waveform_cross, = inset_ax.plot([], [], linestyle='-', color='red', label='$h_\times$')
    inset_ax.legend(loc='upper right', facecolor='white', edgecolor='black')

    def update(frame):
        star_1.center = (data['mass1_x'][frame], data['mass1_y'][frame])
        star_2.center = (data['mass2_x'][frame], data['mass2_y'][frame])
        trajectory1.set_data(data['mass1_x'][:frame+1], data['mass1_y'][:frame+1])
        trajectory2.set_data(data['mass2_x'][:frame+1], data['mass2_y'][:frame+1])
        time_data = time[:frame+1]
        waveform_plus.set_data(time_data, H_plus[:frame+1])
        waveform_cross.set_data(time_data, H_cross[:frame+1])
        return star_1, star_2, trajectory1, trajectory2, waveform_plus, waveform_cross

    animation = FuncAnimation(fig, update, frames=len(time), interval=30, blit=True, save_count=0)
    plt.close()
    animation.save(out_path, writer='pillow', fps=20)

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Run animation analysis on two motion capture CSV files.")
    parser.add_argument('--csv1', type=str, required=True, help='Path to input CSV file for object 1')
    parser.add_argument('--csv2', type=str, required=True, help='Path to input CSV file for object 2')
    parser.add_argument('--out', type=str, default=None, help='Filename for the output animation (default auto-generated)')

    args = parser.parse_args()

    if args.out is None:
        base1 = os.path.splitext(os.path.basename(args.csv1))[0]
        base2 = os.path.splitext(os.path.basename(args.csv2))[0]
        out_path = f"{base1}_{base2}_animation.gif"
    else:
        out_path = args.out

    run_analysis(args.csv1, args.csv2, out_path)
