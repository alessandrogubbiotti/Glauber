import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import json

MARKER = -9999  # special int separating time steps

def load_config(path):
    with open(path, 'r') as f:
        return json.load(f)

def load_binary_spin_data(path, n_simulations, l, n_time_steps):
    with open(path, 'rb') as f:
        raw = np.fromfile(f, dtype=np.int32)

    expected_block = n_simulations * l
    expected_length = n_time_steps * (expected_block + 1)

    if len(raw) != expected_length:
        raise ValueError(f"Unexpected file size. Expected {expected_length} entries, got {len(raw)}")

    data = np.zeros((n_simulations, n_time_steps, l), dtype=np.int32)

    for t in range(n_time_steps):
        start = t * (expected_block + 1)
        end = start + expected_block
        block = raw[start:end]

        for sim in range(n_simulations):
            data[sim, t, :] = block[sim * l : (sim + 1) * l]

        # Optional: check the marker
        marker = raw[end]
        if marker != MARKER:
            raise ValueError(f"Missing marker at timestep {t}, found {marker} instead.")

    return data


def plot_simulation(sim_data, sim_index, outdir, T, L, N, n_micro_evolution,
                    creation_rate, annihilation_rate, resolution):
    params_text = (
        f"N = {N}\n"
        f"Time discretization of\nmicroscopic dynamics = {n_micro_evolution}\n"
        f"Creation rate = {creation_rate}\n"
        f"Annihilation rate = {annihilation_rate}\n"
        f"Resolution = {resolution}"
    )

    n_time_steps, l = sim_data.shape
    extent = (0, L, 0, T)

    plt.figure(figsize=(10, 4))
    plt.imshow(sim_data, cmap='bwr', aspect='auto', interpolation='nearest',
               vmin=-1, vmax=1, extent=extent, origin='lower')
    plt.xlabel("Macroscopic Space Variable")
    plt.ylabel("Macroscopic Time")
    plt.title(f"Simulation {sim_index}")
    plt.colorbar(label='Spin')

    plt.text(
        0.95, 0.95, params_text,
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='top',
        horizontalalignment='right',
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)
    )

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"simulation_{sim_index}.png"))
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", help="Path to the result directory")
    parser.add_argument("--spins_file", default="spins_output.bin", help="Name of the binary spins output file")
    parser.add_argument("--config_file", default="configuration.json", help="Name of the config file")
    args = parser.parse_args()

    config_path = os.path.join(args.results_dir, args.config_file)
    spins_path = os.path.join(args.results_dir, args.spins_file)

    config = load_config(config_path)
    T = config["T"]
    n_simulations = config["N_simulations"]
    N = config["N"]
    L = config["L"]
    Micro_n_steps = config["Micro_n_steps"]
    creation = config["Micro_Creation_Rate"]
    annihilation = config["Micro_Annihilation_Rate"]
    resolution = config["resolution"]

    l = N * L
    total_time_steps = int(resolution * T)

    spins = load_binary_spin_data(spins_path, n_simulations, l, total_time_steps)
    print("spins shape:", spins.shape)

    for i in range(n_simulations):
        plot_simulation(spins[i], i, args.results_dir, T, L, N,  Micro_n_steps, creation,  annihilation, resolution )

    print(f"✅ Plotted {n_simulations} simulations from binary file. Output in: {args.results_dir}")

if __name__ == "__main__":
    main()

