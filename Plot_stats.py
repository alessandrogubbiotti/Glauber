import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os

def parse_header(header_line):
    headers = header_line.strip().split()
    stat_names = []

    # Group columns by statistic name
    sim_stat_labels = headers[1:]  # Skip "time"
    sim_ids = set()
    for label in sim_stat_labels:
        sim_part, stat = label.split('_', 1)
        stat_names.append(stat)
        sim_ids.add(sim_part)
    
    N_simulations = len(sim_ids)
    N_stats = len(stat_names) // N_simulations
    return stat_names[:N_stats], N_simulations

def plot_statistics(df, stat_names, N_simulations, time_col="time", outdir="plots"):
    os.makedirs(outdir, exist_ok=True)
    time = df[time_col]

    for stat in stat_names:
        plt.figure(figsize=(10, 5))
        values = []

        # Plot all simulation curves
        for sim in range(N_simulations):
            col = f"sim{sim}_{stat}"
            if col in df.columns:
                plt.plot(time, df[col], label=f"sim {sim}", alpha=0.5)
                values.append(df[col])
        
        # Plot the mean across simulations
        if values:
            mean_series = pd.concat(values, axis=1).mean(axis=1)
            plt.plot(time, mean_series, label="mean", color='black', linewidth=2)

        plt.title(f"{stat} over time")
        plt.xlabel("Time")
        plt.ylabel(stat)
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{stat}.png"))
        plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Path to the stats .txt file")
    parser.add_argument("--outdir", default="plots", help="Directory to save plots")
    args = parser.parse_args()

    # Read header separately to parse stats
    with open(args.file, "r") as f:
        header_line = f.readline()
    
    stat_names, N_simulations = parse_header(header_line)

    # Now load the full DataFrame (skipping comment prefix)
    df = pd.read_csv(args.file, sep=r"\s+", comment="#")

    plot_statistics(df, stat_names, N_simulations, outdir=args.outdir)

    print(f"✅ Plots saved in '{args.outdir}/'")

if __name__ == "__main__":
    main()

