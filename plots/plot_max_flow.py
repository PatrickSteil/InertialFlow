import io
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_benchmark_results(csv_path_or_data):
    # 1. Load the data
    # Can accept a file path string or an io.StringIO object
    if isinstance(csv_path_or_data, str) and not csv_path_or_data.startswith(
        "bench_name"
    ):
        df = pd.read_csv(csv_path_or_data)
    else:
        df = pd.read_csv(
            io.StringIO(csv_path_or_data)
            if isinstance(csv_path_or_data, str)
            else csv_path_or_data
        )

    # Strip any accidental whitespace from column names and string columns
    df.columns = df.columns.str.strip()
    if "algorithm" in df.columns:
        df["algorithm"] = df["algorithm"].str.strip()

    # Calculate Total Time for comprehensive analysis
    df["total_time"] = df["build_time"] + df["solve_time"]

    # 2. Set up the plotting style
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=False)

    # Sort algorithms by median solve time so the plot looks organized
    order = (
        df.groupby("algorithm")["solve_time"].median().sort_values().index
    )

    # Plot 1: Solve Time Performance
    sns.boxplot(
        ax=axes[0],
        x="algorithm",
        y="solve_time",
        data=df,
        order=order,
        palette="Blues_r",
        hue="algorithm",
        legend=False,
    )
    axes[0].set_title("Algorithm Solve Time Performance", fontsize=14, pad=15)
    axes[0].set_xlabel("Algorithm", fontsize=12)
    axes[0].set_ylabel("Solve Time (seconds)", fontsize=12)
    axes[0].tick_params(axis="x", rotation=45)

    # Plot 2: Total Time Performance (Build + Solve)
    sns.boxplot(
        ax=axes[1],
        x="algorithm",
        y="total_time",
        data=df,
        order=order,
        palette="Oranges_r",
        hue="algorithm",
        legend=False,
    )
    axes[1].set_title(
        "Total Time Performance (Build + Solve)", fontsize=14, pad=15
    )
    axes[1].set_xlabel("Algorithm", fontsize=12)
    axes[1].set_ylabel("Total Time (seconds)", fontsize=12)
    axes[1].tick_params(axis="x", rotation=45)

    # Adjust layout to make room for rotated labels
    plt.tight_layout()

    # Save and show the plot
    plt.savefig("max_flow_benchmark_results.png", dpi=300)
    print("Plot successfully saved to 'max_flow_benchmark_results.png'")
    plt.show()


# --- How to use it ---
if __name__ == "__main__":
    # Replace 'benchmark_output.csv' with the path to your actual log file
    csv_file_path = 'output.multiple.csv.0'
    plot_benchmark_results(csv_file_path)
