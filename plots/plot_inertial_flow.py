import io
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def load_and_clean_log(filepath):
    """Reads the log file, dropping lines that are incomplete or malformed."""
    valid_lines = []

    with open(filepath, "r") as f:
        for line in f:
            # Check if it looks like a header or a valid data row
            if "level" in line or (line.strip() and line.count(",") == 6):
                valid_lines.append(line.strip())

    # Load into data frame
    df = pd.read_csv(io.StringIO("\n".join(valid_lines)))

    # Ensure correct data types
    numeric_cols = [
        "level",
        "lo",
        "hi",
        "num_vertices",
        "best_projection_deg",
        "best_flow",
        "time_ms",
    ]
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors="coerce")
    return df.dropna()


def plot_partition_metrics(df):
    # Set style for cleaner aesthetics
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle(
        "Inertial Flow Algorithm Diagnostics (Germany Dataset)",
        fontsize=18,
        fontweight="bold",
    )

    # 1. Total Time spent per Level
    time_per_level = df.groupby("level")["time_ms"].sum().reset_index()
    sns.barplot(
        ax=axes[0, 0],
        data=time_per_level,
        x="level",
        y="time_ms",
        hue="level",
        palette="viridis",
        legend=False,
    )
    axes[0, 0].set_title("Total Runtime (ms) Spent per Hierarchical Level")
    axes[0, 0].set_xlabel("Tree Level (Depth)")
    axes[0, 0].set_ylabel("Total Time (ms)")

    # 2. Flow Value vs. Number of Vertices
    sns.scatterplot(
        ax=axes[0, 1],
        data=df,
        x="num_vertices",
        y="best_flow",
        hue="level",
        palette="magma",
        alpha=0.7,
        edgecolor=None,
    )
    axes[0, 1].set_title("Best Flow (Cut Size) vs. Problem Size")
    axes[0, 1].set_xlabel("Number of Vertices")
    axes[0, 1].set_ylabel("Best Flow Value")
    axes[0, 1].set_xscale("log")
    axes[0, 1].set_yscale("log")

    # 3. Workload distribution per Level (Boxplot of vertices)
    sns.boxplot(
        ax=axes[1, 0],
        data=df,
        x="level",
        y="num_vertices",
        hue="level",
        palette="cubehelix",
        legend=False,
    )
    axes[1, 0].set_title("Distribution of Subgraph Sizes per Level")
    axes[1, 0].set_xlabel("Tree Level (Depth)")
    axes[1, 0].set_ylabel("Vertices per Subproblem")
    axes[1, 0].set_yscale("log")

    # 4. Runtime Scaling: Time vs Size
    sns.scatterplot(
        ax=axes[1, 1],
        data=df,
        x="num_vertices",
        y="time_ms",
        hue="best_projection_deg",
        palette="crest",
        alpha=0.8,
    )
    axes[1, 1].set_title("Execution Time Scaling per Bisection Step")
    axes[1, 1].set_xlabel("Number of Vertices")
    axes[1, 1].set_ylabel("Step Runtime (ms)")
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")
    axes[1, 1].get_legend().set_title("Projection Deg")

    plt.tight_layout()

    # Save visualization output
    output_filename = "inertial_flow_diagnostics.png"
    plt.savefig(output_filename, dpi=300)
    print(f"[✓] Dashboard visualization saved successfully to {output_filename}")
    plt.show()


if __name__ == "__main__":
    # Point this to your log file path
    log_file_path = "../germany.1024.log"

    try:
        data = load_and_clean_log(log_file_path)
        print(f"Loaded {len(data)} successful bisection steps from log.")
        plot_partition_metrics(data)
    except FileNotFoundError:
        print(
            f"Error: Could not find '{log_file_path}'. Make sure the script is running in the same directory as your log file."
        )
