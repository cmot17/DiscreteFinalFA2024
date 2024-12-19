"""
Main module for processing genomic sequences and generating minimum spanning trees.

This module implements a pipeline for processing genomic sequences using k-mer analysis
and generating minimum spanning trees using both Kruskal's and Prim's algorithms.
"""

import os
from Bio import SeqIO
import pandas as pd
import numpy as np
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import subprocess
import tempfile
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from kruskal_algorithm import Graph
from prims_algorithm import prim_mst

# Configuration constants
KMER_LENGTH = 31
BASE_INPUT_DIR = "unassembled-plants"
BASE_OUTPUT_DIR = "output"
JELLYFISH_MEM = "16G"  # Memory limit for Jellyfish
JELLYFISH_THREADS = max(1, mp.cpu_count() - 1)  # Use all but one CPU core

# Output directories (now top-level)
DISTANCE_DIR = "distance"  # For PHYLIP distance matrices
CSV_DIR = "csvs"  # For CSV files
VIZ_DIR = "viz"  # For visualizations

# Species color mapping for visualization
SPECIES_COLORS = {
    "rapa": "#FF9999",  # Light red
    "sinensis": "#99FF99",  # Light green
    "camaldule": "#9999FF",  # Light blue
    "cacao": "#FFFF99",  # Light yellow
    "rubella": "#FF99FF",  # Light magenta
    "vinifera": "#99FFFF",  # Light cyan
    "thalian": "#FFB366",  # Light orange
    "raimondii": "#B366FF",  # Light purple
    "papaya": "#66FFB3",  # Light mint
    "parvulum": "#66B3FF",  # Light sky blue
    "halophilu": "#FFB3B3",  # Light salmon
    "clementin": "#B3FFB3",  # Light lime
    "grandis": "#B3B3FF",  # Light periwinkle
    "lyrata": "#FFE6B3",  # Light peach
}


def ensure_output_dirs() -> dict[str, str]:
    """
    Create top-level output directories if they don't exist.

    Returns:
        Dictionary containing paths to output directories
    """
    dirs = {"distance": DISTANCE_DIR, "csv": CSV_DIR, "viz": VIZ_DIR}

    for dir_path in dirs.values():
        os.makedirs(dir_path, exist_ok=True)

    return dirs


def find_fasta_files(input_dir: str) -> list[str]:
    """
    Find all FASTA/FASTQ files in the input directory.

    Args:
        input_dir: Directory to search for files

    Returns:
        List of paths to FASTA/FASTQ files
    """
    valid_extensions = (".fasta", ".fa", ".fna", ".fastq", ".fq")
    fasta_files = []
    for root, _, files in os.walk(input_dir):
        for f in files:
            if f.lower().endswith(valid_extensions):
                fasta_files.append(os.path.join(root, f))
    return fasta_files


def get_node_color(name: str) -> str:
    """
    Get color for node based on species name.

    Args:
        name: Name of the species

    Returns:
        Hex color code for the species
    """
    for species, color in SPECIES_COLORS.items():
        if species in name.lower():
            return color
    return "#CCCCCC"  # Default gray


def write_phylip(
    distance_matrix: pd.DataFrame, samples: list[str], output_file: str
) -> None:
    """
    Write the distance matrix in PHYLIP format.

    Args:
        distance_matrix: DataFrame with samples as both rows and columns
        samples: List of sample names
        output_file: Path to output file
    """
    n = len(samples)
    with open(output_file, "w") as f:
        f.write(f"    {n}\n")
        for s in samples:
            name = s[:10].ljust(10)  # PHYLIP format requires 10-char names
            row = distance_matrix.loc[s, :].values
            row_str = " ".join(f"{dist:8.6f}" for dist in row)
            f.write(f"{name} {row_str}\n")


def create_graph_viz(
    mst_edges: list, sample_names: list, algorithm_name: str, coverage_name: str
) -> None:
    """
    Create and save MST visualization.

    Args:
        mst_edges: List of edges in the MST
        sample_names: List of sample names
        algorithm_name: Name of the algorithm used
        coverage_name: Name of the coverage folder
    """
    output_file = os.path.join(
        VIZ_DIR, f"mst_{algorithm_name.lower()}_{coverage_name}.png"
    )
    if os.path.exists(output_file):
        print(f"Visualization for {algorithm_name} already exists, skipping...")
        return

    G = nx.Graph()

    # Add nodes with labels
    for i, name in enumerate(sample_names):
        G.add_node(i, label=name)

    # Add edges with weights
    edge_weights = []
    for u, v, w in mst_edges:
        G.add_edge(u, v, weight=w)
        edge_weights.append(w)

    plt.figure(figsize=(20, 15))

    # Use Kamada-Kawai layout for better node distribution
    pos = nx.kamada_kawai_layout(G, weight="weight")

    # Calculate edge widths based on weights (inverse relationship)
    max_weight = max(edge_weights)
    min_weight = min(edge_weights)
    edge_widths = [
        float(3 * (1 - (w - min_weight) / (max_weight - min_weight)) + 1)
        for w in edge_weights
    ]

    # Draw edges
    nx.draw_networkx_edges(
        G, pos, width=1.0, edge_color="#666666", alpha=0.6, arrows=False
    )

    # Draw nodes
    node_colors = [get_node_color(sample_names[n]) for n in G.nodes()]
    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=2000,
        node_color=node_colors,  # type: ignore[arg-type]
        edgecolors="#333333",
        linewidths=2,
    )

    # Add node labels with white background
    labels = {i: name[:15] for i, name in enumerate(sample_names)}
    for node, (x, y) in pos.items():
        plt.text(
            x,
            y,
            labels[node],
            fontsize=10,
            ha="center",
            va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2),
        )

    # Add edge labels
    edge_labels = nx.get_edge_attributes(G, "weight")
    edge_labels = {k: f"{v:.4f}" for k, v in edge_labels.items()}
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels,
        font_size=8,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=2),
    )

    # Add title
    plt.title(
        f"Minimum Spanning Tree - {algorithm_name}\n{coverage_name}\n"
        f"K-mer Length: {KMER_LENGTH}",
        fontsize=16,
        pad=20,
    )

    # Add legend for species colors
    legend_elements = []
    seen_species = set()
    for name in sample_names:
        species = next(
            (s for s in name.lower().split("_") if s in SPECIES_COLORS), None
        )
        if species and species not in seen_species:
            seen_species.add(species)
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=get_node_color(species),
                    markersize=15,
                    label=species.capitalize(),
                )
            )

    if legend_elements:
        plt.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(1, 0.5))

    plt.axis("off")
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches="tight", pad_inches=0.5)
    plt.close()


def count_kmers_with_jellyfish(
    files: list[str], k: int, tmpdir: str
) -> dict[str, dict[str, int]]:
    """
    Count k-mers using Jellyfish with concurrent file processing.

    Args:
        files: List of input FASTA/FASTQ files
        k: Length of k-mers to count
        tmpdir: Temporary directory path

    Returns:
        Dictionary mapping filenames to their k-mer count dictionaries
    """
    kmer_counts_by_file = {}

    for file_path in files:
        # Create generator file for this file
        generator_file = os.path.join(tmpdir, "generators")
        with open(generator_file, "w") as f:
            f.write(f"cat {file_path}\n")

        jf_db = os.path.join(tmpdir, f"kmers_{os.path.basename(file_path)}.jf")

        # Run Jellyfish count with generators
        count_cmd = [
            "jellyfish",
            "count",
            "-m",
            str(k),
            "-s",
            "100M",  # Initial hash size
            "-t",
            str(JELLYFISH_THREADS),
            "-C",  # Count canonical k-mers
            "-g",
            generator_file,  # Use generator file
            "-G",
            "4",  # Run 4 generators concurrently
            "-o",
            jf_db,
        ]
        with open(os.devnull, "w") as devnull:
            subprocess.run(count_cmd, check=True, stderr=devnull)

        # Dump the counts
        dump_file = os.path.join(tmpdir, f"kmers_{os.path.basename(file_path)}.txt")
        dump_cmd = ["jellyfish", "dump", "-c", jf_db]
        with open(dump_file, "w") as f:
            subprocess.run(dump_cmd, check=True, stdout=f)

        # Read the counts into a dictionary
        kmer_counts = {}
        total_kmers = 0
        unique_kmers = 0

        with open(dump_file) as f:
            for line in f:
                kmer, count = line.strip().split()
                count = int(count)
                kmer_counts[kmer] = count
                total_kmers += count
                unique_kmers += 1

        print(f"\nStats for {os.path.basename(file_path)}:")
        print(f"- Total k-mers: {total_kmers:,}")
        print(f"- Unique k-mers: {unique_kmers:,}")

        kmer_counts_by_file[file_path] = kmer_counts

    return kmer_counts_by_file


def compute_weighted_jaccard(freqA: dict[str, int], freqB: dict[str, int]) -> float:
    """
    Compute weighted Jaccard distance between two k-mer frequency dictionaries.

    The function first normalizes the k-mer counts by dividing by total counts,
    then computes the Jaccard similarity.

    Args:
        freqA: First k-mer frequency dictionary
        freqB: Second k-mer frequency dictionary

    Returns:
        Weighted Jaccard distance (1 - similarity)
    """
    all_kmers = set(freqA.keys()) | set(freqB.keys())
    if not all_kmers:
        return 1.0  # No shared k-mers

    # Normalize frequencies
    total_A = sum(freqA.values())
    total_B = sum(freqB.values())

    norm_A = {k: count / total_A for k, count in freqA.items()}
    norm_B = {k: count / total_B for k, count in freqB.items()}

    intersection_sum = sum(
        min(norm_A.get(kmer, 0), norm_B.get(kmer, 0)) for kmer in all_kmers
    )
    union_sum = sum(max(norm_A.get(kmer, 0), norm_B.get(kmer, 0)) for kmer in all_kmers)

    if union_sum == 0:
        return 1.0

    jaccard = intersection_sum / union_sum
    return 1 - jaccard


def create_mst_visualizations(dist_df: pd.DataFrame, coverage_name: str) -> None:
    """
    Create and visualize MSTs using both Kruskal's and Prim's algorithms.

    Args:
        dist_df: DataFrame containing the distance matrix
        coverage_name: Name of the coverage folder for labeling
    """
    # Ensure output directories exist
    ensure_output_dirs()

    # Check if MST files already exist
    kruskal_png = os.path.join(VIZ_DIR, f"mst_kruskals_algorithm_{coverage_name}.png")
    prim_png = os.path.join(VIZ_DIR, f"mst_prims_algorithm_{coverage_name}.png")

    if os.path.exists(kruskal_png) and os.path.exists(prim_png):
        print(f"MST visualizations already exist for {coverage_name}, skipping...")
        return

    n_vertices = len(dist_df)
    sample_names = list(dist_df.index)

    # Create edges list from distance matrix
    edges = [
        [i, j, dist_df.iloc[i, j]]
        for i in range(n_vertices)
        for j in range(i + 1, n_vertices)
    ]

    # Create MST using Kruskal's algorithm
    g_kruskal = Graph(n_vertices)
    for u, v, w in edges:
        g_kruskal.add_edge(u, v, w)
    kruskal_edges = g_kruskal.kruskal_mst()

    # Create MST using Prim's algorithm
    prim_weight, prim_edges = prim_mst(n_vertices, edges)

    # Create visualizations for both algorithms
    create_graph_viz(kruskal_edges, sample_names, "Kruskal's Algorithm", coverage_name)
    create_graph_viz(prim_edges, sample_names, "Prim's Algorithm", coverage_name)

    # Save edge lists to CSV
    for algorithm, edges in [("kruskal", kruskal_edges), ("prim", prim_edges)]:
        csv_path = os.path.join(CSV_DIR, f"mst_{algorithm}_{coverage_name}.csv")
        if not os.path.exists(csv_path):
            edge_df = pd.DataFrame(edges, columns=["Node1", "Node2", "Weight"])
            edge_df["Node1"] = edge_df["Node1"].map(lambda x: sample_names[x])
            edge_df["Node2"] = edge_df["Node2"].map(lambda x: sample_names[x])
            edge_df.to_csv(csv_path, index=False)


def process_coverage_folder(coverage_folder: str) -> str | None:
    """
    Process a single coverage folder and generate its distance matrix and MSTs.
    Uses Jellyfish's concurrent file processing.

    Args:
        coverage_folder: Name of the coverage folder to process

    Returns:
        Name of the processed folder if successful, None otherwise
    """
    input_dir = os.path.join(BASE_INPUT_DIR, coverage_folder)
    ensure_output_dirs()

    print(f"\nProcessing folder: {coverage_folder}")
    print(f"Input directory: {input_dir}")
    print(f"Using k-mer length: {KMER_LENGTH}")
    print(f"Using Jellyfish with {JELLYFISH_THREADS} threads")

    fasta_files = find_fasta_files(input_dir)
    if not fasta_files:
        print(f"No FASTA/FASTQ files found in {input_dir}")
        return None

    print(f"Found {len(fasta_files)} FASTA files")

    # Process all files at once using a temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Get sample names
        sample_names = [
            os.path.splitext(os.path.relpath(f, input_dir))[0] for f in fasta_files
        ]

        # Count k-mers for all files concurrently
        kmer_counts_by_file = count_kmers_with_jellyfish(
            fasta_files, KMER_LENGTH, tmpdir
        )

        # Sort results to ensure consistent ordering
        sample_names = sorted(sample_names)
        fasta_files = sorted(fasta_files)

        # Compute pairwise distances
        n = len(sample_names)
        dist_matrix = np.zeros((n, n))
        total_comparisons = (n * (n - 1)) // 2

        with tqdm(
            total=total_comparisons, desc=f"Computing distances for {coverage_folder}"
        ) as pbar:
            for i in range(n):
                for j in range(i + 1, n):
                    dist = compute_weighted_jaccard(
                        kmer_counts_by_file[fasta_files[i]],
                        kmer_counts_by_file[fasta_files[j]],
                    )

                    # Log high distances
                    if dist > 0.5:
                        shared = len(
                            set(kmer_counts_by_file[fasta_files[i]].keys())
                            & set(kmer_counts_by_file[fasta_files[j]].keys())
                        )
                        total = len(
                            set(kmer_counts_by_file[fasta_files[i]].keys())
                            | set(kmer_counts_by_file[fasta_files[j]].keys())
                        )
                        print(
                            f"\nHigh distance ({dist:.6f}) between {sample_names[i]} and {sample_names[j]}:"
                        )
                        print(f"- Shared k-mers: {shared:,}")
                        print(f"- Total unique k-mers: {total:,}")
                        print(f"- Raw Jaccard (unweighted): {shared/total:.6f}")

                    dist_matrix[i, j] = dist
                    dist_matrix[j, i] = dist
                    pbar.update(1)

    # Create DataFrame and save results
    dist_df = pd.DataFrame(dist_matrix, index=sample_names, columns=sample_names)

    phylip_path = os.path.join(DISTANCE_DIR, f"{coverage_folder}_k{KMER_LENGTH}.phy")
    csv_path = os.path.join(CSV_DIR, f"{coverage_folder}_k{KMER_LENGTH}.csv")

    write_phylip(dist_df, sample_names, phylip_path)
    dist_df.to_csv(csv_path)

    # Create and save MST visualizations
    create_mst_visualizations(dist_df, coverage_folder)

    print(f"\nResults for {coverage_folder} saved to:")
    print(f"- PHYLIP format: {phylip_path}")
    print(f"- CSV format: {csv_path}")
    print(f"- MST visualizations: {VIZ_DIR}/mst_*.png")

    return coverage_folder


def main() -> None:
    """
    Main function to process all coverage folders and generate MSTs.

    The function performs the following steps:
    1. Creates the output directory structure
    2. Finds all coverage folders in the input directory
    3. For each coverage folder:
       - Processes FASTA/FASTQ files using Jellyfish's multithreading
       - Computes k-mer frequencies
       - Calculates distance matrices
       - Generates MST visualizations
    """
    # Create output directories
    ensure_output_dirs()

    # Find all coverage folders
    coverage_folders = [
        item
        for item in os.listdir(BASE_INPUT_DIR)
        if os.path.isdir(os.path.join(BASE_INPUT_DIR, item))
        and item.startswith("coverage_")
    ]

    if not coverage_folders:
        print("No coverage folders found!")
        return

    print(f"Found {len(coverage_folders)} coverage folders to process")

    # Process each coverage folder
    for folder in tqdm(coverage_folders, desc="Processing coverage folders"):
        csv_path = os.path.join(CSV_DIR, f"{folder}_k{KMER_LENGTH}.csv")

        if os.path.exists(csv_path):
            print(
                f"\nFound existing distance matrix for {folder}, using it to generate MSTs..."
            )
            dist_df = pd.read_csv(csv_path, index_col=0)
            create_mst_visualizations(dist_df, folder)
            print(f"MST visualizations created for {folder}")
        else:
            print(
                f"\nNo existing distance matrix found for {folder}, computing from scratch..."
            )
            result = process_coverage_folder(folder)
            if result:
                print(f"Successfully processed: {result}")

    # Print summary
    print("\nProcessing complete!")
    print("Results saved in:")
    print(f"- Distance matrices: {DISTANCE_DIR}/")
    print(f"- CSV files: {CSV_DIR}/")
    print(f"- Visualizations: {VIZ_DIR}/")

    # List files in each directory
    for dir_name, dir_path in [
        ("Distance matrices", DISTANCE_DIR),
        ("CSV files", CSV_DIR),
        ("Visualizations", VIZ_DIR),
    ]:
        if os.path.exists(dir_path):
            files = os.listdir(dir_path)
            if files:
                print(f"\n{dir_name}:")
                for f in sorted(files):
                    print(f"  - {f}")


if __name__ == "__main__":
    main()
