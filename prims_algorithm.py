from main import read_fasta_file, create_distance_edges, prim_mst

if __name__ == '__main__':
    # Read FASTA sequences from file
    fasta_data = read_fasta_file('output.fasta')

    # Print the sequences from file
    print("Sequences from FASTA file:")
    for header, seq in fasta_data:
        print(f"{header}\n{seq}")

    # Create graph edges using Hamming distance
    edges = create_distance_edges(fasta_data)

    # Apply Prim's algorithm
    mst_weight, mst_edges = prim_mst(len(fasta_data), edges)

    # Display the results
    print(f"\nTotal weight of the MST (using Prim's algorithm): {mst_weight}")
    print("Edges in the MST:")
    for u, v, weight in mst_edges:
        print(f"Sequence {u+1} -- Sequence {v+1} [Weight: {weight}]")
