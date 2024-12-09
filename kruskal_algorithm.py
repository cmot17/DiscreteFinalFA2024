from main import read_fasta_file, create_distance_edges, Graph

if __name__ == '__main__':
    # Read FASTA sequences from file
    fasta_data = read_fasta_file('output.fasta')

    # Print the sequences from file
    print("Sequences from FASTA file:")
    for header, seq in fasta_data:
        print(f"{header}\n{seq}")

    # Create a graph from the sequences
    g = Graph(len(fasta_data))
    edges = create_distance_edges(fasta_data)
    for u, v, w in edges:
        g.addEdge(u, v, w)

    # Apply Kruskal's algorithm to find MST
    kruskal_edges = g.kruskal_mst()
