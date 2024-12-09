import random
import heapq


# FASTA generator function
def generate_fasta_sequences(num_sequences, max_length):
    """
    Generates random DNA sequences.

    Parameters:
        num_sequences (int): Number of sequences to generate.
        max_length (int): Maximum length of the sequences.

    Returns:
        list: A list of tuples with headers and sequences.
    """
    nucleotides = 'ACGT'
    fasta_entries = []

    for i in range(1, num_sequences + 1):
        seq_length = random.randint(1, max_length)
        sequence = ''.join(random.choices(nucleotides, k=seq_length))
        header = f">sequence_{i}"
        fasta_entries.append((header, sequence))

    return fasta_entries


# Function to compute Hamming distance between two sequences
def hamming_distance(seq1, seq2):
    """
    Calculates the Hamming distance between two sequences.
    If sequences are of different lengths, compares up to the shortest length.

    Parameters:
        seq1 (str): First DNA sequence.
        seq2 (str): Second DNA sequence.

    Returns:
        int: Hamming distance between the sequences.
    """
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


# Prim's algorithm for Minimum Spanning Tree
def prim_mst(V, edges):
    """
    Implements Prim's algorithm to find the MST.

    Parameters:
        V (int): Number of vertices.
        edges (list): List of edges, where each edge is (u, v, weight).

    Returns:
        tuple: Total weight of the MST and list of edges in the MST.
    """
    # Create an adjacency list
    adj = [[] for _ in range(V)]
    for u, v, wt in edges:
        adj[u].append((v, wt))
        adj[v].append((u, wt))

    # Priority queue for Prim's algorithm
    pq = []
    visited = [False] * V
    mst_edges = []
    total_weight = 0

    # Start with vertex 0
    heapq.heappush(pq, (0, 0, -1))  # (weight, current_vertex, parent_vertex)

    while pq:
        wt, u, parent = heapq.heappop(pq)
        if visited[u]:
            continue  # Skip if already visited
        total_weight += wt
        visited[u] = True

        # Add the edge to the MST if it has a valid parent
        if parent != -1:
            mst_edges.append((parent, u, wt))

        for v, weight in adj[u]:
            if not visited[v]:
                heapq.heappush(pq, (weight, v, u))

    return total_weight, mst_edges


# Driver Code
if __name__ == '__main__':
    # Generate FASTA sequences
    num_sequences = 5
    max_length = 10
    fasta_data = generate_fasta_sequences(num_sequences, max_length)

    # Print the generated sequences
    print("Generated FASTA Sequences:")
    for header, seq in fasta_data:
        print(f"{header}\n{seq}")

    # Create graph edges using Hamming distance
    edges = []
    for i in range(len(fasta_data)):
        for j in range(i + 1, len(fasta_data)):
            seq1 = fasta_data[i][1]
            seq2 = fasta_data[j][1]
            distance = hamming_distance(seq1, seq2)
            edges.append((i, j, distance))

    # Apply Prim's algorithm
    mst_weight, mst_edges = prim_mst(len(fasta_data), edges)

    # Display the results
    print(f"\nTotal weight of the MST (using Prim's algorithm): {mst_weight}")
    print("Edges in the MST:")
    for u, v, weight in mst_edges:
        print(f"{u} -- {v} [Weight: {weight}]")
