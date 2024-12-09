import heapq
from utils import read_fasta_file, create_distance_edges

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
