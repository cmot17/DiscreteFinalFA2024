import networkx as nx
import matplotlib.pyplot as plt
import heapq

def read_fasta_file(filename):
    """
    Reads sequences from a FASTA file.

    Parameters:
        filename (str): Path to the FASTA file.

    Returns:
        list: A list of tuples with headers and sequences.
    """
    fasta_entries = []
    current_header = None
    current_sequence = []

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    fasta_entries.append((current_header, ''.join(current_sequence)))
                current_header = line
                current_sequence = []
            else:
                current_sequence.append(line)
        
        if current_header:  # Add the last sequence
            fasta_entries.append((current_header, ''.join(current_sequence)))

    return fasta_entries

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

def create_distance_edges(fasta_data):
    """
    Creates edges between sequences based on Hamming distance.

    Parameters:
        fasta_data (list): List of tuples containing (header, sequence).

    Returns:
        list: List of edges where each edge is (u, v, distance).
    """
    edges = []
    for i in range(len(fasta_data)):
        for j in range(i + 1, len(fasta_data)):
            seq1 = fasta_data[i][1]
            seq2 = fasta_data[j][1]
            distance = hamming_distance(seq1, seq2)
            edges.append((i, j, distance))
    return edges

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

class Graph:
    def __init__(self, vertices):
        self.V = vertices
        self.graph = []

    def addEdge(self, u, v, w):
        self.graph.append([u, v, w])

    def find(self, parent, i):
        if parent[i] != i:
            parent[i] = self.find(parent, parent[i])
        return parent[i]

    def union(self, parent, rank, x, y):
        if rank[x] < rank[y]:
            parent[x] = y
        elif rank[x] > rank[y]:
            parent[y] = x
        else:
            parent[y] = x
            rank[x] += 1

    def kruskal_mst(self):
        result = []
        i = 0
        e = 0

        self.graph = sorted(self.graph, key=lambda item: item[2])
        parent = []
        rank = []

        for node in range(self.V):
            parent.append(node)
            rank.append(0)

        while e < self.V - 1 and i < len(self.graph):
            u, v, w = self.graph[i]
            i += 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            if x != y:
                e += 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)

        total_weight = sum(w for _, _, w in result)
        print(f"Total weight of MST (Kruskal's): {total_weight}")
        return result

def visualize_mst(edges, num_vertices, title):
    """
    Visualize the MST using NetworkX.

    Parameters:
        edges (list): List of edges in the MST.
        num_vertices (int): Number of vertices in the graph.
        title (str): Title for the plot.
    """
    G = nx.Graph()
    
    # Add all vertices
    for i in range(num_vertices):
        G.add_node(i, label=f'Seq {i+1}')
    
    # Add edges
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)
    
    # Create the layout
    pos = nx.spring_layout(G)
    
    # Draw the graph
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=500)
    nx.draw_networkx_edges(G, pos)
    
    # Add labels
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels)
    
    # Add edge labels
    edge_labels = nx.get_edge_attributes(G, 'weight')
    nx.draw_networkx_edge_labels(G, pos, edge_labels)
    
    plt.title(title)
    plt.axis('off')
    plt.show()

def run_and_compare_algorithms(fasta_file):
    """
    Run both Prim's and Kruskal's algorithms and visualize their results.

    Parameters:
        fasta_file (str): Path to the FASTA file.
    """
    # Read sequences
    fasta_data = read_fasta_file(fasta_file)
    num_vertices = len(fasta_data)
    
    print("Sequences from FASTA file:")
    for header, seq in fasta_data:
        print(f"{header}\n{seq}")
    
    # Create edges
    edges = create_distance_edges(fasta_data)
    
    # Run Prim's algorithm
    print("\nRunning Prim's Algorithm:")
    prim_weight, prim_edges = prim_mst(num_vertices, edges)
    print(f"Total weight of MST (Prim's): {prim_weight}")
    print("Edges in Prim's MST:")
    for u, v, w in prim_edges:
        print(f"Sequence {u+1} -- Sequence {v+1} [Weight: {w}]")
    
    # Run Kruskal's algorithm
    print("\nRunning Kruskal's Algorithm:")
    g = Graph(num_vertices)
    for u, v, w in edges:
        g.addEdge(u, v, w)
    kruskal_edges = g.kruskal_mst()
    
    # Visualize both MSTs
    visualize_mst(prim_edges, num_vertices, "Minimum Spanning Tree (Prim's Algorithm)")
    visualize_mst(kruskal_edges, num_vertices, "Minimum Spanning Tree (Kruskal's Algorithm)")

if __name__ == '__main__':
    run_and_compare_algorithms('output.fasta')
