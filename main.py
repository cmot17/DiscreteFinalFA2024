import networkx as nx
import matplotlib.pyplot as plt
from prims_algorithm import prim_mst
from kruskal_algorithm import Graph
from utils import read_fasta_file, create_distance_edges

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
