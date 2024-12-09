from utils import read_fasta_file, create_distance_edges

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
        """
        Implements Kruskal's algorithm to find the MST.

        Returns:
            list: List of edges in the MST, where each edge is [u, v, weight].
        """
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
    
    # Display the results
    print("\nEdges in the MST:")
    for u, v, weight in kruskal_edges:
        print(f"Sequence {u+1} -- Sequence {v+1} [Weight: {weight}]")
