import random
import os


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


# Kruskal's Algorithm Implementation
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

    def KruskalMST(self):
        result = []
        i = 0
        e = 0

        self.graph = sorted(self.graph, key=lambda item: item[2])
        parent = []
        rank = []

        for node in range(self.V):
            parent.append(node)
            rank.append(0)

        while e < self.V - 1:
            u, v, w = self.graph[i]
            i += 1
            x = self.find(parent, u)
            y = self.find(parent, v)

            if x != y:
                e += 1
                result.append([u, v, w])
                self.union(parent, rank, x, y)

        print("Edges in the constructed MST:")
        for u, v, weight in result:
            print(f"{u} -- {v} == {weight}")


# Driver Code
if __name__ == '__main__':
    # Generate FASTA sequences
    num_sequences = 25
    max_length = 50
    fasta_data = generate_fasta_sequences(num_sequences, max_length)

    # Print the generated sequences
    print("Generated FASTA Sequences:")
    for header, seq in fasta_data:
        print(f"{header}\n{seq}")

    # Create a graph from the sequences
    g = Graph(len(fasta_data))
    for i in range(len(fasta_data)):
        for j in range(i + 1, len(fasta_data)):
            seq1 = fasta_data[i][1]
            seq2 = fasta_data[j][1]
            distance = hamming_distance(seq1, seq2)
            g.addEdge(i, j, distance)

    # Apply Kruskal's algorithm to find MST
    g.KruskalMST()
