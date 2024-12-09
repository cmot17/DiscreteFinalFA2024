import random
import os


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
    # Read FASTA sequences from file
    fasta_data = read_fasta_file('output.fasta')

    # Print the sequences from file
    print("Sequences from FASTA file:")
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
