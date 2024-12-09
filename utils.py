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