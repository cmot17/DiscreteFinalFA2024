import random
import os

def generate_fasta_file(output_filename, num_sequences, max_length):
    """
    Generates a FASTA file with random DNA sequences.

    Parameters:
        output_filename (str): Name of the output .fasta file.
        num_sequences (int): Number of sequences to generate.
        max_length (int): Maximum length of the sequences.
    """
    nucleotides = 'ACGT'

    # Create an empty list to store sequences
    fasta_entries = []

    # Debugging: Start of sequence generation
    print(f"Generating {num_sequences} sequences, each with a max length of {max_length}")

    for i in range(1, num_sequences + 1):
        # Random length for each sequence
        seq_length = random.randint(1, max_length)
        sequence = ''.join(random.choices(nucleotides, k=seq_length))
        header = f">sequence_{i}"

        # Debugging: Print the header and sequence length
        print(f"Sequence {i}: {header}, Length: {seq_length}")

        fasta_entries.append((header, sequence))

    # Check if the file already exists and create a new file if it doesn't exist
    file_path = os.path.join(os.getcwd(), output_filename)

    # Debugging: Start writing to file
    print(f"Writing sequences to {file_path}")

    try:
        # Open the output file in write mode to create the file (not overwrite if exists)
        with open(file_path, mode='w', encoding='utf-8') as fasta_file:
            for header, sequence in fasta_entries:
                fasta_file.write(f"{header}\n")  # Write the header (sequence identifier)
                # Write the sequence in chunks of 80 characters (FASTA format standard)
                for i in range(0, len(sequence), 80):
                    fasta_file.write(f"{sequence[i:i+80]}\n")
        # Debugging: Successful file writing
        print(f"FASTA file '{file_path}' with {num_sequences} sequences generated successfully!")

    except Exception as e:
        print(f"Error writing to file: {e}")


# Constants (renamed for UPPER_CASE convention)
OUTPUT_FILE = "output.fasta"  # File will be saved in the same folder as the script
NUM_SEQUENCES = 10
MAX_SEQ_LENGTH = 100

# Generate the FASTA file
generate_fasta_file(OUTPUT_FILE, NUM_SEQUENCES, MAX_SEQ_LENGTH)
