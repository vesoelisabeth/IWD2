import sys
import subprocess
from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy

# Get input and output file paths from command line
input_fasta = sys.argv[1]
output_aln = sys.argv[2]
output_plot = sys.argv[3]

# Run Clustal Omega
clustalo_cmd = [
    "clustalo",
    "-i", input_fasta,
    "-o", output_aln,
    "--force",
    "--outfmt=clustal"
]
subprocess.run(clustalo_cmd, check=True)

# Read alignment and compute conservation (simple example)
alignment = AlignIO.read(output_aln, "clustal")
length = alignment.get_alignment_length()
nums_seqs = len(alignment)

# Calculate the Shannon entropy and conservation per column 
conservation = []
entropy_scores = []

for i in range(length):
    column = alignment[:, i]
    unique_chars, counts = np.unique(list(column), return_counts=True)
    if len(unique_chars) == 1 and '-' not in unique_chars:
        conservation.append(1)  # Fully conserved
    elif '-' in unique_chars:
        conservation.append(0)  # Gap
    else:
        conservation.append(0.5)  # Partially conserved

# Output conservation score (could be enhanced)
score = sum(conservation) / length

# Generate conservation plot
plt.figure(figsize=(10, 4))
plt.plot(range(length), entropy_scores, label="Shannon Entropy")
plt.xlabel("Position")
plt.ylabel("Entropy (bits)")
plt.title("Conservation Profile")
plt.legend()
plt.savefig(output_plot)
plt.close()

with open(output_aln, 'a') as f:
    f.write(f"\nConservation Score: {score:.2f}\n")

#BPSM Clustalo information from the Bioinformatic web lecture
#https://www.geeksforgeeks.org/python-numpy-np-unique-method/
#https://stackoverflow.com/questions/63576206/get-unique-values-in-a-list-of-numpy-arrays
#https://biopython.org/docs/dev/Tutorial/chapter_msa.html 
