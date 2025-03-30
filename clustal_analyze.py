#!/usr/bin/python3
import os
import sys
import json
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import entropy

# Set MPLCONFIGDIR to avoid Matplotlib warning
os.environ['MPLCONFIGDIR'] = '/localdisk/home/s2015320/public_html/project2/matplotlib_cache'

# Read JSON input from stdin
input_data = json.load(sys.stdin)
sequences = input_data['sequences']

# Convert sequences to FASTA format in memory
fasta_records = [
    SeqRecord(Seq(seq['sequence']), id=seq['id'], description=seq.get('description', ''))
    for seq in sequences
]

# Write to temporary FASTA file
temp_fasta = "/tmp/clustal_input.fasta"
temp_aln = "/tmp/clustal_output.aln"
output_plot = "/localdisk/home/s2015320/public_html/project2/conservation_plot.png"

with open(temp_fasta, "w") as f:
    SeqIO.write(fasta_records, f, "fasta")

# Check number of sequences
if len(fasta_records) < 2:
    print("Content-type: application/json\n")
    print(json.dumps({
        "result": "Single sequence provided",
        "conservation_score": 1.0,
        "plot_url": None
    }))
    os.remove(temp_fasta)
    sys.exit(0)

# Run Clustal Omega
clustalo_cmd = [
    "clustalo",
    "-i", temp_fasta,
    "-o", temp_aln,
    "--force",
    "--outfmt=clustal"
]
subprocess.run(clustalo_cmd, check=True)

# Read alignment and compute conservation
alignment = AlignIO.read(temp_aln, "clustal")
length = alignment.get_alignment_length()
nums_seqs = len(alignment)

# Calculate Shannon entropy and conservation per column
conservation = []
entropy_scores = []

for i in range(length):
    column = alignment[:, i]
    unique_chars, counts = np.unique(list(column), return_counts=True)
    probs = counts / nums_seqs
    ent = entropy(probs, base=2)
    entropy_scores.append(ent)
    if len(unique_chars) == 1 and '-' not in unique_chars:
        conservation.append(1)  # Fully conserved
    elif '-' in unique_chars:
        conservation.append(0)  # Gap
    else:
        conservation.append(0.5)  # Partially conserved

# Conservation score
score = sum(conservation) / length

# Generate conservation plot
plt.figure(figsize=(10, 4))
sns.lineplot(x=range(length), y=entropy_scores, label="Shannon Entropy")
plt.xlabel("Position")
plt.ylabel("Entropy (bits)")
plt.title("Conservation Profile")
plt.legend()
plt.savefig(output_plot)
plt.close()

# Clean up temporary files
os.remove(temp_fasta)
os.remove(temp_aln)

# Output JSON response
print("Content-type: application/json\n")
print(json.dumps({
    "result": "Alignment completed",
    "conservation_score": round(score, 2),
    "plot_url": "https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/conservation_plot.png"
}))
