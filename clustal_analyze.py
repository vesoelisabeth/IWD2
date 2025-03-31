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

# Define amino acid color scheme (simplified, inspired by DECIPHER/Taylor scheme)
AA_COLORS = {
    'A': '#FF9999', 'C': '#FFFF99', 'D': '#FF99FF', 'E': '#FF66FF',
    'F': '#99FFFF', 'G': '#CCFF99', 'H': '#9999FF', 'I': '#66FF99',
    'L': '#99FFCC', 'K': '#FFCC99', 'M': '#CC99FF', 'N': '#FFCCFF',
    'P': '#99CCFF', 'Q': '#FF99CC', 'R': '#CCFFFF', 'S': '#FFCCCC',
    'T': '#CCFFCC', 'V': '#99FF99', 'W': '#CCCCFF', 'Y': '#FFFFCC',
    '-': '#FFFFFF'  # Gap color
}

# Read JSON input from stdin
input_data = json.load(sys.stdin)
sequences = input_data['sequences']

# Convert sequences to FASTA format in memory
fasta_records = [
    SeqRecord(Seq(seq['sequence']), id=seq['id'], description=seq.get('description', ''))
    for seq in sequences
]

# Write to temporary FASTA file
temp_fasta = "/home/s2015320/public_html/project2/clustal_input_" + os.urandom(8).hex() + ".fasta"
temp_aln = "/home/s2015320/public_html/project2/clustal_output_" + os.urandom(8).hex() + ".aln"
output_plot = os.getenv('PLOT_PATH', '/tmp/conservation_plot.png')
msa_plot = output_plot.replace('conservation_plot', 'msa_plot')

with open(temp_fasta, "w") as f:
    SeqIO.write(fasta_records, f, "fasta")

# Check number of sequences
if len(fasta_records) < 2:
    print("Content-type: application/json\n")
    print(json.dumps({
        "result": "Single sequence provided",
        "conservation_score": 1.0,
        "plot_url": None,
	"msa_plot_url": None,
	"alignment": [{"id": rec.id, "aligned_sequence": str(rec.seq)} for rec in fasta_records]
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
os.chmod(temp_aln, 0o644)

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
plt.savefig(output_plot, 0o644)
plt.close()

# Generate MSA plot
fig, ax = plt.subplots(figsize=(max(10, length / 10), nums_seqs * 0.5))
ax.set_yticks(range(nums_seqs))
ax.set_yticklabels([rec.id for rec in alignment])
ax.set_xticks(range(0, length, max(1, length // 20)))
ax.set_xlabel("Position")

# Create a matrix of colors
msa_matrix = np.array([list(rec.seq) for rec in alignment])
color_matrix = np.vectorize(lambda x: AA_COLORS.get(x, '#FFFFFF'))(msa_matrix)

# Plot each residue as a colored square
for i in range(nums_seqs):
    for j in range(length):
        ax.fill_between([j, j+1], i, i+1, color=color_matrix[i, j], edgecolor='none')

ax.set_ylim(0, nums_seqs)
ax.set_xlim(0, length)
plt.title("Multiple Sequence Alignment")
plt.savefig(msa_plot)
os.chmod(msa_plot, 0o664)
plt.close()

# Prepare alignment data
alignment_data = [{"id": rec.id, "aligned_sequence": str(rec.seq)} for rec in alignment]

# Clean up temporary files
os.remove(temp_fasta)
os.remove(temp_aln)


# Output JSON response
plot_url = f"https://bioinfmsc8.bio.ed.ac.uk/~s2015320/project2/{os.path.basename(output_plot)}"
print("Content-type: application/json\n")
print(json.dumps({
    "result": "Alignment completed",
    "conservation_score": round(score, 2),
    "plot_url": plot_url
    "msa_plot_url": msa_plot_url,
    "alignment": alignment_data
}))
	
