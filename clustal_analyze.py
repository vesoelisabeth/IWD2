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

tmp_dir = "/home/s2015320/public_html/project2/tmp"
os.makedirs(tmp_dir, exist_ok=True)
os.environ['MPLCONFIGDIR'] = tmp_dir

AA_COLORS = {
    'A': '#FF9999', 'C': '#FFFF99', 'D': '#FF99FF', 'E': '#FF66FF',
    'F': '#99FFFF', 'G': '#CCFF99', 'H': '#9999FF', 'I': '#66FF99',
    'L': '#99FFCC', 'K': '#FFCC99', 'M': '#CC99FF', 'N': '#FFCCFF',
    'P': '#99CCFF', 'Q': '#FF99CC', 'R': '#CCFFFF', 'S': '#FFCCCC',
    'T': '#CCFFCC', 'V': '#99FF99', 'W': '#CCCCFF', 'Y': '#FFFFCC',
    '-': '#FFFFFF'
}

try:
    input_data = json.load(sys.stdin)
except json.JSONDecodeError as e:
    print("Content-type: application/json\n")
    print(json.dumps({"error": "Invalid JSON input: " + str(e)}))
    sys.exit(1)

sequences = input_data.get('sequences', [])

fasta_records = [
    SeqRecord(Seq(seq['sequence']), id=seq['id'], description=seq.get('description', ''))
    for seq in sequences
]

temp_fasta = f"{tmp_dir}/input_" + os.urandom(8).hex() + ".fasta"
temp_aln = f"{tmp_dir}/output_" + os.urandom(8).hex() + ".aln"
output_plot = f"{tmp_dir}/conservation_plot_" + os.urandom(8).hex() + ".png"
msa_plot = f"{tmp_dir}/msa_plot_" + os.urandom(8).hex() + ".png"

with open(temp_fasta, "w") as f:
    SeqIO.write(fasta_records, f, "fasta")

if len(fasta_records) < 2:
    print("Content-type: application/json\n")
    print(json.dumps({
        "result": "Single sequence",
        "conservation_score": 1.0,
        "plot_url": os.path.basename(output_plot),
        "msa_plot_url": os.path.basename(msa_plot),
        "alignment": [{"id": rec.id, "aligned_sequence": str(rec.seq)} for rec in fasta_records]
    }))
    os.remove(temp_fasta)
    sys.exit(0)

clustalo_cmd = [
    "clustalo",
    "-i", temp_fasta,
    "-o", temp_aln,
    "--force",
    "--outfmt", "clustal"
]
result = subprocess.run(clustalo_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print("Content-type: application/json\n")
    print(json.dumps({"error": "Clustal failed: " + result.stderr}))
    os.remove(temp_fasta)
    sys.exit(1)
os.chmod(temp_aln, 0o644)

alignment = AlignIO.read(temp_aln, "clustal")
length = alignment.get_alignment_length()
nums_seqs = len(alignment)

for rec in alignment:
    print(f"DEBUG: Aligned {rec.id}: {str(rec.seq)}", file=sys.stderr)

conservation = []
entropy_scores = []

for i in range(length):
    column = alignment[:, i]
    unique_chars, counts = np.unique(list(column), return_counts=True)
    probs = counts / nums_seqs
    ent = entropy(probs, base=2)

    entropy_scores.append(ent)
    if len(unique_chars) == 1 and '-' not in unique_chars:  # Fixed typo: LENGTH -> len
        conservation.append(1)
    elif '-' in unique_chars:
        conservation.append(0)
    else:
        conservation.append(0.5)

score = sum(conservation) / length

plt.figure(figsize=(10, 4))
sns.lineplot(x=range(length), y=entropy_scores, label="Shannon Entropy")
plt.xlabel("Position")
plt.ylabel("Entropy (bits)")
plt.title("Conservation Profile")
plt.legend()
plt.savefig(output_plot)
os.chmod(output_plot, 0o644)
plt.close()

fig, ax = plt.subplots(figsize=(max(10, length / 10), nums_seqs * 0.5))
ax.set_yticks(range(nums_seqs))
ax.set_yticklabels([rec.id for rec in alignment])
ax.set_xticks(range(0, length, max(1, length // 20)))
ax.set_xlabel("Position")

msa_matrix = np.array([list(rec.seq) for rec in alignment])
color_matrix = np.vectorize(lambda x: AA_COLORS.get(x, '#FFFFFF'))(msa_matrix)

for i in range(nums_seqs):
    for j in range(length):
        ax.fill_between([j, j+1], i, i+1, color=color_matrix[i, j], edgecolor='none')

ax.set_ylim(0, nums_seqs)
ax.set_xlim(0, length)
plt.title("Multiple Sequence Alignment")
plt.savefig(msa_plot)
os.chmod(msa_plot, 0o664)
plt.close()

alignment_data = [{"id": rec.id, "aligned_sequence": str(rec.seq)} for rec in alignment]

os.remove(temp_fasta)
os.remove(temp_aln)

print("Content-type: application/json\n")
print(json.dumps({
    "result": "Alignment completed",
    "conservation_score": round(score, 2),
    "plot_url": os.path.basename(output_plot),
    "msa_plot_url": os.path.basename(msa_plot),
    "alignment": alignment_data
}))
