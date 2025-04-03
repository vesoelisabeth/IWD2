#!/usr/bin/python3
import os
import subprocess
from Bio import Entrez, SeqIO

# Configuration
Entrez.email = "vesoelisabeh@gmail.com"
OUTPUT_DIR = "/home/s2015320/public_html/project2/output/"
JOB_ID = "example_g6pase_aves"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Fetch G6Pase sequences from Aves (TaxID: 8782)
print("Fetching glucose-6-phosphatase sequences from Aves...")
query = "glucose-6-phosphatase[Protein Name] txid8782[Organism]"
handle = Entrez.esearch(db="protein", term=query, retmax=10)
record = Entrez.read(handle)
ids = record["IdList"]
if not ids:
    print("No sequences found for glucose-6-phosphatase in Aves.")
    exit(1)
handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
fasta_file = f"{OUTPUT_DIR}input_{JOB_ID}.fasta"
with open(fasta_file, "w") as f:
    f.write(handle.read())
handle.close()
print(f"FASTA file saved: {fasta_file}")

# MSA with Clustal Omega
print("Running Clustal Omega for MSA...")
aln_file = f"{OUTPUT_DIR}input_{JOB_ID}.aln"
clustalo_cmd = ["clustalo", "-i", fasta_file, "-o", aln_file, "--force", "--outfmt=clustal"]
result = subprocess.run(clustalo_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"Clustal Omega failed: {result.stderr}")
    exit(1)
print(f"Alignment file saved: {aln_file}")

# PatMatMotif scan (assuming local install)
print("Running PatMatMotif for motif scanning...")
motif_file = f"{OUTPUT_DIR}motifs_{JOB_ID}.txt"
patmat_cmd = ["patmatmotifs", "-sequence", fasta_file, "-outfile", motif_file, "-auto"]
result = subprocess.run(patmat_cmd, capture_output=True, text=True)
if result.returncode != 0:
    print(f"PatMatMotif failed: {result.stderr}")
    exit(1)
print(f"Motif file saved: {motif_file}")

# Advanced Analysis (reusing your script)
print("Running advanced analysis...")
advanced_cmd = [
    "/home/s2015320/public_html/project2/myenv/bin/python3",
    "/home/s2015320/public_html/project2/advanced_analysis.py",
    aln_file  # Pass alignment file explicitly
]
env = os.environ.copy()
env["MPLCONFIGDIR"] = "/tmp/matplotlib_cache"
env["SESSION_ID"] = JOB_ID
result = subprocess.run(advanced_cmd, env=env, capture_output=True, text=True)
if result.returncode != 0:
    print(f"Advanced analysis failed: {result.stderr}")
    exit(1)
print(f"Advanced analysis output: {result.stdout}")

# Verify outputs
expected_files = [
    f"input_{JOB_ID}.fasta",
    f"input_{JOB_ID}.aln",
    f"motifs_{JOB_ID}.txt",
    f"report_{JOB_ID}.html",
    f"phylogeny_{JOB_ID}.png",
    f"ss_{JOB_ID}.garnier"
]
for file in expected_files:
    full_path = f"{OUTPUT_DIR}{file}"
    if os.path.exists(full_path):
        print(f"Generated: {full_path}")
    else:
        print(f"Missing: {full_path}")

print(f"Example dataset generated for job_id: {JOB_ID}")
