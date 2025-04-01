#!/usr/bin/python3
import os
import subprocess
import sys
import time
from Bio import Entrez, SeqIO, AlignIO
from Bio.PDB import PDBList
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# --- Configuration ---
Entrez.email = "vesoelisabeh@gmail.com"
OUTPUT_DIR = "/tmp/advanced_analysis_" + os.urandom(8).hex() + "/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Secondary Structure Prediction ---
def analyze_secondary_structure(sequence_file):
    ss_file = os.path.join(OUTPUT_DIR, "ss.garnier")
    try:
        result = subprocess.run(["garnier", "-sequence", sequence_file, "-outfile", ss_file], capture_output=True, text=True, check=True)
        print(f"Garnier output: {result.stdout}")
        return parse_garnier(ss_file)
    except subprocess.CalledProcessError as e:
        return [f"Error in garnier: {e.stderr}"]

def parse_garnier(garnier_file):
    if not os.path.exists(garnier_file):
        return ["No output from garnier"]
    with open(garnier_file) as f:
        return [line.strip() for line in f]

# --- 3D Structure Homology ---
def fetch_structures(sequence):
    pdbl = PDBList()
    try:
        time.sleep(1)  # Rate limit delay
        # Search PDB with sequence ID or protein name, link to structure
        handle = Entrez.esearch(db="structure", term="insulin", retmax=3)
        pdb_ids = Entrez.read(handle)["IdList"]
        structures = []
        for pdb_id in pdb_ids:
            time.sleep(1)  # Delay between downloads
            try:
                pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=OUTPUT_DIR, file_format="pdb")
                resolution = get_resolution(pdb_file)
                structures.append({"id": pdb_id, "path": pdb_file, "resolution": resolution})
            except Exception as e:
                structures.append({"id": pdb_id, "path": "", "resolution": f"Download error: {e}"})
        return structures if structures else [{"id": "N/A", "path": "", "resolution": "No structures found"}]
    except Exception as e:
        return [{"id": "N/A", "path": "", "resolution": f"Error: {e}"}]

def get_resolution(pdb_file):
    if not os.path.exists(pdb_file):
        return "N/A"
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("REMARK   2 RESOLUTION."):
                return line.split()[-2] if line.split()[-2] != "ANGSTROMS" else "N/A"
    return "N/A"

# --- Phylogenetic Tree Generation ---
def generate_phylogeny(alignment_file):
    try:
        alignment = AlignIO.read(alignment_file, "clustal")
        if len(alignment) < 2:
            return "Error: Need at least 2 sequences for phylogeny"
        seq_set = {str(record.seq) for record in alignment}
        if len(seq_set) < 2:
            return "Error: All sequences are identical"
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        # Build tree
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree_file = os.path.join(OUTPUT_DIR, "phylogeny.png")
        plt.figure(figsize=(10, 8))
        Phylo.draw(tree, do_show=False)
        plt.savefig(tree_file)
        plt.close()
        return tree_file
    except Exception as e:
        return f"Error: {e}"

# --- Report Generation ---
def generate_report(analysis_results):
    html_lines = [
        "<html><head><title>Protein Insights</title>",
        "<style>body { font-family: Arial; } pre { white-space: pre-wrap; }</style>",
        "</head><body>",
        "<h1>Advanced Protein Insights</h1>",
        "<h2>Secondary Structure</h2><pre>",
        "\n".join(analysis_results['secondary_structure']),
        "</pre>",
        "<h2>3D Structures</h2><ul>",
        "".join(f'<li>{s["id"]} (Resolution: {s["resolution"]})</li>' for s in analysis_results['structures']),
        "</ul>",
        "<h2>Phylogeny</h2>",
        f'<img src="{analysis_results["phylogeny"]}" width="50%">' if os.path.exists(analysis_results["phylogeny"]) else f"No phylogeny: {analysis_results['phylogeny']}",
        "</body></html>"
    ]
    html = "\n".join(html_lines)

    report_file = os.path.join(OUTPUT_DIR, "report.html")
    with open(report_file, "w") as f:
        f.write(html)
    return report_file

# --- Main Execution ---
def main(alignment_file):
    fasta_file = os.path.join(OUTPUT_DIR, "input.fasta")
    AlignIO.convert(alignment_file, "clustal", fasta_file, "fasta")
    sequence = list(SeqIO.parse(fasta_file, "fasta"))[0]
    
    results = {
        "secondary_structure": analyze_secondary_structure(fasta_file),
        "structures": fetch_structures(sequence),
        "phylogeny": generate_phylogeny(alignment_file)
    }
    return generate_report(results)

# --- Command-Line Entry Point ---
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python advanced_analysis.py <alignment.clustal>")
        sys.exit(1)
    report = main(sys.argv[1])
    print(f"Report generated: {report}")
