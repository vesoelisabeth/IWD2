#!/usr/bin/python3
import os
import subprocess
import sys
import time
import urllib.request
from Bio import Entrez, SeqIO, AlignIO
from Bio.PDB import PDBList
from Bio.Blast import NCBIWWW, NCBIXML
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# --- Set MPLCONFIGDIR for Matplotlib ---
if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = "/tmp/matplotlib_cache"
if not os.path.exists(os.environ["MPLCONFIGDIR"]):
    try:
        os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
    except FileExistsError:
        pass

# --- Configuration ---
Entrez.email = "vesoelisabeh@gmail.com"
SESSION_ID = os.environ.get("SESSION_ID", os.urandom(8).hex())
OUTPUT_DIR = "/home/s2015320/public_html/project2/output/"
os.makedirs(OUTPUT_DIR, exist_ok=True)  # Should already exist

# --- Secondary Structure Prediction ---
def analyze_secondary_structure(sequence_file):
    ss_file = os.path.join(OUTPUT_DIR, f"ss_{SESSION_ID}.garnier")
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
        seq_str = str(sequence.seq)
        seq_id = sequence.id
        print(f"Searching structures for sequence ID: {seq_id}")
        blast_handle = NCBIWWW.qblast("blastp", "pdb", seq_str, hitlist_size=3)
	print(f"BLAST alignments found: {len(blast_record.alignments)}")
        blast_record = NCBIXML.read(blast_handle)
        structures = []
        for alignment in blast_record.alignments[:3]:
            pdb_id_full = alignment.accession.split("|")[1] if "|" in alignment.accession else alignment.accession
            pdb_id = pdb_id_full.split("_")[0].upper()
            print(f"Found PDB hit: {pdb_id}")
            time.sleep(1)
            try:
                pdb_file = pdbl.retrieve_pdb_file(pdb_id, pdir=OUTPUT_DIR, file_format="pdb")
                if not os.path.exists(pdb_file):
                    raise Exception("Bio.PDB failed to download")
                resolution = get_resolution(pdb_file)
                print(f"Downloaded PDB {pdb_id} via Bio.PDB")
            except Exception as e:
                print(f"Bio.PDB failed: {e}, trying direct download")
                pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}.pdb")
                url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
                try:
                    urllib.request.urlretrieve(url, pdb_file)
                    resolution = get_resolution(pdb_file)
                    print(f"Downloaded PDB {pdb_id} directly from RCSB")
                except Exception as e2:
                    structures.append({"id": pdb_id, "path": "", "resolution": f"Download error: {e2}"})
                    continue
            structures.append({"id": pdb_id, "path": pdb_file, "resolution": resolution})
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
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        tree_file = os.path.join(OUTPUT_DIR, f"phylogeny_{SESSION_ID}.png")
        plt.figure(figsize=(10, 8))
        Phylo.draw(tree, do_show=False)
        plt.savefig(tree_file)
        plt.close()
        return f"phylogeny_{SESSION_ID}.png"
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
        f'<img src="/~s2015320/project2/output/{analysis_results["phylogeny"]}" width="50%">' if os.path.exists(os.path.join(OUTPUT_DIR, analysis_results["phylogeny"])) else f"No phylogeny: {analysis_results['phylogeny']}",
        "</body></html>"
    ]
    html = "\n".join(html_lines)

    report_file = os.path.join(OUTPUT_DIR, f"report_{SESSION_ID}.html")
    with open(report_file, "w") as f:
        f.write(html)
    return report_file

# --- Main Execution ---
def main(alignment_file):
    fasta_file = os.path.join(OUTPUT_DIR, f"input_{SESSION_ID}.fasta")
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
