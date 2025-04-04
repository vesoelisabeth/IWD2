#!/home/s2015320/public_html/project2/myenv/bin/python3

import sys
import os
import time
import urllib.request
from Bio import AlignIO, SeqIO, Phylo
from Bio.PDB import PDBList
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Align.Applications import ClustalOmegaCommandline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SESSION_ID = os.environ.get("SESSION_ID", "example_g6pase_aves")
TEMP_DIR = os.environ.get("TEMP_DIR", "/tmp")

# --- 3D Structures Section ---
def fetch_structures(sequence):
    pdbl = PDBList()
    try:
        seq_str = str(sequence.seq)
        seq_id = sequence.id
        print(f"Searching structures for sequence ID: {seq_id}")
        blast_handle = NCBIWWW.qblast("blastp", "pdb", seq_str, hitlist_size=3)
        blast_record = NCBIXML.read(blast_handle)
        print(f"BLAST alignments found: {len(blast_record.alignments)}")
        structures = []
        for alignment in blast_record.alignments[:3]:
            pdb_id_full = alignment.accession.split("|")[1] if "|" in alignment.accession else alignment.accession
            pdb_id = pdb_id_full.split("_")[0].upper()
            print(f"Found PDB hit: {pdb_id}")
            time.sleep(1)
            temp_file = os.path.join(TEMP_DIR, f"{pdb_id}_{SESSION_ID}.pdb")
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            try:
                urllib.request.urlretrieve(url, temp_file)
                os.chmod(temp_file, 0o644)
                resolution = get_resolution(temp_file)
                print(f"Downloaded PDB {pdb_id} to {temp_file}")
                structures.append({"id": pdb_id, "path": temp_file, "resolution": resolution})
            except Exception as e:
                print(f"Failed to download {pdb_id}: {e}")
                structures.append({"id": pdb_id, "path": "", "resolution": f"Download error: {e}"})
        return structures if structures else [{"id": "N/A", "path": "", "resolution": "No structures found"}]
    except Exception as e:
        return [{"id": "N/A", "path": "", "resolution": f"Error: {e}"}]

def get_resolution(pdb_file):
    try:
        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith("REMARK   2 RESOLUTION"):
                    return line.split()[-2]
        return "Unknown"
    except Exception as e:
        return f"Resolution error: {e}"

# --- Garnier Section ---
def run_garnier(alignment_file):
    garnier_file = os.path.join(TEMP_DIR, f"ss_{SESSION_ID}.garnier")
    try:
        os.system(f"garnier -sequence {alignment_file} -outfile {garnier_file}")
        print("Garnier output generated")
        return garnier_file
    except Exception as e:
        print(f"Garnier failed: {e}")
        return None

# --- Phylogeny Section ---
def generate_phylogeny(alignment_file):
    phylo_file = os.path.join(TEMP_DIR, f"phylogeny_{SESSION_ID}.png")
    try:
        alignment = AlignIO.read(alignment_file, "clustal")
        tree = Phylo.read(os.path.join(TEMP_DIR, "temp.dnd"), "newick")
        plt.figure(figsize=(10, 5))
        Phylo.draw(tree, do_show=False)
        plt.savefig(phylo_file)
        plt.close()
        return phylo_file
    except Exception as e:
        print(f"Phylogeny generation failed: {e}")
        return None

# --- Report Section ---
def generate_report(alignment_file, garnier_file, structures, phylo_file):
    report_file = os.path.join(TEMP_DIR, "report.html")
    try:
        with open(report_file, "w") as f:
            f.write("<html><head><title>Advanced Analysis</title></head><body>")
            f.write("<h1>Advanced Analysis Report</h1>")
            f.write("<h2>Secondary Structure (Garnier)</h2>")
            if garnier_file and os.path.exists(garnier_file):
                with open(garnier_file, "r") as gf:
                    f.write("<pre>" + gf.read() + "</pre>")
            else:
                f.write("<p>No secondary structure data available.</p>")
            f.write("<h2>3D Structures</h2>")
            for struct in structures:
                f.write(f"<p>{struct['id']} (Resolution: {struct['resolution']})</p>")
            f.write("<h2>Phylogeny</h2>")
            if phylo_file and os.path.exists(phylo_file):
                f.write(f"<img src='{os.path.basename(phylo_file)}' alt='Phylogeny Tree'>")
            else:
                f.write("<p>No phylogeny image available.</p>")
            f.write("</body></html>")
        print(f"Report generated: {report_file}")
        return report_file
    except Exception as e:
        print(f"Report generation failed: {e}")
        return None

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 advanced_analysis.py <alignment_file>")
        sys.exit(1)
    
    alignment_file = sys.argv[1]
    if not os.path.exists(alignment_file):
        print(f"Alignment file {alignment_file} not found")
        sys.exit(1)
    
    try:
        alignment = AlignIO.read(alignment_file, "clustal")
    except Exception as e:
        print(f"Failed to read alignment: {e}")
        sys.exit(1)
    
    garnier_file = run_garnier(alignment_file)
    structures = fetch_structures(alignment[0])
    phylo_file = generate_phylogeny(alignment_file)
    report_file = generate_report(alignment_file, garnier_file, structures, phylo_file)
    
    if not report_file:
        print("Advanced analysis failed")
        sys.exit(1)

if __name__ == "__main__":
    main()
#https://biopython.org/docs/latest/Tutorial/chapter_msa.html

