#!/usr/bin/env python3
# biopython_connect.py

from Bio import Entrez, SeqIO
import sys
import json

Entrez.email = "vesoelisabeth@gmail.com"  

def search_protein(protein_name, taxonomy):
    """Search NCBI for protein sequences and return GenBank data."""
    try:
        # Search NCBI protein database with the Entrez command
        handle = Entrez.esearch(db="protein", term=f"{protein_name}[Protein Name] {taxonomy}[Organism]", retmax=10)
        record = Entrez.read(handle)
        handle.close()
        
        protein_ids = record['IdList']
        if not protein_ids:
            return {"result": "No protein matches found.", "sequences": []}
        
        # Fetch GenBank data
        handle = Entrez.efetch(db="protein", id=protein_ids, rettype="gb", retmode="text")
        records = SeqIO.parse(handle, "genbank") #Processes GenBank data for the next steps
        
        # Extract relevant data
        sequences = []
        for rec in records:
            seq_data = {
                "id": rec.id,                  # Accession number
                "description": rec.description, # Protein name and details
                "sequence": str(rec.seq),      # Amino acid sequence
                "organism": rec.annotations.get("organism", "Unknown")  # Taxonomy info
            }
            sequences.append(seq_data)
        handle.close()
        
        return {"result": f"Found {len(sequences)} sequences.", "sequences": sequences}
    
    except Exception as e:
        return {"error": f"NCBI query error: {e}"} #Error trap

#Outputs a JSON object 
if __name__ == "__main__":
    try:
        input_data = sys.stdin.read()
        if not input_data:
            print("Content-type: application/json\n")
            print(json.dumps({"error": "No input data provided"}))
            sys.exit(0)
        
        data = json.loads(input_data)
        protein_name = data.get('protein_name')
        taxonomy = data.get('taxonomy')
        
        if not protein_name or not taxonomy:
            result = {"error": "Missing protein_name or taxonomy"}
        else:
            result = search_protein(protein_name, taxonomy)
        
        print("Content-type: application/json\n")
        print(json.dumps(result))
    
    except json.JSONDecodeError:
        print("Content-type: application/json\n")
        print(json.dumps({"error": "Invalid JSON input"}))
    except Exception as e:
        print("Content-type: application/json\n")
        print(json.dumps({"error": f"Unexpected error: {e}"}))


#Bipython lecture BPSM
