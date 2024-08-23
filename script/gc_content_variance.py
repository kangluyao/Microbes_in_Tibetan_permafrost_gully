from Bio import SeqIO
from Bio.SeqUtils import GC
import os
import numpy as np

def calculate_gc(fasta_file):
    """Calculate GC content for each sequence in a FASTA file."""
    gc_contents = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gc_contents.append(GC(record.seq))
    return gc_contents

def process_multiple_files(fasta_files):
    """Process multiple FASTA files to calculate GC content and variance."""
    results = {}
    for file in fasta_files:
        gc_contents = calculate_gc(file)
        gc_variance = np.var(gc_contents)
        results[file] = {
            'GC_Content': gc_contents,
            'GC_Variance': gc_variance
        }
    return results

def main(directory):
    """Main function to process all FASTA files in a directory."""
    fasta_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fasta")]
    results = process_multiple_files(fasta_files)
    
    # Output the results
    for file, data in results.items():
        print(f"File: {file}")
        print(f"GC Content: {data['GC_Content']}")
        print(f"GC Variance: {data['GC_Variance']}")
        print("-" * 40)

if __name__ == "__main__":
    # Specify the directory containing your metagenomic FASTA files
    directory = "/data01/kangluyao/thermokarst_gully/temp/qc/mergerd"
    main(directory)
