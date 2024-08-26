from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import numpy as np

def calculate_gc(fasta_file):
    """Calculate GC content for each sequence in a FASTA file."""
    gc_contents = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        gc_contents.append(gc_fraction(record.seq) * 100)  # Convert fraction to percentage
    return gc_contents

def process_multiple_files(fasta_files):
    """Process multiple FASTA files to calculate GC content and variance."""
    results = []
    for file in fasta_files:
        gc_contents = calculate_gc(file)
        mean_gc = np.mean(gc_contents)
        gc_variance = np.var(gc_contents)
        results.append({
            'File': os.path.basename(file),
            'Mean_GC_Content': mean_gc,
            'GC_Variance': gc_variance
        })
    return results

def write_results_to_file(results, output_file):
    """Write the GC content results to a text file."""
    with open(output_file, 'w') as f:
        f.write("Metagenomics Name\tMean GC Content (%)\tGC Variance\n")
        for result in results:
            f.write(f"{result['File']}\t{result['Mean_GC_Content']:.2f}\t{result['GC_Variance']:.2f}\n")

def main(directory, output_file):
    """Main function to process all FASTA files in a directory and write results to a file."""
    fasta_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fasta")]
    results = process_multiple_files(fasta_files)
    write_results_to_file(results, output_file)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    # Specify the directory containing your metagenomic FASTA files
    directory = "/data01/kangluyao/thermokarst_gully/temp/qc/merged"
    # Specify the output file path
    output_file = "gc_content_results.txt"
    main(directory, output_file)
