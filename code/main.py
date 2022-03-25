## Libraries 
import random
import time
import re
from Bio import Entrez, SeqIO
from packages import orf_scan_analyze, query_proteins, blast_proteins


if __name__ == "__main__":    
    
    ## Set seed
    random.seed(time.time())

    ## Use Entrez for NCBI querying
    Entrez.email = "neel.mehtani@gmail.com"

    ## Handle to extract sequence data
    handle = Entrez.efetch(db="nucleotide", rettype="gb", id="QGKT01000001")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    ## Extract sequence information from seq object
    seq = record.seq

    ## Parameters 
    min_pro_len = 50 # ORF threshold
    table = 11 # translation reference table
    threshold = 0.3 # BLAST E-value threshold for output logs

    ## Output Directory
    log_dir = "../data"

    ## ORF Scanning, codon translation, and molar mass calculation
    orf_proteins_arr, protein_lens = orf_scan_analyze(seq, log_dir, table=table)

    ## Query Protein Indexing
    q_proteins = query_proteins(orf_proteins_arr, protein_lens)

    ## BLAST Analysis
    blast_proteins(q_proteins, log_dir, e_value_thresh=threshold)