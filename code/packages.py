## Libraries
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import random
import time


## Helper Functions
def orf_scan_analyze(seq, log_dir, min_pro_len=50, table=11):

    ## Lists for ORF protein details
    orf_proteins = []
    protein_wts = []
    protein_lens = []
    
    ## Start scan of both genomic strands
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
    
        ## Scan each possible frame
        for frame in range(3):

            ## Transalte the frame
            for pro in nuc[frame:].translate(table).split("*"):

                if len(pro) >= min_pro_len:

                    ## Run a protein "param" analysis -- returns molecular weight/amino acid freqs etc.
                    analysis = ProteinAnalysis(str(pro))

                    ## Add protein details 
                    orf_proteins.append(pro)
                    protein_wts.append(round(analysis.molecular_weight(), 3))
                    protein_lens.append(len(pro))
                    
                    # Output to screen
                    # print("%s...%s - length %i, strand %i, frame %i" % (pro[:30], pro[-3:], len(pro), strand, frame))

    protein_wts = np.array(protein_wts)
    protein_lens = np.array(protein_lens)
    orf_proteins_arr = np.array(orf_proteins, dtype='object')

    save_file_path = log_dir + "/" + "protein_weights" + '.txt'
    save_file = open(save_file_path, "w")
    save_file.write("Protein Molar Mass Results\n" + "--------------------------------------------" + "\n\n")
    
    for i in range(len(protein_wts)):
        save_file.write("Protein " + str(orf_proteins[i]) + ":  " + str(protein_wts[i]) + "\n\n")

    return orf_proteins_arr, protein_lens

def query_proteins(orf_proteins, protein_lens):
    
    
    query_protein_idxs = random.sample(list(np.where(protein_lens < 1000)[0]), 5)
    query_proteins = orf_proteins[query_protein_idxs]

    return query_proteins

def blast_proteins(query_proteins, log_dir, e_value_thresh=0.04):
    
    #Entrez.email = "neel.mehtani@gmail.com"
    query_ct = 1
    
    for q in query_proteins:
        
        out_file_path = log_dir + "/" + "blast_{}.xml".format(query_ct) 
        out_file = open(out_file_path, "w") 
        
        blast_handle = NCBIWWW.qblast("blastp", "nr", q)
        out_file.write(blast_handle.read())
        out_file.close()
        blast_handle.close()

        result_handle = open(out_file_path)
        blast_records = NCBIXML.read(result_handle)
        
        save_file_path = log_dir + "/" + "BLAST_query_{}_result_log".format(query_ct) + '.txt'
        save_file = open(save_file_path, "w")
        save_file.write("BLAST Query #{} Results -- Sequence Example\n".format(query_ct) + "----------------------------------------" + "\n\n")
        
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < e_value_thresh:

                    save_file.write("sequence: {} \n".format(alignment.title))
                    save_file.write("length: {} \n".format(alignment.length))
                    save_file.write("E Value: {} \n".format(hsp.expect))
                    save_file.write("{} ...  \n".format(hsp.query[0:20]))
                    save_file.write("{} ...  \n".format(hsp.match[0:20]))
                    save_file.write("{} ...  \n".format(hsp.sbjct[0:20]))
                    save_file.write("\n\n")
                else:
                    save_file.write("NULL RESULT -- No hits returned!!")

        result_handle.close()
        query_ct += 1

    return 