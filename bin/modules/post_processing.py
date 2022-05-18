import os
from re import T
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
import pandas as pd
import numpy as np
from io import StringIO

pd.options.mode.chained_assignment = None

def get_total_read_count(out_dir,prefix):

    flagstat_file =  os.path.join(out_dir,prefix+'_flagstat.txt')

    # read first line 
    with open(flagstat_file) as f:
        first_line = f.readline()
    # split and strip
    spl_list = first_line.split("+")
    total_reads = spl_list[0].strip()

    return total_reads


def parse_bam(out_dir, prefix):
    sam_file = (os.path.join(out_dir,prefix+'_mapped.sam'))
    with open(sam_file, 'r') as temp_f:
        # get No of columns in each line
        col_count = [ len(l.split(",")) for l in temp_f.readlines() ]
    
    # ### Generate column names  (names will be 0, 1, 2, ..., maximum columns - 1)
    column_names = [i for i in range(0, max(col_count))]

    df = pd.read_csv(sam_file, header=None, delimiter="\t", names=column_names)

    print(df)



  
def create_gff(phanotate_mmseqs_df, length_df, fasta_input, out_dir):
    # write the headers of the gff file
    with open(os.path.join(out_dir, "phrokka.gff"), 'w') as f:
        f.write('##gff-version 3\n')
        for index, row in length_df.iterrows():
            f.write('##sequence-region ' + row['contig'] + ' 1 ' + str(row['length']) +'\n')
  
    # rearrange start and stop so that start is always less than stop for gff
    cols = ["start","stop"]
    #indices where start is greater than stop
    ixs = phanotate_mmseqs_df['frame'] == '-'
    # Where ixs is True, values are swapped
    phanotate_mmseqs_df.loc[ixs,cols] = phanotate_mmseqs_df.loc[ixs, cols].reindex(columns=cols[::-1]).values
    
    phanotate_mmseqs_df['phase'] = 0
    phanotate_mmseqs_df['attributes'] = "phrog=" + phanotate_mmseqs_df["phrog"] + ";" + "top_hit=" + phanotate_mmseqs_df["top_hit"] + ";" + "annotation=" + phanotate_mmseqs_df["annot"] + ";" + "function=" + phanotate_mmseqs_df["category"]

    # get gff dataframe in correct order 
    gff_df = phanotate_mmseqs_df[["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]]

    with open(os.path.join(out_dir, "phrokka.gff"), 'a') as f:
        gff_df.to_csv(f, sep="\t", index=False, header=False)

      
    ### trnas

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
    trna_df = pd.read_csv(os.path.join(out_dir,"trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
    # keep only trnas
    trna_df = trna_df[trna_df['Region'] == 'tRNA']
    trna_df.start = trna_df.start.astype(int)
    trna_df.stop = trna_df.stop.astype(int)
    with open(os.path.join(out_dir, "phrokka.gff"), 'a') as f:
        trna_df.to_csv(f, sep="\t", index=False, header=False)


    # write fasta on the end 

    ##FASTA
    with open(os.path.join(out_dir, "phrokka.gff"), 'a') as f:
        f.write('##FASTA\n')
        fasta_sequences = SeqIO.parse(open(fasta_input),'fasta')
        SeqIO.write(fasta_sequences, f, "fasta")

def create_tbl(phanotate_mmseqs_df, length_df, out_dir):

    ### readtrnas

    col_list = ["contig", "Method", "Region", "start", "stop", "score", "frame", "phase", "attributes"]
     # check if no trnas
    empty = False
    if os.stat(os.path.join(out_dir, "trnascan_out.gff")).st_size == 0:
        empty = True
    if empty == False:    
        trna_df = pd.read_csv(os.path.join(out_dir, "trnascan_out.gff"), delimiter= '\t', index_col=False, names=col_list ) 
        # keep only trnas
        trna_df = trna_df[trna_df['Region'] == 'tRNA']
        trna_df.start = trna_df.start.astype(int)
        trna_df.stop = trna_df.stop.astype(int)

        trna_df[['attributes','isotypes']] = trna_df['attributes'].str.split(';isotype=',expand=True)
        trna_df[['isotypes','anticodon']] = trna_df['isotypes'].str.split(';anticodon=',expand=True)
        trna_df[['anticodon','rest']] = trna_df['anticodon'].str.split(';gene_biotype',expand=True)
        trna_df['trna_product']='tRNA-'+trna_df['isotypes']+"("+trna_df['anticodon']+")"

    with open( os.path.join(out_dir, "phrokka.tbl"), 'w') as f:
        for index, row in length_df.iterrows():
            contig = row['contig']
            f.write('>' + contig + '\n')
            subset_df = phanotate_mmseqs_df[phanotate_mmseqs_df['contig'] == contig]
            for index, row in subset_df.iterrows():
                f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "PHANOTATE\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "phrog=" + str(row['phrog']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['annot']) + "\n")
                f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")
            if empty == False:
                subset_trna_df = trna_df[trna_df['contig'] == contig]
                for index, row in subset_trna_df.iterrows():
                    f.write(str(row['start']) + "\t" + str(row['stop']) + "\t" + row['Region'] + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"inference" + "\t"+ "tRNAscan-SE")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"product" + "\t"+ str(row['trna_product']) + "\n")
                    f.write(""+"\t"+""+"\t"+""+"\t"+"transl_table" + "\t"+ "11" + "\n")







