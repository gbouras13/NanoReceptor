#!/usr/bin/env python3
import os
import sys
import subprocess as sp
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def instantiate_install(db_dir):
    instantiate_dir(db_dir)
    get_IMGT_fa(db_dir)
    get_rat_receptors(db_dir)

# create the dir

def instantiate_dir(db_dir):
    if os.path.isdir(db_dir) == False:
        os.mkdir(db_dir)
 
# gets the IMGT.fa file 
# wget "http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP" -P . -O total_fasta.fa

def get_IMGT_fa(db_dir):
    print("Getting IMGT Database")
    filepath = "http://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP"
    outfile = os.path.join(db_dir,"total_fasta.fa")
    
    # get tarball if not already present
    if os.path.isfile(os.path.join(db_dir,"total_fasta.fa")) == True: 
         print("IMGT Database already downloaded")
        # download 
    else:
        try:
            sp.call(["wget", filepath, "-P", db_dir, "-O", outfile])
        except:
            sys.stderr.write("Error: IMGT Database not found - link likely broken\n")  
            return 0

def get_rat_receptors(db_dir):
    print("Getting Rat IMGT File")

    # getting entries with multiple entries

    headers = []
    
    for dna_record in SeqIO.parse(os.path.join(db_dir, "total_fasta.fa"), "fasta"):
        # split the list
        spl_list = dna_record.id.split("|")
        # get second item (the receptor)
        header = spl_list[1]
        # get third item - species
        species = spl_list[2]
        # write new record if has 'Rattus'
        if 'Rattus' in species:
            headers.append(header)

    # get all dupes     

    seen = set()
    dupes = []

    for x in headers:
        if x in seen:
            dupes.append(x)
        else:
            seen.add(x)


    with open(os.path.join(db_dir, "rat_IMGT+C.fa"), 'w') as fa_out:
        # need description for entire header
        for dna_record in SeqIO.parse(os.path.join(db_dir, "total_fasta.fa"), "fasta"):
            # split the list
            spl_list = dna_record.description.split("|")
            # split to get after g,
            # get second item (the receptor)
            header = spl_list[1]
            # add append to id if a dupe - choose the chromosome (I think)
            if header in dupes:
                header = header + "_" + spl_list[4]
            # get third item - species
            species = spl_list[2]
            # write new record if has 'Rattus norvegicus'
            if 'Rattus norvegicus' in species:
                fa_record = SeqRecord(dna_record.seq, id=header, description="")
                SeqIO.write(fa_record, fa_out, 'fasta')




