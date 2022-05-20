import os
import pandas as pd
import pysam

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
    imported = pysam.AlignmentFile(sam_file)
    bam_it = imported.fetch(until_eof = True)
    # Use head(n) instead of fetch(), if you only want to retrieve the first 'n' reads
    qnames = []
    flags = []
    rnames = []

    for read in bam_it:
        qnames.append(read.query_name)
        flags.append(read.flag)
        rnames.append(read.reference_name)

    match_df = pd.DataFrame(list(zip(qnames, flags, rnames)),
               columns =['qname', 'flag', 'rname'])
    return(match_df)

def pivot_df(out_dir, df, total_read_count, prefix):
    # get counts for each igg
    primary_align_df = df[df["flag"] == 16] 
    # remove chromosome (might address this in database)
    primary_align_df[['rname','chr']] = primary_align_df['rname'].str.split('_',expand=True)
    counts = primary_align_df['rname'].value_counts()
    tpm = counts/int(total_read_count)*1000000
    tpm_df =tpm.to_frame()
    tpm_df = tpm_df.rename_axis('IG').reset_index()
    tpm_df = tpm_df[tpm_df['IG'].str.contains("IGH")]
    #print(tpm_df)
    # summing 
    IGA = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA"), 'rname'].sum(),1)
    IGG = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG"), 'rname'].sum(),1)
    IGE = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHE"), 'rname'].sum(),1)
    IGM = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHM"), 'rname'].sum(),1)
    IGD = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHD"), 'rname'].sum(),1)
    ig_dict = {IGA: 'igA', IGG: 'igG', IGE: 'igE', IGM: 'igM', IGD: 'igD'}
    # convert to dataframe
    ig_items = ig_dict.items()
    ig_list = list(ig_items)
    ig_df = pd.DataFrame(ig_list) 
    ig_df.columns = ['TPM','IG']
    cols = ig_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    ig_df = ig_df[cols]
    print(ig_df)
    ig_df.to_csv( os.path.join(out_dir, prefix + "_ig_summary.csv"), sep=",", index=False)

    ### subsets 
    # https://stackoverflow.com/questions/26577516/how-to-test-if-a-string-contains-one-of-the-substrings-in-a-list-in-pandas
    #IGHA
    IGA1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA\*01"), 'rname'].sum(),1)
    IGA2 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA\*02"), 'rname'].sum(),1)
    IGG1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG1"), 'rname'].sum(),1)
    IGG2A = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2A"), 'rname'].sum(),1)
    IGG2B = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2B"), 'rname'].sum(),1)
    IGG2C = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2C"), 'rname'].sum(),1)
    iga_dict = {IGA1: 'igA1', IGA2: 'igA2',IGG1: 'igG1',   IGG2A: 'igG2A',  IGG2B: 'igG2B',  IGG2C: 'igG2C'}
    # convert to dataframe
    iga_items = iga_dict.items()
    iga_list = list(iga_items)
    iga_df = pd.DataFrame(iga_list) 
    iga_df.columns = ['TPM','IG']
    cols = iga_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    iga_df = iga_df[cols]
    print(iga_df)
    iga_df.to_csv( os.path.join(out_dir, prefix + "_iga_subtypes_summary.csv"), sep=",", index=False)

  




