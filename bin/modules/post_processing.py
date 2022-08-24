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

def pivot_df(out_dir, df, total_read_count, prefix, kit, species):
    # get counts for each igg
    # cdna keeps all the reverse strand mapping
    if kit == "cdna":
        numbers = [16, 2048]
    # direct
    else:
        numbers = [0, 2048]
    primary_align_df = df[df["flag"].isin(numbers)] 
    # remove chromosome (might address this in database)
    print(primary_align_df['rname'].unique())
    primary_align_df[['rname','chr']] = primary_align_df['rname'].str.split('_',expand=True)
    counts = primary_align_df['rname'].value_counts()
    tpm = counts/int(total_read_count)*1000000
    tpm_df =tpm.to_frame()
    tpm_df = tpm_df.rename_axis('IG').reset_index()
    tpmh_df = tpm_df[tpm_df['IG'].str.contains("IGH")]
    tpml_df = tpm_df[tpm_df['IG'].str.contains("IGL")]
    tpmk_df = tpm_df[tpm_df['IG'].str.contains("IGK")]
    #print(tpm_df)
    # summing 
    IGA = round(tpmh_df.loc[tpmh_df['IG'].str.contains("IGHA"), 'rname'].sum(),1)
    IGG = round(tpmh_df.loc[tpmh_df['IG'].str.contains("IGHG"), 'rname'].sum(),1)
    IGE = round(tpmh_df.loc[tpmh_df['IG'].str.contains("IGHE"), 'rname'].sum(),1)
    IGM = round(tpmh_df.loc[tpmh_df['IG'].str.contains("IGHM"), 'rname'].sum(),1)
    IGD = round(tpmh_df.loc[tpmh_df['IG'].str.contains("IGHD"), 'rname'].sum(),1)
    IGKC = round(tpmk_df.loc[tpmk_df['IG'].str.contains("IGKC"), 'rname'].sum(),1)
    IGLC = round(tpml_df.loc[tpml_df['IG'].str.contains("IGLC"), 'rname'].sum(),1)
    ig_dict = {IGA: 'igA', IGG: 'igG', IGE: 'igE', IGM: 'igM', IGD: 'igD', IGKC: 'igK', IGLC: 'igL'}
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


    # rats

    if species == "rat":
    #IGHA
        IGA1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA\*01"), 'rname'].sum(),1)
        IGA2 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA\*02"), 'rname'].sum(),1)
        IGG1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG1"), 'rname'].sum(),1)
        IGG2 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2"), 'rname'].sum(),1)
        IGG2A = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2A"), 'rname'].sum(),1)
        IGG2B = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2B"), 'rname'].sum(),1)
        IGG2C = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2C"), 'rname'].sum(),1)
        ig_sub_dict = {IGA1: 'igA1', IGA2: 'igA2',IGG1: 'igG1', IGG2:'igG2',   IGG2A: 'igG2A',  IGG2B: 'igG2B',  IGG2C: 'igG2C'} 
    elif species == "human":
        IGA1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA1"), 'rname'].sum(),1)
        IGA2 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHA2"), 'rname'].sum(),1)
        IGG1 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG1"), 'rname'].sum(),1)
        IGG2 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG2"), 'rname'].sum(),1)
        IGG3 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG3"), 'rname'].sum(),1)
        IGG4 = round(tpm_df.loc[tpm_df['IG'].str.contains("IGHG4"), 'rname'].sum(),1)
        ig_sub_dict = {IGA1: 'igA1', IGA2: 'igA2',IGG1: 'igG1', IGG2:'igG2', IGG3:'igG3', IGG4: 'igG4'}
    
    # convert to dataframe
    ig_sub_items = ig_sub_dict.items()
    iga_sub_list = list(ig_sub_items)
    ig_sub_df = pd.DataFrame(iga_sub_list) 
    ig_sub_df.columns = ['TPM','IG']
    cols = ig_sub_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    ig_sub_df = ig_sub_df[cols]
    print(ig_sub_df)
    ig_sub_df.to_csv( os.path.join(out_dir, prefix + "_ig_subtypes_summary.csv"), sep=",", index=False)

  




