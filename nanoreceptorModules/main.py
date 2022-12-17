from audioop import ratecv
import os
import sys
# from input_commands import get_input, validate_fastq, instantiate_dirs
# from processes import run_minimap2, keep_primary_supplementary_mappings_convert_sam, keep_primary_supplementary_mappings_convert_sam_direct, run_flagstat
# from post_processing import get_total_read_count, parse_bam, pivot_df


DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

def main(argv):  

    args = get_input()
    validate_fastq(args.infile)
    out_dir = instantiate_dirs(args.outdir) # incase there is already an outdir
    run_minimap2(args.infile, out_dir, DBDIR, args.prefix, args.species)
    if args.kit == "cdna":
        keep_primary_supplementary_mappings_convert_sam(out_dir, args.prefix)
    elif args.kit == "direct":
        keep_primary_supplementary_mappings_convert_sam_direct(out_dir, args.prefix)
    run_flagstat(out_dir, args.prefix)
    # post processing
    total_read_count = get_total_read_count(out_dir,args.prefix)
    match_df = parse_bam(out_dir, args.prefix)
    pivot_df(out_dir, match_df, total_read_count, args.prefix, args.kit, args.species)
   
    sys.exit("nanoreceptor has finished")  


def run():
    main(sys.argv)


