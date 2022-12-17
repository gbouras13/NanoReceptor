from audioop import ratecv
import os
import sys
import nanoreceptorModules

DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

def main(argv):  

    args = nanoreceptorModules.get_input()
    nanoreceptorModules.validate_fastq(args.infile)
    out_dir = nanoreceptorModules.instantiate_dirs(args.outdir) # incase there is already an outdir
    nanoreceptorModules.run_minimap2(args.infile, out_dir, DBDIR, args.prefix, args.species)
    if args.kit == "cdna":
        nanoreceptorModules.keep_primary_supplementary_mappings_convert_sam(out_dir, args.prefix)
    elif args.kit == "direct":
        nanoreceptorModules.keep_primary_supplementary_mappings_convert_sam_direct(out_dir, args.prefix)
    nanoreceptorModules.run_flagstat(out_dir, args.prefix)
    # post processing
    total_read_count = nanoreceptorModules.get_total_read_count(out_dir,args.prefix)
    match_df = nanoreceptorModules.parse_bam(out_dir, args.prefix)
    nanoreceptorModules.pivot_df(out_dir, match_df, total_read_count, args.prefix, args.kit, args.species)
   
    print("nanoreceptor has finished")


def run():
    main(sys.argv)


