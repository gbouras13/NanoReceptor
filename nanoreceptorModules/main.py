from audioop import ratecv
import os
import sys
from nanoreceptorModules import input_commands
from nanoreceptorModules import processes
from nanoreceptorModules import post_processing

DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

def main(argv):  

    args = input_commands.get_input()
    input_commands.validate_fastq(args.infile)
    out_dir = input_commands.instantiate_dirs(args.outdir) # incase there is already an outdir
    processes.run_minimap2(args.infile, out_dir, DBDIR, args.prefix, args.species)
    if args.kit == "cdna":
        processes.keep_primary_supplementary_mappings_convert_sam(out_dir, args.prefix)
    elif args.kit == "direct":
        processes.keep_primary_supplementary_mappings_convert_sam_direct(out_dir, args.prefix)
    processes.run_flagstat(out_dir, args.prefix)
    # post processing
    total_read_count = post_processing.get_total_read_count(out_dir,args.prefix)
    match_df = post_processing.parse_bam(out_dir, args.prefix)
    post_processing.pivot_df(out_dir, match_df, total_read_count, args.prefix, args.kit, args.species)
   
    sys.exit("nanoreceptor has finished")  


def run():
    main(sys.argv)


