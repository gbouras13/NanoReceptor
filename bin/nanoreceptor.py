#!/usr/bin/env python3
import sys
from modules import input_commands
from modules import processes
from modules import post_processing
import os

DBDIR = os.path.join(os.path.dirname(__file__),'../',"databases/")  

if __name__ == "__main__":
    args = input_commands.get_input()
    input_commands.validate_fastq(args.infile)
    out_dir = input_commands.instantiate_dirs(args.outdir) # incase there is already an outdir
    processes.run_minimap2(args.infile, out_dir, DBDIR, args.prefix)
    processes.keep_primary_supplementary_mappings_convert_sam(out_dir, args.prefix)
    processes.run_flagstat(out_dir, args.prefix)
    # post processing
    total_read_count = post_processing.get_total_read_count(out_dir,args.prefix)
    match_df = post_processing.parse_bam(out_dir, args.prefix)
    post_processing.pivot_df(out_dir, match_df, total_read_count, args.prefix)
   
    sys.exit("nanoreceptor has finished")  



