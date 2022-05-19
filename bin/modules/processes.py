import os
import sys
import subprocess as sp


def run_minimap2(input_file, out_dir, db_dir, prefix):
    print("Beginning Minimap")
    try:
        # minimap2 -ax splice IMGT/rat_IMGT+C.fa test.fasta | samtools view -S -b > sample.bam
        p1 = sp.Popen(["minimap2", "-ax", "splice", os.path.join(db_dir, "rat_IMGT+C.fa"), input_file],stdout=sp.PIPE)
        fout = open(os.path.join(out_dir,prefix+ '_sample.bam'), 'wb')
        sp.run(['samtools', 'view', '-S', '-b'], stdin=p1.stdout, stdout=fout)
    except:
        sys.exit("Error with mapping step \n")  


def keep_primary_supplementary_mappings_convert_sam(out_dir,prefix):
    try:
        # samtools view -b -f 16   sample.bam | samtools view -b -F 256 | samtools view -h > mapped.sam
        p1 = sp.Popen(["samtools", "view", "-b", "-f", "16", os.path.join(out_dir, prefix+"_sample.bam")],stdout=sp.PIPE)
        p2 = sp.Popen(['samtools', 'view', '-b', '-F', "256"],stdin=p1.stdout, stdout=sp.PIPE)
        fout = open(os.path.join(out_dir,prefix+'_mapped.sam'), 'wb')
        sp.run(['samtools', 'view', '-h'], stdin=p2.stdout, stdout=fout)
    except:
        sys.exit("Error with filtering step \n")  


def run_flagstat(out_dir,prefix):
    try:
        # samtools flagstat  sample.bam > sample_flagstat.txt
        fout = open(os.path.join(out_dir,prefix+'_flagstat.txt'), 'wb')
        sp.run(["samtools", "flagstat",os.path.join(out_dir, prefix+"_sample.bam")], stdout=fout)
    except:
        sys.exit("Error with flagstat step \n")  




