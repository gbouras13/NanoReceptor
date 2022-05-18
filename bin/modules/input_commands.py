import argparse
import os
import sys
from argparse import RawTextHelpFormatter
import datetime
from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial

v = '1.0.0'


### GLOBAL VARIABLES


def get_input():
	usage = 'NanoReceptor ...'
	parser = argparse.ArgumentParser(description='NanoReceptor: Program to infer IG and TRA quantities from Long Read RNA-Seq Data', formatter_class=RawTextHelpFormatter)
	parser.add_argument('-i', '--infile', action="store", help='input file in fastq format',  required=True)
	parser.add_argument('-o', '--outdir', action="store", help='directory to write the output', default=os.path.join(os.getcwd(), "output/") )
	parser.add_argument('-t', '--threads', action="store", help='number of threads', default=1 )
	parser.add_argument('-V', '--version', action='version', version=v)
	args = parser.parse_args()

	return args

def instantiate_dirs(output_dir):
	# if the output directory already exists add the date and time on the end to make unique dir
	# comment out if running tests
	# if os.path.isdir(output_dir) == True:
	# 	output_dir = str(output_dir).rstrip('/') + str(datetime.datetime.now().strftime('_%Y%m%d_%H%M_%S%f'))
	if os.path.isdir(output_dir) == False:
		os.mkdir(output_dir)
	return output_dir

# https://stackoverflow.com/questions/42757283/seqio-parse-on-a-fasta-gz

def validate_fastq(filename):
	encoding = guess_type(filename)[1]  # uses file extension
	_open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open
	with _open(filename) as f:
		fastq = SeqIO.parse(f, "fastq")
		if any(fastq):
			print("FASTQ checked")
		else:
			sys.exit("Error: Input file is not in the FASTQ format.\n")  





