import argparse
import glob
import os
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import random

"""
Ian Brennan
created 23 March 2022
Written assuming:
	* macse XX
"""

def get_args():
        parser = argparse.ArgumentParser(
                description="Align, possibly trim, and infer gene trees for UCE loci.",
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )

        parser.add_argument(
                "--dir",
                type=str,
                default=None,
                help="The base directory if running the pipeline."
        )

	parser.add_argument(
		"--program",
		type=str,
		default=None,
		help="Settings for aligning with MACSE. " +
			"Options are: refine, refineLemmon, align, export."
	)

        parser.add_argument(
                "--outdir",
                type=str,
                default=None,
                help="The directory with the phylogeny, if "
                         "not using in context of a pipeline."
        )

	parser.add_argument(
		"--trimmed",
		action="store_true",
		default=False,
		help="Were the input alignments trimmed with Gblocks?"
	)

        parser.add_argument(
                "--CPU",
                type=int,
                default=1,
                help="""Process alignments in parallel using --CPU for alignment. """ +
                """This is the number of PHYSICAL CPUs."""
        )

	return parser.parse_args()

def get_dir(args):
	# make a new directory to hold the macse alignments
	basedir = os.path.join(args.dir, 'phylogeny', 'alignments')
	outdir = os.path.join(basedir, 'macse')
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

	# copy all the AHE alignments to the macse directory
	if args.trimmed:
		subprocess.call('cp %s/AHE-L*-gb.fasta %s' % (basedir, outdir), shell=True)
	else:
                subprocess.call('cp %s/AHE-L*-fasta.aln.fasta %s' % (basedir, outdir), shell=True)

	return basedir, outdir

def macse(params):
	# set the arguments for running macse
	file, macse = params
	
	proc = subprocess.call("java -jar ./macse_v2.06.jar -prog refineAlignment -align %s -stop 10" % file, shell=True)
#	if args.program == refine:
#		proc = subprocess.call("java -jar ./macse_v2.06.jar -prog refineAlignment -align %s -stop 10" % file, shell=True)
#	if args.program == refineLemmon:
#		proc = subprocess.call("java -jar ./macse_v2.06.jar -prog refineAlignment -align %s -optim 1 -local_realign_init 0.1 -local_realign_dec 0.1 -fs 10" % file, shell=True)
#	if args.program == export:
#		proc == subprocess.call("java -jar ./macse_v2.06.jar -prog exportAlignment -align %s -stop 10" % file, shell=True)
#	if args.program == align:
#		proc == subprocess.call("java -jar ./macse_v2.06.jar -prog alignSequences -seq_lr %s -stop_lr 10" % file, shell=True)
	
	os.remove(file)


def run_macse(outdir, args, basedir):
	# get the list of files to realign
	#if args.trimmed:
	#	files = glob.glob(outdir + '/AHE-L*-gb.fasta')
	#else:
	#	files = glob.glob(outdir + '/AHE-L*.fasta.aln.fasta')
	
	files = glob.glob(outdir + '/AHE-L*-gb.fasta')	

	if len(files) > 0:
		params = zip(files, ['java -jar ./macse_v2.06.jar'] * len(files))
		
		if args.CPU > 1:
			pool = mp.Pool(args.CPU)
			alns = pool.map(macse, params)
	return alns

def shell_macse(outdir, args, basedir):
	
	align_shell = os.path.join(outdir, "macse_shell.txt")
	ashell = open(align_shell, "w")
	
        if args.trimmed:
                al_files = glob.glob(outdir + '/*fasta.aln-gb.fasta')
        else:
                al_files = glob.glob(outdir + '/*fasta.aln.fasta')

	# run a loop across all the files and add them to the shell script
	for z in al_files:
		if args.program == 'refine':
               		proc = 'java -jar ./macse_v2.06.jar -prog refineAlignment -align %s -stop 10' % z
			ashell.writelines(proc + "\n")
		if args.program == 'refineLemmon':
	                proc = "java -jar ./macse_v2.06.jar -prog refineAlignment -align %s -optim 1 -local_realign_init 0.1 -local_realign_dec 0.1 -fs 10" % z
                        ashell.writelines(proc + "\n")
	        if args.program == 'export':
	                proc == "java -jar ./macse_v2.06.jar -prog exportAlignment -align %s -stop 10" % z
                        ashell.writelines(proc + "\n")
	        if args.program == 'align':
	                proc == "java -jar ./macse_v2.06.jar -prog alignSequences -seq_lr %s -stop_lr 10" % z
                        ashell.writelines(proc + "\n")
	ashell.close()

	macse_run = subprocess.call('parallel -j %s --bar :::: %s/macse_shell.txt' % (args.CPU, outdir), shell=True)

def cleanup(outdir):
	subprocess.call('perl -pi -w -e "s/!/N/g;" %s/*NT.fasta' % outdir, shell=True)

def main():
	# get arguments
	args = get_args()
	# get output directory
	basedir, outdir = get_dir(args)
	# run macse
	#run_macse(outdir, args, basedir)
	# run shell_macse
	shell_macse(outdir, args, basedir)
	# replace any ! with N
	cleanup(outdir)
if __name__ == "__main__":
	main()
