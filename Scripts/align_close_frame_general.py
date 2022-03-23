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
                "--indir",
                type=str,
                default=None,
                help="The directory with the alignments, if "
                         "not using in context of a pipeline."
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
        basedir = args.indir
        outdir = os.path.join(basedir, 'macse_NT')
	AAdir = os.path.join(basedir, 'macse_AA')

        if not os.path.isdir(outdir):
                os.mkdir(outdir)
	
	if not os.path.isdir(AAdir):
		os.mkdir(AAdir)


        return basedir, outdir, AAdir

def shell_macse(outdir, args, basedir):
        
        align_shell = os.path.join(outdir, "macse_shell.txt")
        ashell = open(align_shell, "w")
        
        al_files = glob.glob(basedir + '*.fasta')

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

def cleanup(args):
        subprocess.call('perl -pi -w -e "s/!/N/g;" %s/*NT.fasta' % args.indir, shell=True)

def moveNT(args, outdir, AAdir):
        subprocess.call('mv %s/*NT.fasta %s' % (args.indir, outdir), shell=True) 
	subprocess.call('mv %s/*AA.fasta %s' % (args.indir, AAdir), shell=True)

def main():
        # get arguments
        args = get_args()
        # get output directory
        basedir, outdir, AAdir = get_dir(args)
        # run macse
        #run_macse(outdir, args, basedir)
        # run shell_macse
        shell_macse(outdir, args, basedir)
        # replace any ! with N
        cleanup(args)
        # mmove the NT alignments to the new directory
        moveNT(args, outdir, AAdir)

if __name__ == "__main__":
        main()

