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
                "--indir",
                type=str,
                default=None,
                help="Full path to directory with the *.fasta alignments. " +
                "e.g. if wd is SqCL_Pipeline: --indir Frogs/Alignments"
        )

        parser.add_argument(
                "--CPU",
                type=int,
                default=1,
                help="Generate genetrees in parallel using --CPU for alignment. " +
                "This is the number of PHYSICAL CPUs."
        )

        return parser.parse_args()

def get_dir(args):
        # make a new directory to hold the macse alignments
        basedir = args.indir
        outdir = os.path.join(basedir, 'contrees')
	IQdir = os.path.join(basedir, 'iqtree_files')

        if not os.path.isdir(outdir):
                os.mkdir(outdir)
	
	if not os.path.isdir(IQdir):
		os.mkdir(IQdir)


        return basedir, outdir, IQdir

def shell_iqtree(outdir, IQdir, args):

        tree_shell = os.path.join(outdir,"genetree_shell.txt")
        tshell = open(tree_shell, "w")

        al_files = glob.glob(args.indir + '/*.fasta')
        
        for z in al_files:
                sbpc = 'singularity exec ./iqtree_2.1.2--h56fc30b_0.sif iqtree -s %s --quiet -T AUTO -B 1000' % z
                tshell.writelines(sbpc + "\n")
        tshell.close()

#        if args.nohup:
#                iqtree_run = subprocess.call('nohup parallel -j %s --bar :::: %s/genetree_shell.txt > genetrees_nohup.out &' % (args.CPU, treedir), shell=True)       

#        else:
        iqtree_run = subprocess.call('parallel -j %s --bar :::: %s/genetree_shell.txt' % (args.CPU, outdir), shell=True)       


def moveTREE(args, outdir, IQdir):
        subprocess.call('mv %s/*contree %s' % (args.indir, outdir), shell=True) 
	subprocess.call('mv %s/*.model.gz %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.log %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.ckp.gz %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.mldist %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.bionj %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.splits.nex %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.treefile %s' % (args.indir, IQdir), shell=True)
        subprocess.call('mv %s/*.iqtree %s' % (args.indir, IQdir), shell=True)


def main():
        # get arguments
        args = get_args()
        # get output directory
        basedir, outdir, IQdir = get_dir(args)
        # run macse
        #run_macse(outdir, args, basedir)
        # run shell_macse
        shell_iqtree(outdir, IQdir, args)
        # move the contree alignments to the new directory; and all other iqtree files to IQtree dir
        moveTREE(args, outdir, IQdir)

if __name__ == "__main__":
        main()

