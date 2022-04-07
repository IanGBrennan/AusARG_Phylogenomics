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
                help="Process alignments in parallel using --CPU for alignment. " +
                "This is the number of PHYSICAL CPUs."
        )

        parser.add_argument(
                "--b1",
                type=float,
                default=0.5,
                help="GBLOCKS -b1 proportion; min # of seqs" +
                 " required for a conserved position"
        )

        parser.add_argument(
                "--b2",
                type=float,
                default=0.85,
                help="GBLOCKS -b2 proportion; min # of seqs " +
                 " required to be at a flanking position"
        )

        parser.add_argument(
                "--b3",
                type=int,
                default=8,
                help="GBLOCKS -b3 proportion; max number of" +
                         " contiguous nonconserved positions"
        )

        parser.add_argument(
                "--b4",
                type=int,
                default=10,
                help="GBLOCKS -b4 proportion;" +
                         " minimum block length"
                )


        return parser.parse_args()

def get_dir(args):
        # make a new directory to hold the macse alignments
        outdir = os.path.join(args.indir, 'gb_alignments')
        finaldir = os.path.join(args.indir, 'clean_gb_alignments')

        if not os.path.isdir(outdir):
                os.mkdir(outdir)

        if not os.path.isdir(finaldir):
                os.mkdir(finaldir)


        return outdir, finaldir

def trim_align(args, outdir):

        alns = glob.glob(args.indir + '/*fasta')

        align_shell = os.path.join(outdir,"alignment_shell.txt")        
        ashell = open(align_shell, "w")
        
        for pp in alns:

                num_seq = 0
                f = open(pp, 'r')
                for l in f:
                        if re.search('>', l):
                                num_seq += 1
                f.close()
        # could be a good place to fix the '_R_' problem

                b1 = int(round(args.b1 * num_seq)) + 1
                b2 = int(round(args.b2 * num_seq))

                if b2 < b1:
                        b2 = b1

                sbpca = 'singularity exec gblocks_0.91b--h9ee0642_2.sif Gblocks %s -t=DNA -b1=%s -b2=%s -b3=%s -b4=%s -b5=h -p=n' % (pp, b1, b2, args.b3, args.b4)
                ashell.writelines(sbpca + "\n")
        ashell.close()
#        if args.nohup:
#                align_em = subprocess.call('nohup parallel -j %s --bar :::: %s/alignment_shell.txt > aligns_nohup.out &' % (args.CPU, outdir), shell=True)

#        else :
        align_em = subprocess.call('parallel -j %s --bar :::: %s/alignment_shell.txt' % (args.CPU, outdir), shell=True)

def moveGB(args, outdir):
        subprocess.call('mv %s/*-gb %s' % (args.indir, outdir), shell=True) 

def drop_baddies(args, outdir, finaldir):

        subprocess.call("sed -r -i 's/\s+//g' %s/*-gb" % (outdir), shell=True)

        drop_shell = os.path.join(finaldir,"dropbadseq_shell.txt")
        dshell = open(drop_shell, "w")

        al_files = glob.glob(outdir + '/*-gb')

        for z in al_files:
                sbpcd = 'singularity exec ./bbmap_38.90--he522d1c_3.sif reformat.sh in=%s out=%s.fasta minconsecutivebases=100 dotdashxton=true fastawrap=32000' % (z,z)
                dshell.writelines(sbpcd + "\n")
        dshell.close()

        dropbaddies_run = subprocess.call('parallel -j %s --bar :::: %s/dropbadseq_shell.txt' % (args.CPU, finaldir), shell=True)

def moveCLEAN(args, outdir, finaldir):
        subprocess.call('mv %s/*-gb.fasta %s' % (outdir, finaldir), shell=True) 

def main():
        # get arguments
        args = get_args()
        # get output directory
        outdir, finaldir = get_dir(args)
        # run gblocks
        trim_align(args, outdir)
        # move the *gb files
        moveGB(args, outdir)
        # drop bad sequences
        drop_baddies(args, outdir, finaldir)
        # move the final seqs
        moveCLEAN(args, outdir, finaldir)

if __name__ == "__main__":
        main()

