import argparse
import glob
import os
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 23 June 2016
Written assuming:
	* mafft 7.294
	* RAxML 8.2.4

Adapted for Singularity by Ian Brennan
3 July 2021

This script borrows heavily from:
https://github.com/faircloth-lab/phyluce/blob/master/bin/align/phyluce_align_seqcap_align
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
		"--outdir",
		type=str,
		default=None,
		help="The directory with the phylogeny, if "
			 "not using in context of a pipeline."
	)

	parser.add_argument(
		"--raxmltrees",
		action="store_true",
		default=False,
		help="Will infer RAxML gene trees for each alignment if "
		 "flagged."
	)

	parser.add_argument(
		"--phymltrees",
		action="store_true",
		default=False,
		help="Will infer PhyML gene trees for each alignment if "
			 "flagged."
	)

	parser.add_argument(
		"--iqtrees",
		action="store_true",
		default=False,
		help="Will infer gene trees with IQTREE for each alignment "
			"flagged."
	)

	parser.add_argument(
		"--trim",
		action="store_true",
		default=False,
		help="Will trim alignments using gblocks if flagged."
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

	parser.add_argument(
		"--CPU",
		type=int,
		default=1,
		help="""Process alignments in parallel using --CPU for alignment. """ +
		"""This is the number of PHYSICAL CPUs."""
		)

	#parser.add_argument(
	#	"--mafft",
	#	type=str,
	#	default=None,
	#	help="Full path to mafft executable."
	#	)

	#parser.add_argument(
	#	"--gblocks",
	#	type=str,
	#	default=None,
	#	help="Full path to GBLOCKS executable."
	#	)

	parser.add_argument(
		"--raxml",
		type=str,
		default=None,
		help="Full path to RAxML executable."
		)

	parser.add_argument(
		"--tree_method",
		type=str,
		default='iqtree',
		help="Your choice of phylogenetic reconstruction program. " +
				"Options are: iqtree, raxml, phyml."
		)

	parser.add_argument(
		"--jmodel",
		type=str,
		default=None,
		help="Full path to jModelTest jar."
		)

	return parser.parse_args()


def get_dir(args):
	if not args.outdir:
		outdir = os.path.join(args.dir, 'phylogeny', 'alignments')
		treedir = os.path.join(args.dir, 'phylogeny', 'gene_trees')
	else:
		outdir = os.path.join(args.outdir, 'alignments')
		treedir = os.path.join(args.outdir, 'gene_trees')	
	
	return outdir, treedir


def align(params):
	file, mafft = params

	aln_out = file.replace('.fasta', '.fasta.aln')
	proc = subprocess.call("singularity exec mafft_7.480--h779adbc_0.sif mafft --maxiterate 1000 --globalpair "
			"--adjustdirection --quiet %s > %s" %
			(file, aln_out), shell=True)

	os.remove(file)
	return aln_out


def run_alignments(outdir, args):
	files = glob.glob(outdir + '/*fasta')
	
	if len(files) > 0:
		params = zip(files, ['singularity exec mafft_7.480--h779adbc_0.sif mafft'] * len(files))
	
		if args.CPU > 1:
			pool = mp.Pool(args.CPU)
			alns = pool.map(align, params)
	alns = glob.glob(outdir + '/*fasta.aln')	
		
	return alns


def trim_align(args, alns, outdir):
	#(aln, gblocks, b1, b2, b3, b4) = params

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

	align_em = subprocess.call('parallel -j %s --bar :::: %s/alignment_shell.txt' % (args.CPU, outdir), shell=True)

# let's remove samples that have little sequence data in a given alignment
def drop_baddies(args, outdir):
	# Gblocks outputs alignments as 10bp chunks separated by a space (and line wrapped at 60bp)
	# use sed to remove all those spaces, because BBmap won't handle them
	if args.trim:
		subprocess.call("sed -r -i 's/\s+//g' %s/*fasta.aln-gb" % (outdir), shell=True)
	else:
		subprocess.call("sed -r -i 's/\s+//g' %s/*fasta.aln" % (outdir), shell=True)

	drop_shell = os.path.join(outdir,"dropbadseq_shell.txt")
	dshell = open(drop_shell, "w")

	if args.trim:
		al_files = glob.glob(outdir + '/*fasta.aln-gb')
	else:
		al_files = glob.glob(outdir + '/*fasta.aln')
	
	# use reformat.sh from BBmap to remove sequences without at least 100 consecutive basepairs of data
	# we have to append ".fasta" to the end or BBmap will output the file as a fastq
	for z in al_files:
		sbpcd = 'singularity exec ./bbmap_38.90--he522d1c_3.sif reformat.sh in=%s out=%s.fasta minconsecutivebases=100 dotdashxton=true fastawrap=32000' % (outdir,z,outdir,z)
         	dshell.writelines(sbpcd + "\n")
	dshell.close()

	dropbaddies_run = subprocess.call('parallel -j %s --bar :::: %s/dropbadseq_shell.txt' % (args.CPU, outdir), shell=True)


# remove the prepended _R_ from reversed seqs
def fix_reverse(args, outdir):
	if args.trim:
		subprocess.call("grep --include=\*.fasta.aln-gb.fasta -rl '%s' -e '_R_' >> %s/R_files.txt" % (outdir, outdir), shell=True)
	else:
		subprocess.call("grep --include=\*.fasta.aln.fasta -rl '%s' -e '_R_' >> %s/R_files.txt" % (outdir, outdir), shell=True)

# Need to think of a way to read in the file and remove those things
	rgb = open('%s/R_files.txt' % outdir, 'r')
	r_aln = rgb.readlines()
	r_aln = [x.strip() for x in r_aln]	

	#for zz in r_aln.tolist():
	#	print(zz)	

	for kk in r_aln:
		subprocess.call("sed -i 's/_R_//g' %s" % (kk), shell=True)

#def run_trimming(alns, args, outdir):
	 
	params = []
 #       align_shell = os.path.join(outdir,"genetree_shell.txt")
 #       ashell = open(align_shell, "w")

	#al_files = glob.glob(outdir + '/*fasta.aln')
	
	
#	for k in al_files:
#		sbpc_al = '

#	for aln in alns:
#		param = [aln, args.gblocks, args.b1, args.b2, args.b3, args.b4]
#		params.append(param)
		
		#if args.CPU > 1:
				#pool = mp.Pool(args.CPU)
				#trim = pool.map(trim_align, params)
		#		print('trimming using parallel')
	#trim = subprocess.call('parallel -j %s --bar :::: params' % args.CPU, shell=True)
		#return trim


def convert_phyml(locus_file):
	f = open(locus_file, 'r')
	phy_file = re.sub('.fasta.*$', '.aln.phy', locus_file)
	o = open(phy_file, 'w')

	seq = {}
	id = ''
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			if re.search('^_R_', id):
				id = re.sub('^_R_', '', id)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	for sp, s in seq.items():
		# get rid of white spaces from gblocks
		s = re.sub('\s+', '', s)
		seq[sp] = s

	o.write(' %s %s\n' % (len(seq), len(seq.values()[0])))
	for sp, s in seq.items():
		o.write('%s   %s\n' % (sp, s))
	o.close()

	return phy_file

def sub_raxml(file, outdir, raxml):

	locus = re.sub('^.*/', '', file)
	locus = re.sub('\.aln.*', '', locus)

	os.chdir(outdir)
	
	orig_boot = 'RAxML_bootstrap.%s' % locus
	orig_tree = 'RAxML_bipartitions.%s' % locus

	new_boot = '%s.bootstrap.trees' % locus
	new_tree = '%s.bestTree.tre' % locus

	if not os.path.isfile(new_tree):
		subprocess.call('%s -x %s -# 100 -p %s -m GTRCAT -f a -n %s -s %s' % 
							(raxml, random.randint(0,1000), random.randint(0,1000), 
							locus, file), shell=True)

		os.rename(orig_boot, new_boot)
		os.rename(orig_tree, new_tree)

		subprocess.call("rm RAxML_*%s" % locus, shell=True) 

	return new_tree

def shell_iqtree(file, outdir, treedir, args):

	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	tree_shell = os.path.join(treedir,"genetree_shell.txt")
	tshell = open(tree_shell, "w")

	if args.trim:
		al_files = glob.glob(outdir + '/*fasta.aln-gb.fasta')
	else:
		al_files = glob.glob(outdir + '/*fasta.aln.fasta')
	
	for z in al_files:
		sbpc = 'singularity exec /home/ian/SqCL_Pipeline/iqtree_2.1.2--h56fc30b_0.sif iqtree -s %s --quiet -T AUTO -B 1000' % z
         	tshell.writelines(sbpc + "\n")
	tshell.close()

	iqtree_run = subprocess.call('parallel -j %s --bar :::: %s/genetree_shell.txt' % (args.CPU, treedir), shell=True)	

result_list = []
def log_result(result):
	# This is called whenever foo_pool(i) returns a result.
	# result_list is modified only by the main process, not the pool workers.
	result_list.append(result)


def run_raxml(outdir, treedir, alns, args):	
	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	if args.CPU > 1:
		pool = mp.Pool(args.CPU)
		phys = pool.map(convert_phyml, alns)

		dirs = [treedir] * len(phys)
		raxml = [args.raxml] * len(phys)
		for i in range(len(phys)):
			pool.apply_async(sub_raxml, args=(phys[i], treedir, args.raxml, ), callback=log_result)
		pool.close()
		pool.join()

def run_iqtree(outdir, treedir, alns, args):
	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	if args.CPU > 1:
		pool = mp.Pool(args.CPU)
		phys = pool.map(convert_phyml, alns)

		for i in range(len(phys)):
			pool.apply_async(sub_iqtree, args=(phys[i], treedir, ), callback=log_result)
		pool.close()
		pool.join()


def sub_phyml(file, outdir, jmodel):

	os.chdir(outdir)

	locus = re.sub('^.*/', '', file)
	locus = re.sub('\.aln.*', '', locus)

	out1 = os.path.join(outdir, '%s.jmodel.txt' % locus)
	out2 = os.path.join(outdir, '%s.jmodel.tre' % locus)

	subprocess.call('java -jar %s -d %s -g 4 -i -f -AIC -dLRT -o %s' %
		(jmodel, file, out1), shell=True)

	f = open(out1, 'r')
	o = open(out2, 'w')
	for l in f:
		if re.search('Tree for the best AIC model', l):
			tree = re.search('=\s+(.*)$', l).group(1)
			break
	f.close()
	o.write(tree)
	o.close()

	return out2


def run_phyml(outdir, treedir, alns, args):
	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	if args.CPU > 1:
		pool = mp.Pool(args.CPU)
		phys = pool.map(convert_phyml, alns)

		for i in range(len(phys)):
				pool.apply_async(sub_phyml, args=(phys[i], treedir, args.jmodel, ), callback=log_result)
		pool.close()
		pool.join()

#def move_trees(outdir, treedir):
#	subprocess.call('mv %s/*.contree %s' % (outdir, treedir), shell=True)
	

def main():
	args = get_args()
	outdir, treedir = get_dir(args)	
	alns = run_alignments(outdir, args)

	tree_alns = alns
	if args.trim:
		trim_align(args, alns, outdir)
		drop_baddies(args, outdir)
		fix_reverse(args, outdir)
		#trims = run_trimming(alns, args, outdir)
		#tree_alns = trims
	if args.tree_method == 'raxml':
		run_raxml(outdir, treedir, tree_alns, args)
	if args.tree_method == 'phyml':
		run_phyml(outdir, treedir, tree_alns, args)
	if args.tree_method == 'iqtree':
		#run_iqtree(outdir, treedir, tree_alns, args)
		shell_iqtree(file, outdir, treedir, args)
	if args.tree_method == 'none':
		print('** NOTE: we have not built trees **')
	#clean_up = move_trees(outdir, treedir)

if __name__ == "__main__":
	main()
