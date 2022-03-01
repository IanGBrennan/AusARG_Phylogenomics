require(stringr)
require(dplyr)

# A series of functions to help process data generated as part of the AusARG phylogenomics project
## most functions produce a shell script that can be loaded and run on a separate machine. 
## For more information see the 'AusARG Phylogenomics Workflow' file.

# generate a sample file for downstream purposes
concatenate_collate <- function(sample.dir, metadata.file, outfile="samples.csv",
                                out.path = "/home/ian/SqCL_Pipeline/Optimization/",
                                adaptor1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG",
                                adaptor2 = "AATGATACGGCGACCACCGAGATCTACAC*ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
                                total.readno = c(2,6,8)){
  
  # set working directory to sample.dir
  setwd(sample.dir)
  
  # read in the names of the sequence files
  seq.files <- dir(sample.dir, pattern=".fastq.gz")
  
  # make a dataframe after cleaning the information
  file.ids <- data.frame(filename = seq.files,
                         sample.id = sapply(seq.files, function(x) strsplit(x, "_")[[1]][1]),
                         barcode = sapply(seq.files, function(x) strsplit(x, "_")[[1]][5]), # was [3]
                         sample.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][6]), # was [4]
                         l.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][7]), # was [5]
                         direction = sapply(seq.files, function(x) strsplit(x, "_")[[1]][8])) # was [8]
  rownames(file.ids) <- NULL
  
  # read in the sample metadata file
  info <- read.csv(paste0(sample.dir,"/",metadata.file), header=T)
  # choose the appropriate data
  info <- dplyr::select(info, "library_id", "genus", "species", "specimen_id", 
                 #"library_index_seq_P7", "library_index_seq_P5")
                 "library_index_seq", "library_index_seq_dual")
  # make a column for the sample info
  info <- dplyr::mutate(info, sample = paste0(genus, "_", species, "_", specimen_id))
  # add the adaptor information to the file
  info$adaptor1 <- adaptor1
  info$adaptor2 <- adaptor2
  
  # identify all forward files
  forwards <- filter(file.ids, direction == "R1")
  # identify all reverse files
  reverses <- filter(file.ids, direction == "R2")
  
  concatenated.files <- NULL
  for (j in 1:nrow(info)){
    # choose current sample
    curr.sample <- info$library_id[[j]]
    
    # filter the appropriate fwd/rev files
    curr.fwd <- filter(forwards, curr.sample == sample.id)
    curr.rev <- filter(reverses, curr.sample == sample.id)
    
    # make the outfile names
    out.fwd <- paste0(info[j,"sample"], "_", curr.sample, "_R1_concat.fastq.gz")
    out.rev <- paste0(info[j,"sample"], "_", curr.sample, "_R2_concat.fastq.gz")
    
    # make the bash call to concatenate the files (depends on the number of read files per sample!)
    if(total.readno == 2) {
          cat.fwd <- paste("mv", curr.fwd$filename[[1]], out.fwd)
          cat.rev <- paste("mv", curr.rev$filename[[1]], out.rev)
    }
    if(total.readno == 8) {
          cat.fwd <- paste("cat", curr.fwd$filename[[1]],
                     curr.fwd$filename[[2]],
                     curr.fwd$filename[[3]],
                     curr.fwd$filename[[4]],
                     ">>", out.fwd)
          cat.rev <- paste("cat", curr.rev$filename[[1]],
                     curr.rev$filename[[2]],
                     curr.rev$filename[[3]],
                     curr.rev$filename[[4]],
                     ">>", out.rev)
    }

    
    # do the concatenating
    system(cat.fwd)
    system(cat.rev)
    
    # add the file names to a dataframe
    concatenated.files <- rbind(concatenated.files, data.frame(read1 = paste0(out.path, out.fwd), 
                                                               read2 = paste0(out.path, out.rev)))
  }
  
  # create the SqCL_Pipeline sample file
  pipeline.in <- data.frame(sample = paste0(info$sample, "_", info$library_id),
                            read1 = concatenated.files$read1,
                            read2 = concatenated.files$read2,
                            adaptor1 = info$adaptor1,
                            adaptor2 = info$adaptor2,
                            barcode1 = info$library_index_seq,
                            barcode2 = info$library_index_seq_dual,
                            lineage = paste0(info$genus,"_",info$species,"_",info$specimen_id))
  
  write.csv(pipeline.in, paste0(sample.dir,"/",outfile), row.names = F)
}
# sample.dir: full path to the folder holding your 'fastq.gz' read files
# metadata.file: file name for the csv metadata sequencing file (assumes it's in the 'sample.dir')
# outfile: name the output file
# adaptor1: sequence info for the first adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)
# adaptor2: sequence info for the second adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)


# generate a shell script for DEDUPE CLEAN FILTER
dcf_shell <- function(sample.file.path, no.parallel=2, 
                      memory=40, cpu=40, sample.dir, source.dir,
                      reference, minid=0.25, project.name){
  
  # read in the sample file csv
  sample.file <- read.csv(sample.file.path, header=T)
  
  # create a loop to make the dcf call for each sample
  for (k in 1:nrow(sample.file)){
    # identify curreny sample
    curr.sample <- sample.file[k,]
    # write the script call
    py.call <- paste("cd", source.dir, ";",
                     "python Scripts/dedupe_clean_filter_reads.py",
                     "--dir", sample.dir,
                     "--file", paste0(source.dir,"/",sample.dir,"samples.csv"), # this is a bad line, only works if file is named samples.csv. think of something better!
                     "--sample", curr.sample$sample,
                     "--mem", memory,
                     "--CPU", cpu,
                     "--minid", minid,
                     "--ref", reference)
    write(py.call, file=paste0("~/Desktop/", project.name, "_dcf_parallel.txt"), append=T, sep="\n")
    #print(py.call)
  }
  # exit message
  write("echo 'You should be all done now.'", file=paste0("~/Desktop/", project.name, "_dcf_parallel.txt"), append=T, sep="\n")
  
  # call to screen telling you the file name and where to run it
  cat("Your shell script for running 'dedupe_clean_filter_reads.py' in parallel is written to:\n", 
      paste0("~/Desktop/", project.name, "_dcf_parallel.txt"),"\n")
  cat("Execute the command in parallel by copy/paste to your terminal:\n", 
      paste(paste0("parallel ", "-j ", no.parallel, " --bar ::::"), paste0(sample.dir, project.name, "_dcf_parallel.txt")))
}
# sample.file.path: full path to the samples.csv file on your local machine
# no.parallel: number of samples to process in parallel (probably <1/2 the 'cpu' command)
# memory: amount of RAM to make available
# cpu: number of available cores/cpu
# sample.dir: directory name that holds all your read files (end with "/")
# source.dir: directory name of the parent directory to sample.dir (full path!)
# reference: file name for the fasta of reference sequences (give path if not in source.dir)
# minid: bbmap 'minimum identity' value for mapping raw reads to the reference sequences (be lenient)


# generate a shell script for TRINITY ASSEMBLY 
trinity_shell <- function(sample.file.path, no.parallel=2,
                          memory=40, cpu=40, sample.dir, 
                          source.dir, project.name){
  
  # read in the sample file csv
  sample.file <- read.csv(sample.file.path, header=T)
  
  # create a loop to make the assembly call for each sample
  for (k in 1:nrow(sample.file)){
    # identify current sample
    curr.sample <- sample.file[k,]
    # write the script call
    py.call <- paste("cd", source.dir, ";",
                     "python Scripts/trinity_filtered_assembly.py",
                     "--dir", sample.dir,
                     "--sample", curr.sample$sample,
                     "--mem", memory,
                     "--CPU", cpu)
    write(py.call, file=paste0("~/Desktop/", project.name, "_assembly_parallel.txt"), append=T, sep="\n")
  }
  # exit message
  write('echo "You should be all done now."', file=paste0("~/Desktop/", project.name, "_assembly_parallel.txt"), append=T, sep="\n")
  
  # call to screen telling you the file name and where to run it
  cat("Your shell script for running 'dedupe_clean_filter_reads.py' in parallel is written to:\n", 
      paste0("~/Desktop/", project.name, "_assembly_parallel.txt"),"\n")
  cat("Execute the command in parallel by copy/paste to your terminal:\n", 
      paste(paste0("parallel ", "-j ", no.parallel, " --bar ::::"), paste0(sample.dir, project.name, "_assembly_parallel.txt")))
}
# sample.file.path: full path to the samples.csv file on your local machine
# no.parallel: number of samples to process in parallel (probably <1/2 the 'cpu' command)
# memory: amount of RAM to make available
# cpu: number of available cores/cpu
# sample.dir: directory name that holds all your read files (end with "/")
# source.dir: directory name of the parent directory to sample.dir (full path!)
# project.name: simple project title to track output files

# python Scripts/trinity_filtered_assembly.py --sample Pletholax_gracilis_WAM_R154023_350790 --dir OpTest_PhyFilter/ --mem 40 --CPU 40


# generate a shell script for MATCH CONTIGS to TARGETS
match2targets_shell <- function(sample.file.path, no.parallel=10,
                                sample.dir, source.dir, project.name, 
                                evalue=1e-30, targets.file.path){
  
  # read in the sample file csv
  sample.file <- read.csv(sample.file.path, header=T)
  
  # start the file by making the output directory
  dir.call <- paste0("mkdir ",source.dir,"/",sample.dir,"matches")
  write(dir.call, file=paste0("~/Desktop/", project.name, "_match2targets_parallel.txt"), append=T, sep="\n")
  
  # create a loop to make the assembly call for each sample
  for (k in 1:nrow(sample.file)){
    # identify current sample
    curr.sample <- sample.file[k,]
    # write the script call
    py.call <- paste("cd", source.dir, ";",
                     "python Scripts/match_contigs_to_probes.py",
                     "--dir", sample.dir,
                     "--sample", curr.sample$sample,
                     "--evalue", evalue,
                     "--db", targets.file.path)
    write(py.call, file=paste0("~/Desktop/", project.name, "_match2targets_parallel.txt"), append=T, sep="\n")
  }
  # exit message
  #write('echo "You should be all done now."', file=paste0("~/Desktop/", project.name, "_match2targets_parallel.txt"), append=T, sep="\n")
  
  # call to screen telling you the file name and where to run it
  cat("Your shell script for running 'match_contigs_to_probes.py' in parallel is written to:\n", 
      paste0("~/Desktop/", project.name, "_match2targets_parallel.txt"),"\n")
  cat("Execute the command in parallel by copy/paste to your terminal:\n", 
      paste(paste0("parallel ", "-j ", no.parallel, " --bar ::::"), paste0(sample.dir, project.name, "_match2targets_parallel.txt")))
}
# sample.file.path: full path to the samples.csv file on your local machine
# no.parallel: number of samples to process in parallel (probably <1/2 the 'cpu' command)
# sample.dir: directory name that holds all your read files (end with "/")
# source.dir: directory name of the parent directory to sample.dir (full path!)
# project.name: simple project title to track output files
# evalue: accuracy score to distinguish genuine matches from ones ocurring due to chance
# targets.file.path: file name for the fasta of target sequences (give path if not in source.dir)

# python ~/squamateUCE/match_contigs_to_probes.py --blat ~/bin/blat --sample Anolis_carolinensis --dir /scratch/drabosky_flux/sosi/uce_test/ --evalue 1e-30 --db squamate_AHE_UCE_genes_loci.fasta


# generate a shell script to GENERATE PSEUDO-REFERENCE GENOMES 
prg_shell <- function(sample.file.path, no.parallel=2,
                          sample.dir, source.dir, project.name, 
                          keep = "easy_recip_match"){
  
  # read in the sample file csv
  sample.file <- read.csv(sample.file.path, header=T)
  
  # start the file by making the output directory
  dir.call <- paste0("mkdir ",source.dir,"/",sample.dir,"PRG")
  write(dir.call, file=paste0("~/Desktop/", project.name, "_PRG_parallel.txt"), append=T, sep="\n")
  
  # create a loop to make the assembly call for each sample
  for (k in 1:nrow(sample.file)){
    # identify current sample
    curr.sample <- sample.file[k,]
    # write the script call
    py.call <- paste("cd", source.dir, ";",
                     "python Scripts/make_PRG.py",
                     "--dir", sample.dir,
                     "--lineage", curr.sample$lineage,
                     "--file", paste0(source.dir,"/",sample.dir,"samples.csv"), # this is a bad line, only works if file is named samples.csv. think of something better!
                     "--keep", keep)
    write(py.call, file=paste0("~/Desktop/", project.name, "_PRG_parallel.txt"), append=T, sep="\n")
  }
  # exit message
  #write('echo "You should be all done now."', file=paste0("~/Desktop/", project.name, "_PRG_parallel.txt"), append=T, sep="\n")
  
  # remove any duplicated lines
  dd.call <- paste("awk '!x[$0]++'", paste0("~/Desktop/", project.name, "_PRG_parallel.txt"), "| tee", 
                     paste0("~/Desktop/", project.name, "_PRG_parallel.txt"), "> /dev/null")
  system(dd.call)
  #cat(dd.call)

  # call to screen telling you the file name and where to run it
  cat("Your shell script for running 'make_PRG.py' in parallel is written to:\n", 
      paste0("~/Desktop/", project.name, "_PRG_parallel.txt"),"\n")
  cat("Execute the command in parallel by copy/paste to your terminal:\n", 
      paste(paste0("parallel ", "-j ", no.parallel, " --bar ::::"), paste0(sample.dir, project.name, "_PRG_parallel.txt")))
}
# sample.file.path: full path to the samples.csv file on your local machine
# no.parallel: number of samples to process in parallel (probably <1/2 the 'cpu' command)
# sample.dir: directory name that holds all your read files (end with "/")
# source.dir: directory name of the parent directory to sample.dir (full path!)
# project.name: simple project title to track output files
# keep: which matches to keep? see https://github.com/singhal/SqCL section '5. Match assemblies to original targets' for explanation
       # my suggestion is to leave on 'easy_recip_match' unless you're doing some deep digging.

# python ~/squamateUCE/make_PRG.py --lineage l1 --file /scratch/drabosky_flux/sosi/uce_test/samples.csv --dir /scratch/drabosky_flux/sosi/uce_test/ --keep easy_recip_match


# generate a shell script for QUALITY CONTROL
qc2_shell <- function(sample.file.path, no.parallel=1,
                      sample.dir, source.dir, project.name, 
                      output.dir, sample.or.lineage=c("lineage", "sample")){
  
  # read in the sample file csv
  sample.file <- read.csv(sample.file.path, header=T)
  
  # start the file by making the output directory
  dir.call <- paste0("mkdir ",source.dir,"/",sample.dir,output.dir)
  write(dir.call, file=paste0("~/Desktop/", project.name, "_qc2_parallel.txt"), append=T, sep="\n")
  
  column <- which(colnames(sample.file)==sample.or.lineage)
  
  # create a loop to make the assembly call for each sample
  for (k in 1:nrow(sample.file)){
    # identify current sample
    curr.sample <- sample.file[k,]
    # write the script call
    py.call <- paste("cd", source.dir, ";",
                     "python Scripts/quality_2_assembly.py",
                     "--dir", sample.dir,
                     "--ind", curr.sample[column],
                     "--file", paste0(source.dir,"/",sample.dir,"samples.csv"), # this is a bad line, only works if file is named samples.csv. think of something better!
                     "--outdir", paste0(source.dir,"/",sample.dir,output.dir))
    write(py.call, file=paste0("~/Desktop/", project.name, "_qc2_parallel.txt"), append=T, sep="\n")
  }
  # exit message
  #write('echo "You should be all done now."', file=paste0("~/Desktop/", project.name, "_qc2_parallel.txt"), append=T, sep="\n")
  
  # remove any duplicated lines
  dd.call <- paste("awk '!x[$0]++'", paste0("~/Desktop/", project.name, "_qc2_parallel.txt"), "| tee", 
                   paste0("~/Desktop/", project.name, "_qc2_parallel.txt"), "> /dev/null")
  system(dd.call)
  
  # call to screen telling you the file name and where to run it
  cat("Your shell script for running 'quality_2_assembly.py' in parallel is written to:\n", 
      paste0("~/Desktop/", project.name, "_qc2_parallel.txt"),"\n")
  cat("Execute the command by copy/paste to your terminal:\n", "(don't run in parallel, keep -j 1)\n", 
      paste(paste0("parallel ", "-j ", no.parallel, " --bar ::::"), paste0(sample.dir, project.name, "_qc2_parallel.txt")))
}
# sample.file.path: full path to the samples.csv file on your local machine
# no.parallel: number of samples to process in parallel (probably <1/2 the 'cpu' command)
# sample.dir: directory name that holds all your read files (end with "/")
# source.dir: directory name of the parent directory to sample.dir (full path!)
# project.name: simple project title to track output files
# output.dir: name the output directory for these results
# sample.or.lineage: have you assembled each individual separately (yes: "lineage")
      # or have you assembled all individuals of a species to a single sample (yes: "sample")

# python Scripts/quality_2_assembly.py --file Optimization/samples.csv --ind Pletholax_gracilis_WAM_R154023_350790 --dir Optimization/ --outdir Optimization/quality_scores


