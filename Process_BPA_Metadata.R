#library(dplyr)
process_bpa_metadata <- function(package.metadata, resource.metadata){
  resource.metadata$samp.no <- sapply(resource.metadata$Name, function(x) strsplit(x, "_")[[1]][1])
  package.metadata$samp.no  <- sapply(package.metadata$library_id, function(x) strsplit(x, "/")[[1]][2])
  
  resource.metadata$File.size..mb. <- resource.metadata$File.size..bytes./1000000
  
  agg.rm <- aggregate(resource.metadata$File.size..mb., list(resource.metadata$samp.no), sum)
  colnames(agg.rm) <- c("samp.no", "total.filesize.mb")
  
  pm.rm <- dplyr::inner_join(package.metadata, agg.rm)
  pm.rm$total_filesize_gb <- pm.rm$total.filesize.mb/1000
  
  return(pm.rm)
}

# e.g.

# pm <- read.csv("~/Downloads/bpa_7d47cf9a_20220729T0618/package_metadata/package_metadata_bpa_7d47cf9a_20220729T0618_ausarg-illumina-fastq.csv")
# rm <- read.csv("~/Downloads/bpa_7d47cf9a_20220729T0618/resource_metadata/resource_metadata_bpa_7d47cf9a_20220729T0618_ausarg-illumina-fastq.csv")
# process_bpa_metadata(package.metadata = pm, resource.metadata = rm)



# generate a sample file for downstream purposes
make_sample_info <- function(package.metadata, resource.metadata,
                             sample.dir, metadata.file, outfile="samples.csv",
                             out.path = "/home/ian/SqCL_Pipeline/Optimization",
                             adaptor1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG",
                             adaptor2 = "AATGATACGGCGACCACCGAGATCTACAC*ACACTCTTTCCCTACACGACGCTCTTCCGATCT"){
  
  # combine the two metadata files
  pm.rm <- process_bpa_metadata(package.metadata = package.metadata,
                                resource.metadata = resource.metadata)
  
  # read in the names of the sequence files
  seq.files <- resource.metadata$Name
  
  # make a dataframe after cleaning the information
  file.ids <- data.frame(filename = seq.files,
                         sample.id = sapply(seq.files, function(x) strsplit(x, "_")[[1]][1]),
                         barcode = sapply(seq.files, function(x) strsplit(x, "_")[[1]][5]), # was [3]
                         sample.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][6]), # was [4]
                         l.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][7]), # was [5]
                         direction = sapply(seq.files, function(x) strsplit(x, "_")[[1]][8])) # was [8]
  rownames(file.ids) <- NULL
  
  # read in the sample metadata file
  info <- pm
  info$samp.no <- sapply(info$library_id, function(x) strsplit(x, "/")[[1]][2])
  # make sure there aren't spaces in the specimen_id field
  info$specimen_id <- sub(" ", "_", info$specimen_id)
  # choose the appropriate data
  info <- dplyr::select(info, "library_id", "genus", "species", "specimen_id", 
                        #"library_index_seq_P7", "library_index_seq_P5")
                        "library_index_seq", "library_index_seq_dual", "samp.no")
  # FYI 'library_index_seq' = p7/i7 barcode; 'library_index_seq_dual' = p5/i5 barcode
  
  # make a column for the sample info
  info <- dplyr::mutate(info, sample = paste0(genus, "_", species, "_", specimen_id))
  # add the adaptor information to the file
  info$adaptor1 <- adaptor1
  info$adaptor2 <- adaptor2
  
  # identify all forward files
  forwards <- filter(file.ids, direction == "R1")
  # identify all reverse files
  reverses <- filter(file.ids, direction == "R2")
  
  r1.concat <- sapply(1:nrow(info), function(x) paste0(info$sample[[x]], "_", info$samp.no[[x]], "_R1_concat.fastq.gz"))
  
  
  concatenated.files <- NULL
  for (j in 1:nrow(info)){
    # choose current sample
    #curr.sample <- info$library_id[[j]]
    curr.sample <- info$samp.no[[j]]
    
    # filter the appropriate fwd/rev files
    curr.fwd <- filter(forwards, curr.sample == sample.id)
    curr.rev <- filter(reverses, curr.sample == sample.id)
    
    # make the outfile names
    out.fwd <- paste0(info[j,"sample"], "_", curr.sample, "_R1_concat.fastq.gz")
    out.rev <- paste0(info[j,"sample"], "_", curr.sample, "_R2_concat.fastq.gz")
    
    # identify the number of read file pairs
    n.reads <- nrow(curr.fwd)
    
    # make the bash call to concatenate the files (depends on the number of read files per sample!)
    if(n.reads == 1) {
      cat.fwd <- paste("mv", curr.fwd$filename[[1]], out.fwd)
      cat.rev <- paste("mv", curr.rev$filename[[1]], out.rev)
    }
    if(n.reads == 2) {
      cat.fwd <- paste("cat", curr.fwd$filename[[1]],
                       curr.fwd$filename[[2]],
                       ">>", out.fwd)
      cat.rev <- paste("cat", curr.rev$filename[[1]],
                       curr.rev$filename[[2]],
                       ">>", out.rev)    }
    if(n.reads == 3) {
      cat.fwd <- paste("cat", curr.fwd$filename[[1]],
                       curr.fwd$filename[[2]],
                       curr.fwd$filename[[3]],
                       ">>", out.fwd)
      cat.rev <- paste("cat", curr.rev$filename[[1]],
                       curr.rev$filename[[2]],
                       curr.rev$filename[[3]],
                       ">>", out.rev)    }
    if(n.reads == 4) {
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
    concatenated.files <- rbind(concatenated.files, data.frame(read1 = paste0(out.path,"/",out.fwd), 
                                                               read2 = paste0(out.path,"/",out.rev)))
  }
  
  # create the SqCL_Pipeline sample file
  pipeline.in <- data.frame(sample = paste0(info$sample, "_", info$samp.no), # previously info$library_id
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
# outpath: the path on the machine where your files will be stored
# adaptor1: sequence info for the first (P7/i7) adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)
# adaptor2: sequence info for the second (P5/i5) adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)
# total.readno: the total number of read files per sample (probably an even number between 2-8)

# takes as input the output of process_bpa_metadata




# pm <- read.csv("~/Downloads/bpa_b1d97a74_20230509T0956/package_metadata/package_metadata_bpa_b1d97a74_20230509T0956_ausarg-exon-capture.csv")
# rm <- read.csv("~/Downloads/bpa_b1d97a74_20230509T0956/resource_metadata/resource_metadata_bpa_b1d97a74_20230509T0956_ausarg-exon-capture.csv")

# generate a sample file for downstream purposes
BPA_sample_info <- function(package.metadata, resource.metadata,
                            adaptor1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG",
                            adaptor2 = "AATGATACGGCGACCACCGAGATCTACAC*ACACTCTTTCCCTACACGACGCTCTTCCGATCT"){
  
  # combine the two metadata files
  pm.rm <- process_bpa_metadata(package.metadata = package.metadata,
                                resource.metadata = resource.metadata)
  
  # read in the names of the sequence files
  seq.files <- resource.metadata$Name
  
  # make a dataframe after cleaning the information
  file.ids <- data.frame(filename = seq.files,
                         sample.id = sapply(seq.files, function(x) strsplit(x, "_")[[1]][1]),
                         barcode = sapply(seq.files, function(x) strsplit(x, "_")[[1]][5]), # was [3]
                         sample.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][6]), # was [4]
                         l.no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][7]), # was [5]
                         direction = sapply(seq.files, function(x) strsplit(x, "_")[[1]][8])) # was [8]
  rownames(file.ids) <- NULL
  
  # read in the sample metadata file
  info <- pm
  info$samp.no <- sapply(info$library_id, function(x) strsplit(x, "/")[[1]][2])
  # make sure there aren't spaces in the specimen_id field
  info$specimen_id <- sub(" ", "_", info$specimen_id)
  # choose the appropriate data
  info <- dplyr::select(info, "library_id", "genus", "species", "specimen_id", 
                        #"library_index_seq_P7", "library_index_seq_P5")
                        "library_index_seq", "library_index_seq_dual", "samp.no")
  # FYI 'library_index_seq' = p7/i7 barcode; 'library_index_seq_dual' = p5/i5 barcode
  
  # make a column for the sample info
  info <- dplyr::mutate(info, sample = paste0(genus, "_", species, "_", specimen_id))
  # add the adaptor information to the file
  info$adaptor1 <- adaptor1
  info$adaptor2 <- adaptor2
  
  # create the SqCL_Pipeline sample file
  pipeline.in <- data.frame(BPA_ID = info$samp.no,
                            sample = paste0(info$sample, "_", info$samp.no), # previously info$library_id
                            read1 = "",
                            read2 = "",
                            adaptor1 = info$adaptor1,
                            adaptor2 = info$adaptor2,
                            barcode1 = info$library_index_seq,
                            barcode2 = info$library_index_seq_dual,
                            lineage = paste0(info$genus,"_",info$species,"_",info$specimen_id))
  
  write.csv(pipeline.in, paste0(getwd(),"/BPA_SampleInfo.csv"), row.names = F)
}
# sample.dir: full path to the folder holding your 'fastq.gz' read files
# metadata.file: file name for the csv metadata sequencing file (assumes it's in the 'sample.dir')
# outfile: name the output file
# outpath: the path on the machine where your files will be stored
# adaptor1: sequence info for the first (P7/i7) adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)
# adaptor2: sequence info for the second (P5/i5) adaptor, with barcode replaced by '*' (don't touch this unless you really mean it)
# total.readno: the total number of read files per sample (probably an even number between 2-8)


BPA_sample_info <- function(package.metadata, resource.metadata,
                            adaptor1 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG",
                            adaptor2 = "AATGATACGGCGACCACCGAGATCTACAC*ACACTCTTTCCCTACACGACGCTCTTCCGATCT"){
  
  # read in the names of the sequence files
  seq.files <- resource.metadata$Name
  
  # make a dataframe after cleaning the information
  file.ids <- data.frame(filename = seq.files,
                         sample_id = sapply(seq.files, function(x) strsplit(x, "_")[[1]][1]),
                         barcode = sapply(seq.files, function(x) strsplit(x, "_")[[1]][5]), # was [3]
                         sample_no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][6]), # was [4]
                         lane_no = sapply(seq.files, function(x) strsplit(x, "_")[[1]][7]), # was [5]
                         direction = sapply(seq.files, function(x) strsplit(x, "_")[[1]][8])) # was [8]
  rownames(file.ids) <- NULL
  
  # read in the sample metadata file
  info <- package.metadata
  info$sample_id <- sapply(info$library_id, function(x) strsplit(x, "/")[[1]][2])
  # make sure there aren't spaces in the specimen_id field
  info$specimen_id <- sub(" ", "_", info$specimen_id)
  # choose the appropriate data
  info <- dplyr::select(info, "library_id", "genus", "species", "specimen_id", 
                        #"library_index_seq_P7", "library_index_seq_P5")
                        "library_index_seq", "library_index_seq_dual", "sample_id")
  # FYI 'library_index_seq' = p7/i7 barcode; 'library_index_seq_dual' = p5/i5 barcode
  
  # make a column for the sample info
  info <- dplyr::mutate(info, lineage = paste0(genus, "_", species, "_", specimen_id))
  # add the adaptor information to the file
  info$adaptor1 <- adaptor1
  info$adaptor2 <- adaptor2
  
  # split files into two dataframes we'll merge later
  file.ids.f <- dplyr::filter(file.ids, direction=="R1")
  file.ids.r <- dplyr::filter(file.ids, direction=="R2")
  
  # grab the appropriate barcodes
  file.ids.f$barcode <- sapply(file.ids.f$sample_id, function(x) info[which(info$sample_id==x),"library_index_seq"])
  file.ids.r$barcode <- sapply(file.ids.r$sample_id, function(x) info[which(info$sample_id==x),"library_index_seq_dual"])
  
  # grab the appropriate adaptors
  file.ids.f$adaptor <- adaptor1
  file.ids.r$adaptor <- adaptor2
  
  # grab the appropriate barcodes
  file.ids.f$lineage <- sapply(file.ids.f$sample_id, function(x) info[which(info$sample_id==x),"lineage"])
  file.ids.r$lineage <- sapply(file.ids.r$sample_id, function(x) info[which(info$sample_id==x),"lineage"])
  
  # combine the f/r file info together
  file.id.combo <- rbind(file.ids.f, file.ids.r)
  file.id.combo <- dplyr::select(file.id.combo, filename, sample_id, direction, lane_no, barcode, adaptor, lineage)
  
  # write the file.id.combo to file
  write.csv(file.id.combo, paste0(getwd(),"/BPA_ReadFileInfo.csv"), row.names = F)
  
  # determine the number of lanes and samples
  no.lanes <- unique(file.id.combo$lane_no)
  no.samps <- unique(file.id.combo$sample_id)
  # create an empty sample info dataframe
  samp.info <- NULL
  # run a loop across the file.id.combo file making a traditional sample.info.csv file
  for (k in 1:length(no.lanes)){
    curr.lane <- dplyr::filter(file.id.combo, lane_no == no.lanes[[k]])
    for (j in 1:length(no.samps)){
      curr.samp <- dplyr::filter(curr.lane, sample_id == no.samps[[j]])
      samp.info <- rbind(samp.info, data.frame(sample_id = no.samps[[j]],
                                               read1 = curr.samp[which(curr.samp$direction=="R1"),"filename"],
                                               read2 = curr.samp[which(curr.samp$direction=="R2"),"filename"],
                                               barcode1 = curr.samp[which(curr.samp$direction=="R1"),"barcode"],
                                               barcode2 = curr.samp[which(curr.samp$direction=="R2"),"barcode"],
                                               adaptor1 = adaptor1,
                                               adaptor2 = adaptor2,
                                               lineage = curr.samp$lineage[[1]]))
    }
  }
  # write the sample info dataframe to file
  write.csv(samp.info, paste0(getwd(),"/BPA_SampleInfo.csv"), row.names = F)
}








