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