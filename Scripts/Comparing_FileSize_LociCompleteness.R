
setwd("~/Desktop/contrees")
tree.files <- dir(pattern=".contree")

contrees <- read.tree("~/Desktop/contrees/UCE_loci.contrees")

all.names <- lapply(contrees, function(x) x$tip.label)
all.names <- unlist(all.names)
all.number <- table(all.names)
all.df <- data.frame(sample_name = names(all.number), loci = all.number)
write.csv(all.df, file="~/Desktop/contrees/sample_UCE_loci.csv", row.names = F)


metadata <- read.csv("/Users/ianbrennan/Downloads/Anilios_AusARG_Download/Project_Metadata_combined.csv")
trees <- read.csv("/Users/ianbrennan/Desktop/contrees/sample_ALL_loci.csv")

setdiff(metadata$sample_name, trees$sample_name)
setdiff(trees$sample_name, metadata$sample_name)

md2 <- left_join(metadata, trees)
md2[is.na(md2)] <- 0


ggplot(md2) +
  geom_point(aes(x=total_filesize_gb, y=loci_assembled_all), color="darkGrey") + 
  geom_point(aes(x=total_filesize_gb, y=loci_assembled_AHE), color="#4894c7") +
  geom_point(aes(x=total_filesize_gb, y=loci_assembled_UCE), color= "#ed622f") + 
  theme_bw() + geom_vline(xintercept=0.1, color="#d9d443", lwd=2, alpha=0.5) + 
  scale_x_continuous(limits=c(0, 6), breaks=c(seq(0,1,0.1),seq(1,6,1)))


