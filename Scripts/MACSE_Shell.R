setwd("~/Desktop/Individual_FASTA")
#alignment.files <- align.files
alignment.files <- dir(getwd(), pattern =".fasta")
#iq.path <- "/Applications/iqtree-1.7.11/bin/iqtree"
macse.path <- "java -jar /Applications/macse_v2.03.jar"

## Create a function to split up the alignment files into a series of shell scripts to run through MACSE
macse.shell <- function(in.files, batch.size=10, prog=c("refine", "align", "export", "refineLemmon")){
  curr.dir <- getwd()
  dir.create(paste0(curr.dir,"/MACSE_Shells")) # make a directory for the new files
  
  splits <- split(in.files, ceiling(seq_along(in.files)/batch.size))
  for (k in 1:length(splits)){
    cd.call <- paste("cd", curr.dir) 
    write.table(cd.call, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
    for (g in 1:length(splits[[k]])){
      if (prog == "refine") {curr.call <- paste(macse.path, "-prog", "refineAlignment", "-align", paste0(curr.dir,"/",splits[[k]][g]), "-stop 10")}
      else if (prog == "refineLemmon") {curr.call <- paste(macse.path, "-prog", "refineAlignment", "-align", paste0(curr.dir,"/",splits[[k]][g]),
                                                           "-optim 1 -local_realign_init 0.1 -local_realign_dec 0.1 -fs 10")}
      else if (prog == "align") {curr.call <- paste(macse.path, "-prog", "alignSequences", "-seq_lr", paste0(curr.dir,"/",splits[[k]][g]), "-stop_lr 10")}
      else if (prog == "export") {curr.call <- paste(macse.path, "-prog", "exportAlignment", "-align", paste0(curr.dir,"/",splits[[k]][g]), "-stop 10")}
      write.table(curr.call, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
      permission.call <- paste("chmod =rwx,g+s", paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"))
      system(permission.call)
    }
    exit.msg <- paste("echo finished batch", k, "of", length(splits))
    write.table(exit.msg, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
    fix.exclam <- paste("perl -pi -w -e 's/!/N/g;' *.fasta")
    write.table(fix.exclam, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
    exit.call <- paste("exit 0")
    write.table(exit.call, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
  }
  shell.files <- dir(paste0(curr.dir,"/MACSE_Shells"), pattern=".sh")
  cd.call <- paste("cd", curr.dir) 
  shell.calls <- NULL
  for (q in 1:length(shell.files)){shell.calls <- append(shell.calls, paste("source", paste0(curr.dir, "/MACSE_Shells/", shell.files[q]), "&"))}
  combo <- c(cd.call, shell.calls, "exit 0")
  write.table(combo, file=paste0(curr.dir,"/MACSE_Shells/", "MASTER_MACSE_Shell.sh"), append=F, row.names=F, col.names=F, quote=F)
  #exit.msg <- paste("echo all done, you are a rockstar!")
  #write.table(exit.msg, file=paste0(curr.dir,"/MACSE_Shells/", "MACSE_IQTREE_Shell.sh"), append=T, row.names=F, col.names=F, quote=F)
  #exit.call <- paste("exit 0")
  #write.table(exit.call, file=paste0(curr.dir,"/MACSE_Shells/", "MASTER_MACSE_Shell.sh"), append=T, row.names=F, col.names=F, quote=F)
  permission.call <- paste("chmod =rwx,g+s", paste0(curr.dir,"/MACSE_Shells/", "MASTER_MACSE_Shell.sh"))
  system(permission.call)
}
macse.shell(alignment.files, batch.size=20, prog=c("refineLemmon"))

### IMPORTANT! MAKE SURE TO DO THE STEP BELOW! ###

# when finished remove the '!' character and replace it with 'N'
getwd() # check to make sure you're in the right folder
rename <- paste0("perl -pi -w -e 's/", "!", "/", "N", "/g;' *.fasta")
system(rename)




