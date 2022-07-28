# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Rsamtools")
# library(Rsamtools)
# install.packages("data.table", type="source", dependencies=TRUE)

find_stat_values=function(overlapIDs, reducedPeaks){
  logicValue=c()
  i=1
  for (p in 1:length(reducedPeaks)){
    if (is.na(overlapIDs[i])){
      logicValue=c(logicValue, FALSE)
    }
    else if (p == overlapIDs[i]){
      logicValue=c(logicValue, TRUE)
      i=i+1
    }
    else{
      logicValue=c(logicValue, FALSE)
    }
  }
  
  tableValid=data.table(FPKM, TPM, isCommon=logicValue)
  
  true.table=tableValid[isCommon==TRUE]
  false.table=tableValid[isCommon==FALSE]
  
  trueFPKM=mean(true.table$FPKM)
  falseFPKM=mean(false.table$FPKM)
  
  trueTPM=mean(true.table$TPM)
  falseTPM=mean(false.table$TPM)
  stat_numbers=c(trueFPKM, falseFPKM, trueTPM ,falseTPM)
  return(stat_numbers)
}

source("/Users/aleksander.ostruk/Documents/asia_master/peaks_analysis_VennDiag.R")

read_bam_file <- function(exp_type, cell_line){
  dir_path=file.path(path, exp_type, cell_line)
  peaks_file=list.files(dir_path, pattern='.for.BROWSER.bam', 
                        recursive=FALSE, full.names = TRUE)[1]
  return(peaks_file)
}
  
path="/Users/aleksander.ostruk/Documents/asia_master"

#Data for the triple Venn Diagrams
experiment="chiapet"
cell="GM19238"

peakFile=read_peak_file(experiment, cell)
bamFile=read_bam_file(experiment, cell)

#Read and Reduce peaks
reducedPeaks=filter_peaks(peakFile)

#Get the number of reads for each peak 
readsCount=c()
for (i in 1:length(reducedPeaks)){
  which <- GRanges(seqnames(reducedPeaks)[i], IRanges(start(reducedPeaks)[i], end(reducedPeaks)[i]))
  what <- c("rname", "pos")
  param <- ScanBamParam(which=which, what=what)
  # bam <- scanBam(bamFile, param=param)
  summaryCount=countBam(bamFile, param=param)
  readCount=summaryCount[["records"]]
  readsCount=c(readsCount, readCount)
}

#Calculate the FPKM value
RC=readsCount
SF=1000000
PL=width(reducedPeaks)
RPM=RC*SF/sum(RC)
FPKM=RPM/PL*1000

#Calculate the TPM value
RPK=RC/PL*1000
TPM=RPK*SF/sum(RPK)

#Overlaps between clean experiment and common between 3 technologies
overlapIDs.123=unique(queryHits(findOverlaps(reducedPeaks, common_overlaps$com123)))

#Overlaps between clean experiment and common between 2 and 3 technology 
overlap24=subsetByOverlaps(common_overlaps$com23, common_overlaps$com123, invert = TRUE)
overlapIDs.24=unique(queryHits(findOverlaps(reducedPeaks, overlap24)))

#Overlaps between clean experiment and common between 1 and 2 technology 
overlap14=subsetByOverlaps(common_overlaps$com12, common_overlaps$com123, invert = TRUE)
overlapIDs.14=unique(queryHits(findOverlaps(reducedPeaks, overlap12)))

#Overlapping ids of unique peaks
overlap2=subsetByOverlaps(reducedPeaks, common_overlaps$com12, invert = TRUE)
overlap_raw2=subsetByOverlaps(overlap2, overlap24, invert = TRUE)
overlapIDs.raw=unique(queryHits(findOverlaps(reducedPeaks, overlap_raw2)))


numbers123=find_stat_values(overlapIDs.123, reducedPeaks)
numbers24=find_stat_values(overlapIDs.24, reducedPeaks)
numbers14=find_stat_values(overlapIDs.14, reducedPeaks)
numbers.raw=find_stat_values(overlapIDs.raw, reducedPeaks)











