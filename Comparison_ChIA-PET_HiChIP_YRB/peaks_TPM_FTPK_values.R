# BiocManager::install("Rsamtools")
# install.packages("data.table", type="source", dependencies=TRUE)
library(data.table)
suppressMessages(library(Rsamtools))

source("/Users/aleksander.ostruk/Documents/asia_master/functions_ChIA-PET_HiChIP_YRB.R")
  
#Get the number of reads for each peak from BAM file
calculate_reads_count=function(peaks, bam_file){
  reads_count=c()
  for (i in 1:NROW(peaks)){
    which <- GRanges(seqnames(peaks)[i], IRanges(start(peaks)[i], 
                                                 end(peaks)[i]))
    what = c("rname", "pos")
    param = ScanBamParam(which=which, what=what)
    # bam = scanBam(bamFile, param=param)
    summary_count=countBam(bam_file, param=param)
    read_count=summary_count[["records"]]
    reads_count=c(reads_count, read_count)
  }
  return(reads_count)
}

#Calculate the FPKM value
calculate_FPKM_TPM=function(peaks, reads_count){
  RC=reads_count
  SF=1000000
  PL=width(reduced_peaks)
  RPM=RC*SF/sum(RC)
  FPKM=RPM/PL*1000
  
  #Calculate the TPM value
  RPK=RC/PL*1000
  TPM=RPK*SF/sum(RPK)
  
  return(list(FPKM, TPM))
}

find_stat_values=function(overlapIDs, peaks, stat_values){
  logicValue=c()
  i=1
  for (p in 1:length(peaks)){
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
  
  tableValid=data.table(FPKM=stat_values[[1]], TPM=stat_values[[2]], isCommon=logicValue)
  
  true.table=tableValid[isCommon==TRUE]
  false.table=tableValid[isCommon==FALSE]
  
  trueFPKM=mean(true.table$FPKM)
  falseFPKM=mean(false.table$FPKM)
  
  trueTPM=mean(true.table$TPM)
  falseTPM=mean(false.table$TPM)
  stat_numbers=c(trueFPKM, falseFPKM, trueTPM ,falseTPM)
  return(stat_numbers)
}

#Data for the triple Venn Diagrams
experiment_1="chiapet"
experiment_2="hichip"
experiment_3="chipseq"
cell_1="GM19238"
cell_2="GM19238"
cell_3="GM19238"

main_exp="chipseq"

#Read and Reduce peaks
peak_file=read_peak_file(experiment_1, cell_1)
peaks_1=filter_peaks(peak_file)

peak_file=read_peak_file(experiment_2, cell_2)
peaks_2=filter_peaks(peak_file)

peak_file=read_peak_file(experiment_3, cell_3)
peaks_3=filter_peaks(peak_file)

#Read the bam file and calculate reads count for every peak
bam_file="/Users/aleksander.ostruk/Documents/asia_master/chipseq/GM19238/GM19238_CTCFall_1_subset_check_samtools_dedup.bam"
reads_count=calculate_reads_count(reduced_peaks, bam_file)
stat_values = calculate_FPKM_TPM(reduced_peaks, reads_count)

all_common_overpals=find_all_relative_overlaps(peaks_1, peaks_2, peaks_3, main_exp)

if (main_exp=="chiapet"){
  reduced_peaks=peaks_1
  
  # 10 region
  numbers_vs_all=find_stat_values(all_common_overpals$vs_all, reduced_peaks, stat_values)
  # 8 reion
  numbers_vs_ChIPSeq=find_stat_values(all_common_overpals$vs_ChIPSeq, reduced_peaks, stat_values)
  # 7 region
  numbers_vs_HiChIP=find_stat_values(all_common_overpals$vs_HiChIP, reduced_peaks, stat_values)
  # 4 region
  numbers_raw_1=find_stat_values(all_common_overpals$raw, reduced_peaks, stat_values)
  return(numbers_vs_all, numbers_vs_ChIPSeq, numbers_vs_HiChIP, numbers_raw_1)
} else if (main_exp=="hichip"){
  reduced_peaks=peaks_2
  # 10 region
  numbers_vs_all=find_stat_values(all_common_overpals$vs_all, reduced_peaks, stat_values)
  # 8 reion
  numbers_vs_ChIPSeq=find_stat_values(all_common_overpals$vs_ChIPSeq, reduced_peaks, stat_values)
  # # 9 region
  numbers_vs_ChIAPET=find_stat_values(all_common_overpals$vs_ChIAPET, reduced_peaks, stat_values)
  # # 5 region
  numbers_raw_2=find_stat_values(all_common_overpals$raw, reduced_peaks,stat_values)
  return(numbers_vs_all, numbers_vs_ChIPSeq, numbers_vs_ChIAPET, numbers_raw_2)
} else if (main_exp=="chipseq"){
  reduced_peaks=peaks_3
  # 10 region
  numbers_vs_all=find_stat_values(all_common_overpals$vs_all, reduced_peaks, stat_values)
  # # 9 region
  numbers_vs_ChIAPET=find_stat_values(all_common_overpals$vs_ChIAPET, reduced_peaks, stat_values)
  # 8 reion
  numbers_vs_HiChIP=find_stat_values(all_common_overpals$vs_HiChIP, reduced_peaks, stat_values)
  # # 6 region
  numbers_raw_3=find_stat_values(all_common_overpals$raw, reduced_peaks, stat_values)
  return(numbers_vs_all, numbers_vs_ChIAPET, numbers_vs_HiChIP, numbers_raw_3)
}














