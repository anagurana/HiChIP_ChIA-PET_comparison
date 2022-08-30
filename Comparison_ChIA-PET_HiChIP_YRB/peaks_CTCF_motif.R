suppressMessages(library(Rsubread))
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(ggridges)

# Function that adds column about the CTCF motif fo the peaks data

add_motif_info_to_peaks=function(peaks, CTCF_data){
  peaks.ctcfoverlap = findOverlaps(peaks, CTCF_data)
  
  numbers_CTCF=c()
  max_motif_CTCF=c()
  average_motif_CTCF=c()
  v=1
  
  npeaks=NROW(peaks)
  noverlap=NROW(peaks.ctcfoverlap)
  
  for (p in 1:npeaks){
    id_motifs_CTCF_per_peak=c()
    while (v <= noverlap && queryHits(peaks.ctcfoverlap)[v] == p){
      id_motifs_CTCF_per_peak=c(id_motifs_CTCF_per_peak, subjectHits(peaks.ctcfoverlap[v]))
      v=v+1
    }
    
    numbers_CTCF=c(numbers_CTCF, NROW(id_motifs_CTCF_per_peak))
    
    if (NROW(id_motifs_CTCF_per_peak) == 0){
      max_motif_CTCF=c(max_motif_CTCF, 0)
      average_motif_CTCF=c(average_motif_CTCF, 0)
    }
    else{
      motifs_CTCF_per_peak=CTCF_data[id_motifs_CTCF_per_peak]
      max_motif_CTCF=c(max_motif_CTCF, max(motifs_CTCF_per_peak$score))
      average_motif_CTCF=c(average_motif_CTCF, sum(motifs_CTCF_per_peak$score)/NROW(motifs_CTCF_per_peak))
    }
  }
  
  peaks$numbers_CTCF=numbers_CTCF
  peaks$max_motif_CTCF=max_motif_CTCF
  peaks$average_motif_CTCF=average_motif_CTCF
  return(peaks)
}

get_peaks_with_motif=function(peaks_with_mowif){
  k=1
  numbers=c()
  for (i in peaks_with_mowif$numbers_CTCF){
    if (i != 0){
      numbers=c(numbers, k)
      k=k+1
    }
    else{
      k=k+1
    }
  }
  peaks_non_empty=peaks_with_mowif[numbers]
  # max_CTCF=sum(peaks_non_empty$max_motif_CTCF)/NROW(peaks_non_empty)
  # average_CTCF=sum(peaks_non_empty$average_motif_CTCF)/NROW(peaks_non_empty)
  # return(c(max_CTCF, average_CTCF))
  return(peaks_non_empty)
}

# source("/Users/aleksander.ostruk/Documents/asia_master/peaks_analysis_VennDiag.R")
# 
# experiment="hichip"
# # cell="GM19238"

# # Read CTCF motif peaks
# CTCF.hg38 = read.table(file.path(path, "CTCF_hg38.txt"))
# CTCF.hg38 = GRanges(seqnames = CTCF.hg38$V1, IRanges(start = CTCF.hg38$V3, width = CTCF.hg38$V4), strand = ifelse(CTCF.hg38$V5=='p', "+", "-"), score=CTCF.hg38$V6)
# CTCF.hg38 = CTCF.hg38[CTCF.hg38$score >0]

# Reduce the number of peaks
# peaks_gm19238_chiapet=filter_peaks(read_peak_file('chiapet', cell))
# peaks_gm19238_hichip=filter_peaks(read_peak_file('hichip', cell))
# peaks_gm19238_chipseq=filter_peaks(read_peak_file('chipseq', cell))

# peaks_gm19238_chiapet=filter_peaks(read_peak_file(experiment, 'GM19238'))
# peaks_gm19238_hichip=filter_peaks(read_peak_file(experiment, 'GM19239'))
# peaks_gm19238_chipseq=filter_peaks(read_peak_file(experiment, 'GM19240'))

# # Find how many peaks in each technology contain CTCF motif
# npeaks_chiapet=NROW(peaks_gm19238_chiapet)
# npeaks_with_ctcf_chiapet=NROW(subsetByOverlaps(peaks_gm19238_chiapet, CTCF.hg38))
# frac_with_motif_chiapet=npeaks_with_ctcf_chiapet/npeaks_chiapet*100

# npeaks_hichip=NROW(peaks_gm19238_hichip)
# npeaks_with_ctcf_hichip=NROW(subsetByOverlaps(peaks_gm19238_hichip, CTCF.hg38))
# frac_with_motif_hichip=npeaks_with_ctcf_hichip/npeaks_hichip*100
# 
# npeaks_chipseq=NROW(peaks_gm19238_chipseq)
# npeaks_with_ctcf_chipseq=NROW(subsetByOverlaps(peaks_gm19238_chipseq, CTCF.hg38))
# frac_with_motif_chipseq=npeaks_with_ctcf_chipseq/npeaks_chipseq*100

# ctcf_motif_counts=data.frame(Technology=c("ChIA-PET", "HiChIP", "ChIP-Seq"), 
#                              Total_peaks=c(npeaks_with_ctcf_chiapet, npeaks_with_ctcf_hichip, npeaks_with_ctcf_chipseq),
#                              Peaks_with_CTCF=c(npeaks_with_ctcf_chiapet, npeaks_with_ctcf_hichip, npeaks_with_ctcf_chipseq),
#                              "Per_peaks_with_CTCF"=c(frac_with_motif_chiapet, frac_with_motif_hichip, frac_with_motif_chipseq))


# ## Add motif info to peaks 
# # for all peaks in each technology
# peaks_chiapet_with_motif_info=add_motif_info_to_peaks(peaks_gm19238_chiapet, CTCF.hg38)
# peaks_hichip_with_motif_info=add_motif_info_to_peaks(peaks_gm19238_hichip, CTCF.hg38)
# peaks_chipseq_with_motif_info=add_motif_info_to_peaks(peaks_gm19238_chipseq, CTCF.hg38)
# 
# # for common peaks for 3 technologies
# peaks_between_3_repl_with_motif_info=add_motif_info_to_peaks(common_overlaps$com123, CTCF.hg38)
# 
# # for common peaks between 2 technologies
# overlap14=subsetByOverlaps(common_overlaps$com13, common_overlaps$com123, invert = TRUE)
# peaks_14_with_motif_info=add_motif_info_to_peaks(overlap14, CTCF.hg38)
# overlap24=subsetByOverlaps(common_overlaps$com12, common_overlaps$com123, invert = TRUE)
# peaks_24_with_motif_info=add_motif_info_to_peaks(overlap24, CTCF.hg38)
# overlap34=subsetByOverlaps(common_overlaps$com23, common_overlaps$com123, invert = TRUE)
# peaks_34_with_motif_info=add_motif_info_to_peaks(overlap34, CTCF.hg38)
# 
# # for unique peaks for each technology
# overlap1=subsetByOverlaps(peaks_gm19238_chiapet, common_overlaps$com13, invert = TRUE)
# overlap_raw1=subsetByOverlaps(overlap1, overlap24, invert = TRUE)
# peaks_unique_chiapet_with_motif_info=add_motif_info_to_peaks(overlap_raw1, CTCF.hg38)
# 
# overlap2=subsetByOverlaps(peaks_gm19238_hichip, common_overlaps$com12, invert = TRUE)
# overlap_raw2=subsetByOverlaps(overlap2, overlap34, invert = TRUE)
# peaks_unique_hichip_with_motif_info=add_motif_info_to_peaks(overlap_raw2, CTCF.hg38)
# 
# overlap3=subsetByOverlaps(peaks_gm19238_chipseq, common_overlaps$com23, invert = TRUE)
# overlap_raw3=subsetByOverlaps(overlap3, overlap14, invert = TRUE)
# peaks_unique_chipseq_with_motif_info=add_motif_info_to_peaks(overlap_raw3, CTCF.hg38)
# 
# # Make a GRangesList with all regions included
# all_peaks_with_motifs_info=GRangesList(
#   peaks_chiapet_with_motif_info, peaks_hichip_with_motif_info, peaks_chipseq_with_motif_info,
#   peaks_unique_chiapet_with_motif_info, peaks_unique_hichip_with_motif_info, peaks_unique_chipseq_with_motif_info, 
#   peaks_14_with_motif_info, peaks_24_with_motif_info, peaks_34_with_motif_info, 
#   peaks_between_3_repl_with_motif_info)

# i=1
# peaks_with_motifs_allregions=c()
# for (i in 1:NROW(all_peaks_with_motifs_info)){
#   peaks_with_motifs_allregions=c(peaks_with_motifs_allregions, get_peaks_with_motif(all_peaks_with_motifs_info[[i]]))}

## Create box plots, 1 - total number of peaks for 3 tech, 2 - common regions
# make a list with CTCF strongest motid value
a1=peaks_with_motifs_allregions[[1]]$max_motif_CTCF
a2=peaks_with_motifs_allregions[[2]]$max_motif_CTCF
a3=peaks_with_motifs_allregions[[3]]$max_motif_CTCF
a4=peaks_with_motifs_allregions[[4]]$max_motif_CTCF
a5=peaks_with_motifs_allregions[[5]]$max_motif_CTCF
a6=peaks_with_motifs_allregions[[6]]$max_motif_CTCF
a7=peaks_with_motifs_allregions[[7]]$max_motif_CTCF
a8=peaks_with_motifs_allregions[[8]]$max_motif_CTCF
a9=peaks_with_motifs_allregions[[9]]$max_motif_CTCF
a10=peaks_with_motifs_allregions[[10]]$max_motif_CTCF

all_peaks = list(a1, a2, a3)
common_peaks <- list(a4,a8,a5,a7, a10, a9, a6)

# change the names of the elements of the list :
names(all_peaks) <- c(paste("All ChIA-PET peaks"),
              paste("All HiChIP peaks"),
              paste("All ChIP-Seq peaks"))
names(common_peaks) <- c(paste("\nUnique \npeaks\n ChIA-PET"),
              paste("\nCommon \npeaks\nChIA-PET & HiChIP"),
              paste("\nUnique \npeaks\n HiChIP"),
              paste("\nCommon \npeaks\n ChIP-Seq & ChIA-PET"),
              paste("\nCommon \npeaks\n for 3 technologies"),
              paste("\nCommon \npeaks\n ChIP-Seq & HiChIP"),
              paste("\nUnique \npeaks\n ChIP-Seq"))

# final Boxplots
values=c("blue", "green","red", "#9496ff","#73cb91", "#b1fe88","#db4e97", "#d36f4f","#ed8f48", "#ff8f92")
boxplot(all_peaks , aes(x=names(all_peaks)), col=values[1:3] , ylab="Strongest CTCF motif score",
        main="Comparison of distribution of Strongest CTCF Motif Score for 3 technologies)")
par(mgp = c(6, 3.5,0), mar=c(5, 8, 5, 5))
boxplot(common_peaks , aes(x=names(C)), col=values[4:10] , ylab="Strongest CTCF motif score", 
        main="Distribution of Strongest CTCF motif score for GM19238", cex.lab=1.5, cex.axis=1.1)
dev.off()

## Create histograms, 1 - total number of peaks for 3 tech, 2 - common regions
# define the grid and margins
par(
  mfrow=c(3,3),
  mar=c(6,6,6,6)
)
line<-par(lwd=5)
p1 <- hist(peaks_with_motifs_allregions[[1]]$max_motif_CTCF)                  
p2 <- hist(peaks_with_motifs_allregions[[2]]$max_motif_CTCF)   
p3 <- hist(peaks_with_motifs_allregions[[3]]$max_motif_CTCF)  
p4 <- hist(peaks_with_motifs_allregions[[4]]$max_motif_CTCF) 
p5 <- hist(peaks_with_motifs_allregions[[5]]$max_motif_CTCF) 
p6 <- hist(peaks_with_motifs_allregions[[6]]$max_motif_CTCF) 
p7 <- hist(peaks_with_motifs_allregions[[7]]$max_motif_CTCF) 
p8 <- hist(peaks_with_motifs_allregions[[8]]$max_motif_CTCF) 
p9 <- hist(peaks_with_motifs_allregions[[9]]$max_motif_CTCF) 
p10 <- hist(peaks_with_motifs_allregions[[10]]$max_motif_CTCF) 

plot( p1, col="blue", main="All ChIA-PET peaks", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=3)  
plot( p2, col="green", main="All HiChIP peaks", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=3)
plot( p3, col="red", main="All ChIP-Seq peaks", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=3) 

plot( p4, col="#9496ff", main="Unique peaks in ChIA-PET", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4) 
plot( p8, col="#73cb91", main="Common peaks ChIA-PET & HiChIP", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4)
plot( p5, col="#b1fe88", main="Unique peaks in HiChIP", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4) 
plot( p7, col="#db4e97", main="Common peaks ChIP-Seq & ChIA-PET", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4)
plot( p10, col="#d36f4f", main="Common peaks for 3 technologies", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4)
plot( p9, col="#ed8f48", main="Common peaks ChIP-Seq & HiChIP", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4)
plot.new()
plot( p6, col="#ff8f92", main="Unique peaks in ChIP-Seq", xlab="Strongest CTCF motif score", cex.lab=2.5, cex.axis=2.5, cex.main=4)

dev.off()
































