
library(GenomicRanges)
# BiocManager::install("Rsubread")
library(Rsubread)

# FUNCTIONS
read_peaks=function(peaks_path){
  peaks.tb <- read.table(peaks_path)
  peaks.gr <- GRanges(seqnames = peaks.tb$V1, IRanges(
    start = peaks.tb$V2, end = peaks.tb$V3), strand = '*',
    pval=10**-peaks.tb$V8, qval = 10^-peaks.tb$V9)
  peaks.gr=sort(peaks.gr)
  return(peaks.gr)
}

convert_peak_to_df_featcnt=function(peak){
  peak_df=data.frame(
    Chr = seqnames(peak),
    Start = start(peak),
    End = end(peak),
    Strand = '*',
    GeneID = paste('Peak', 1:NROW(peak), sep = ''))
  return(peak_df)
}



# Read CTCF motif peaks
CTCF.hg38 = read.table(file.path("../data/CTCF_hg38.txt"))
CTCF.hg38 = GRanges(seqnames = CTCF.hg38$V1, 
                    IRanges(start = CTCF.hg38$V3, 
                            width = CTCF.hg38$V4), 
                    strand = ifelse(CTCF.hg38$V5=='p', "+", "-"), 
                    score=CTCF.hg38$V6)
CTCF.hg38 = CTCF.hg38[CTCF.hg38$score >0]

### HICHIP ###
# UPLOAD THE HICHIP PEAKS
peak_hichip_gm38=read_peaks("../out_chiapipe/hichip/GM19238_merged_replicates/GM19238_merged_replicates.for.BROWSER.spp.z6.broadPeak")
peak_hichip_gm38=subsetByOverlaps(peak_hichip_gm38, CTCF.hg38)
peak_hichip_gm38=IRanges::reduce(peak_hichip_gm38+1000)

peak_hichip_gm39=read_peaks("../out_chiapipe/hichip/gm19239/gm19239.for.BROWSER.spp.z6.broadPeak")
peak_hichip_gm39=subsetByOverlaps(peak_hichip_gm39, CTCF.hg38)
peak_hichip_gm39=IRanges::reduce(peak_hichip_gm39+1000)

peak_hichip_gm40=read_peaks("../out_chiapipe/hichip/gm19240/gm19240.for.BROWSER.spp.z6.broadPeak")
peak_hichip_gm40=subsetByOverlaps(peak_hichip_gm40, CTCF.hg38)
peak_hichip_gm40=IRanges::reduce(peak_hichip_gm40+1000)

peaks_reduced_with_motif_HiChIP=GRangesList("GM19238"=peak_hichip_gm38, "GM19239"=peak_hichip_gm39, "GM19240"=peak_hichip_gm40)

# ADD COVERAGE TRACK

reads_count_hichip_GM19238=featureCounts(
  "../out_chiapipe/hichip/GM19238_merged_replicates/GM19238_merged_replicates.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_HiChIP$GM19238), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_HiChIP$GM19238$coverage=unname(reads_count_hichip_GM19238$counts[,1])

reads_count_hichip_GM19239=featureCounts(
  "../out_chiapipe/hichip/gm19239/gm19239.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_HiChIP$GM19239), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_HiChIP$GM19239$coverage=unname(reads_count_hichip_GM19239$counts[,1])

reads_count_hichip_GM19240=featureCounts(
  "../out_chiapipe/hichip/gm19240/gm19240.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_HiChIP$GM19240), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_HiChIP$GM19240$coverage=unname(reads_count_hichip_GM19240$counts[,1])

trios=c("GM19238", "GM19239", "GM19240")
experiment="hichip"
for (cell_line in trios){
  write.table(unname(data.frame(seqnames(peaks_reduced_with_motif_HiChIP[[cell_line]]), 
                                start(peaks_reduced_with_motif_HiChIP[[cell_line]]), 
                                end(peaks_reduced_with_motif_HiChIP[[cell_line]]),
                                peaks_reduced_with_motif_HiChIP[[cell_line]]$coverage)),
              paste0("~/ctcf_motif_",experiment,"_",cell_line, ".txt"),
              sep="\t",row.names=FALSE, quote = FALSE)
}

# CHIAPET
# UPLOAD THE CHIAPET PEAKS
peak_chiapet_gm38=read_peaks("../out_chiapipe/chia_pet/gm19239/GM19238_merged_replicates.for.BROWSER.spp.z6.broadPeak")
peak_chiapet_gm38=subsetByOverlaps(peak_chiapet_gm38, CTCF.hg38)
peak_chiapet_gm38=IRanges::reduce(peak_chiapet_gm38+1000)

peak_chiapet_gm39=read_peaks("../out_chiapipe/chia_pet/gm19239/GM19239_merged_replicates.for.BROWSER.spp.z6.broadPeak")
peak_chiapet_gm39=subsetByOverlaps(peak_chiapet_gm39, CTCF.hg38)
peak_chiapet_gm39=IRanges::reduce(peak_chiapet_gm39+1000)

peak_chiapet_gm40=read_peaks("../out_chiapipe/chia_pet/gm19240/GM19240_merged_replicates.for.BROWSER.spp.z6.broadPeak")
peak_chiapet_gm40=subsetByOverlaps(peak_chiapet_gm40, CTCF.hg38)
peak_chiapet_gm40=IRanges::reduce(peak_chiapet_gm40+1000)

peaks_reduced_with_motif_ChIAPET=GRangesList("GM19238"=peak_chiapet_gm38, "GM19239"=peak_chiapet_gm39, "GM19240"=peak_chiapet_gm40)

# ADD COVERAGE TRACK

reads_count_chiapet_GM19238=featureCounts(
  "../out_chiapipe/chia_pet/gm19238/GM19238_merged_replicates.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIAPET$GM19238), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIAPET$GM19238$coverage=unname(reads_count_chiapet_GM19238$counts[,1])

reads_count_chiapet_GM19239=featureCounts(
  "../out_chiapipe/chia_pet/gm19239/GM19239_merged_replicates.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIAPET$GM19239), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIAPET$GM19239$coverage=unname(reads_count_chiapet_GM19239$counts[,1])

reads_count_chiapet_GM19240=featureCounts(
  "../out_chiapipe/chia_pet/gm19240/GM19240_merged_replicates.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIAPET$GM19240), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIAPET$GM19240$coverage=unname(reads_count_chiapet_GM19240$counts[,1])


experiment="chiapet"
for (cell_line in trios){
  write.table(unname(data.frame(seqnames(peaks_reduced_with_motif_ChIAPET[[cell_line]]), 
                                start(peaks_reduced_with_motif_ChIAPET[[cell_line]]), 
                                end(peaks_reduced_with_motif_ChIAPET[[cell_line]]),
                                peaks_reduced_with_motif_ChIAPET[[cell_line]]$coverage)),
              paste0("~/ctcf_motif_",experiment,"_", cell_line, ".txt"),
              sep="\t",row.names=FALSE, quote = FALSE)
}

### CHIPSEQ ###
# UPLOAD THE CHIPSEQ PEAKS
peak_chipseq_gm38=read_peaks("../out_chiapipe/hichip/GM19238_merged_replicates/GM19238_merged_replicates.for.BROWSER.spp.z6.broadPeak")
peak_chipseq_gm38=subsetByOverlaps(peak_chipseq_gm38, CTCF.hg38)
peak_chipseq_gm38=IRanges::reduce(peak_chipseq_gm38+1000)

peak_chipseq_gm39=read_peaks("../out_chiapipe/hichip/gm19239/gm19239.for.BROWSER.spp.z6.broadPeak")
peak_chipseq_gm39=subsetByOverlaps(peak_chipseq_gm39, CTCF.hg38)
peak_chipseq_gm39=IRanges::reduce(peak_chipseq_gm39+1000)

peak_chipseq_gm40=read_peaks("../out_chiapipe/hichip/gm19240/gm19240.for.BROWSER.spp.z6.broadPeak")
peak_chipseq_gm40=subsetByOverlaps(peak_chipseq_gm40, CTCF.hg38)
peak_chipseq_gm40=IRanges::reduce(peak_chipseq_gm40+1000)

peaks_reduced_with_motif_ChIPSeq=GRangesList("GM19238"=peak_chipseq_gm38, "GM19239"=peak_chipseq_gm39, "GM19240"=peak_chipseq_gm40)

# ADD COVERAGE TRACK

reads_count_chipseq_GM19238=featureCounts(
  "../out_chiapipe/hichip/GM19238_merged_replicates/GM19238_merged_replicates.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19238), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19238$coverage=unname(reads_count_chipseq_GM19238$counts[,1])

reads_count_chipseq_GM19239=featureCounts(
  "../out_chiapipe/hichip/gm19239/gm19239.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19239), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19239$coverage=unname(reads_count_chipseq_GM19239$counts[,1])

reads_count_chipseq_GM19240=featureCounts(
  "../out_chiapipe/hichip/gm19240/gm19240.for.BROWSER.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19240), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19240$coverage=unname(reads_count_chipseq_GM19240$counts[,1])

trios=c("GM19238", "GM19239", "GM19240")
experiment="chipseq"
for (cell_line in trios){
  write.table(unname(data.frame(seqnames(peaks_reduced_with_motif_ChIPSeq[[cell_line]]), 
                                start(peaks_reduced_with_motif_ChIPSeq[[cell_line]]), 
                                end(peaks_reduced_with_motif_ChIPSeq[[cell_line]]),
                                peaks_reduced_with_motif_ChIPSeq[[cell_line]]$coverage)),
              paste0("~/ctcf_motif_",experiment,"_",cell_line, ".txt"),
              sep="\t",row.names=FALSE, quote = FALSE)
}

### CHIPSEQ ###
# UPLOAD THE CHIPSEQ PEAKS
peak_chipseq_gm38=read_peaks("../out_chiapipe/chipseq/GM19238/GM19238_CTCF_dedup.broadPeak")
peak_chipseq_gm38=subsetByOverlaps(peak_chipseq_gm38, CTCF.hg38)
peak_chipseq_gm38=IRanges::reduce(peak_chipseq_gm38+1000)

peak_chipseq_gm39=read_peaks("../out_chiapipe/chipseq/GM19239/GM19239_CTCF_dedup.broadPeak")
peak_chipseq_gm39=subsetByOverlaps(peak_chipseq_gm39, CTCF.hg38)
peak_chipseq_gm39=IRanges::reduce(peak_chipseq_gm39+1000)

peak_chipseq_gm40=read_peaks("../out_chiapipe/chipseq/GM19240/GM19240_CTCF_dedup.broadPeak")
peak_chipseq_gm40=subsetByOverlaps(peak_chipseq_gm40, CTCF.hg38)
peak_chipseq_gm40=IRanges::reduce(peak_chipseq_gm40+1000)

peaks_reduced_with_motif_ChIPSeq=GRangesList("GM19238"=peak_chipseq_gm38, "GM19239"=peak_chipseq_gm39, "GM19240"=peak_chipseq_gm40)

# ADD COVERAGE TRACK

reads_count_chipseq_GM19238=featureCounts(
  "../out_chiapipe/chipseq/GM19238/GM19238_CTCFall_1_subset_check_samtools_dedup_sorted.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19238), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19238$coverage=unname(reads_count_chipseq_GM19238$counts[,1])

reads_count_chipseq_GM19239=featureCounts(
  "../out_chiapipe/chipseq/GM19239/GM19239_CTCFall_1_subset_check_samtools_dedup_sorted.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19239), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19239$coverage=unname(reads_count_chipseq_GM19239$counts[,1])

reads_count_chipseq_GM19240=featureCounts(
  "../out_chiapipe/chipseq/GM19240/GM19240_CTCFall_1_subset_check_samtools_dedup_sorted.bam", 
  annot.ext = convert_peak_to_df_featcnt(peaks_reduced_with_motif_ChIPSeq$GM19240), 
  isPairedEnd = T, 
  nthreads = 100)
peaks_reduced_with_motif_ChIPSeq$GM19240$coverage=unname(reads_count_chipseq_GM19240$counts[,1])

trios=c("GM19238", "GM19239", "GM19240")
experiment="chipseq"
for (cell_line in trios){
  write.table(unname(data.frame(seqnames(peaks_reduced_with_motif_ChIPSeq[[cell_line]]), 
                                start(peaks_reduced_with_motif_ChIPSeq[[cell_line]]), 
                                end(peaks_reduced_with_motif_ChIPSeq[[cell_line]]),
                                peaks_reduced_with_motif_ChIPSeq[[cell_line]]$coverage)),
              paste0("~/ctcf_motif_",experiment,"_",cell_line, ".txt"),
              sep="\t",row.names=FALSE, quote = FALSE)
}


save.image("~/ctcf_motif_chipseq.RData")

print("end")











