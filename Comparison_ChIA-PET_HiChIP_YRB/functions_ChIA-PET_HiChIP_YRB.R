"
Here you can find all important function that were invented for the analysis 
between ChIA-PET and HiChIP
"

suppressMessages(library(stringr))
suppressMessages(library(ChIPseeker))
suppressMessages(library(GenomicRanges))
suppressMessages(library(VennDiagram))
suppressMessages(library(gridExtra))
suppressMessages(library(Rsubread))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(magrittr))
suppressMessages(library(hrbrthemes))
suppressMessages(library(ggridges))
suppressMessages(library(data.table))
suppressMessages(library(diagram))
suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(IRanges))
suppressMessages(library(patchwork))
suppressMessages(library(caret))
suppressMessages(library(MAnorm2))
suppressMessages(library(stringr))

get_peaks_path <- function(exp_type, cell_line){
  dir_path=file.path("data", exp_type, cell_line)
  peaks_file_path=list.files(dir_path, pattern='.broadPeak', 
                             recursive=FALSE, full.names = TRUE)
  return(peaks_file_path)
}


read_bam_file <- function(exp_type, cell_line){
  dir_path=file.path("data", exp_type, cell_line)
  bam_file_path=list.files(dir_path, pattern='\\.for.BROWSER.bam$', 
                           recursive=FALSE, full.names = TRUE)
  return(bam_file_path)
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

#Calculate the FPKM value
calculate_FPKM_TPM=function(peaks){
  RC=peaks$coverage
  SF=1000000
  PL=width(peaks)
  RPM=RC*SF/sum(RC)
  FPKM=RPM/PL*1000
  
  #Calculate the TPM value
  RPK=RC/PL*1000
  TPM=RPK*SF/sum(RPK)
  
  return(data.frame(FPKM=FPKM, TPM=TPM))
}


read_peaks=function(peaks_path){
  peaks.tb <- read.table(peaks_path)
  peaks.gr <- GRanges(seqnames = peaks.tb$V1, IRanges(
    start = peaks.tb$V2, end = peaks.tb$V3), strand = '*',
    pval=10**-peaks.tb$V8, qval = 10^-peaks.tb$V9)
  peaks.gr=sort(peaks.gr)
  return(peaks.gr)
}

find_all_relative_overlaps=function(peaks1, peaks2, peaks3, relative_experiment){
  common_overlaps=find_main_overlaps(peaks1, peaks2, peaks3)
  
  if (relative_experiment=="chiapet"){
    reduced_peaks=peaks1
    
    #Overlaps between clean experiment and common between 3 technologies
    overlapIDs.123=unique(queryHits(findOverlaps(reduced_peaks,
                                                 common_overlaps$com123)))
    
    #Overlaps between clean experiment and common between 1 and 3 technology
    # 7 region
    overlap7=subsetByOverlaps(common_overlaps$com13, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.7=unique(queryHits(findOverlaps(reduced_peaks, overlap7)))
    
    #Overlaps between clean experiment and common between 1 and 2 technology
    # region 8
    overlap8=subsetByOverlaps(common_overlaps$com12, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.8=unique(queryHits(findOverlaps(reduced_peaks, overlap8)))
    
    #Overlapping ids of unique peaks of first experiment
    overlap1=subsetByOverlaps(reduced_peaks, common_overlaps$com12, invert = TRUE)
    overlap_raw1=subsetByOverlaps(overlap1, overlap34, invert = TRUE)
    overlapIDs.raw.1=unique(queryHits(findOverlaps(reduced_peaks, overlap_raw1)))
    
    chiapet_overlaps=(GRangesList(vs_all=overlapIDs.123, vs_ChIPSeq=overlapIDs.7, vs_HiChIP=overlapIDs.8, raw=overlapIDs.raw.1))
  }
  else if (relative_experiment=="hichip"){
    reduced_peaks=peaks2
    
    #Overlaps between clean experiment and common between 3 technologies
    overlapIDs.123=unique(queryHits(findOverlaps(reduced_peaks,
                                                 common_overlaps$com123)))
    
    #Overlaps between clean experiment and common between 2 and 3 technology
    # 9 region
    overlap9=subsetByOverlaps(common_overlaps$com23, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.9=unique(queryHits(findOverlaps(reduced_peaks, overlap9)))
    
    #Overlaps between clean experiment and common between 1 and 2 technology
    # region 8
    overlap8=subsetByOverlaps(common_overlaps$com12, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.8=unique(queryHits(findOverlaps(reduced_peaks, overlap8)))
    
    # #Overlapping ids of unique peaks of second experiment
    overlap2=subsetByOverlaps(reduced_peaks, common_overlaps$com12, invert = TRUE)
    overlap_raw2=subsetByOverlaps(overlap2, overlap24, invert = TRUE)
    overlapIDs.raw.2=unique(queryHits(findOverlaps(reduced_peaks, overlap_raw2)))
    
    return(GRangesList(vs_all=overlapIDs.123, vs_ChIPSeq=overlapIDs.9, vs_ChIAPET=overlapIDs.8, raw=overlapIDs.raw.2))
  }
  else if (relative_experiment=="chipseq"){
    reduced_peaks=peaks3
    
    #Overlaps between clean experiment and common between 3 technologies
    overlapIDs.123=unique(queryHits(findOverlaps(reduced_peaks,
                                                 common_overlaps$com123)))
    
    #Overlaps between clean experiment and common between 2 and 3 technology
    # 9 region
    overlap9=subsetByOverlaps(common_overlaps$com23, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.9=unique(queryHits(findOverlaps(reduced_peaks, overlap9)))
    
    #Overlaps between clean experiment and common between 1 and 3 technology
    # 7 region
    overlap7=subsetByOverlaps(common_overlaps$com13, common_overlaps$com123,
                              invert = TRUE)
    overlapIDs.7=unique(queryHits(findOverlaps(reduced_peaks, overlap7)))
    
    # #Overlapping ids of unique peaks of third experiment
    overlap3=subsetByOverlaps(reduced_peaks, common_overlaps$com23, invert = TRUE)
    overlap_raw3=subsetByOverlaps(overlap3, overlap34, invert = TRUE)
    overlapIDs.raw.3=unique(queryHits(findOverlaps(reduced_peaks, overlap_raw3)))
    
    return(GRangesList(vs_all=overlapIDs.123, vs_ChIAPET=overlapIDs.7, vs_HiChIP=overlapIDs.9, raw=overlapIDs.raw.3))
  }
}

# Function that looks for all overlaps in technology, in final peaks are not overlaping with each other. 

find_true_all_relative_overlaps=function(peaks1, peaks2, peaks3, rel_param){
  
  peaks1$common_peaks2=1:NROW(peaks1) %in% unique(queryHits(findOverlaps(peaks1, peaks2)))
  peaks1$common_peaks3=1:NROW(peaks1) %in% unique(queryHits(findOverlaps(peaks1, peaks3)))
  
  peaks2$common_peaks1=1:NROW(peaks2) %in% unique(queryHits(findOverlaps(peaks2, peaks1)))
  peaks2$common_peaks3=1:NROW(peaks2) %in% unique(queryHits(findOverlaps(peaks2, peaks3)))
  
  peaks3$common_peaks1=1:NROW(peaks3) %in% unique(queryHits(findOverlaps(peaks3, peaks1)))
  peaks3$common_peaks2=1:NROW(peaks3) %in% unique(queryHits(findOverlaps(peaks3, peaks2)))
  
  id_common_peaks_1=id_common_peaks_2=id_common_peaks_3=
    peaks1_common_peaks2_id=peaks1_common_peaks3_id=peaks2_common_peaks1_id=
    peaks2_common_peaks3_id=peaks3_common_peaks1_id=peaks3_common_peaks2_id=
    peaks1_unique_id=peaks2_unique_id=peaks3_unique_id=c()
  
  id_common_peaks_1=which(peaks1$common_peaks2 == TRUE & peaks1$common_peaks3 == TRUE)
  peaks1_common_peaks2_id=which(peaks1$common_peaks2 == TRUE & peaks1$common_peaks3 == FALSE)
  peaks1_common_peaks3_id=which(peaks1$common_peaks2 == FALSE & peaks1$common_peaks3 == TRUE)
  peaks1_unique_id=which(peaks1$common_peaks2 == FALSE & peaks1$common_peaks3 == FALSE)
  
  id_common_peaks_2=which(peaks2$common_peaks1 == TRUE & peaks2$common_peaks3 == TRUE)
  peaks2_common_peaks1_id=which(peaks2$common_peaks1 == TRUE & peaks2$common_peaks3 == FALSE)
  peaks2_common_peaks3_id=which(peaks2$common_peaks1 == FALSE & peaks2$common_peaks3 == TRUE)
  peaks2_unique_id=which(peaks2$common_peaks1 == FALSE & peaks2$common_peaks3 == FALSE)
  
  id_common_peaks_3=which(peaks3$common_peaks1 == TRUE & peaks3$common_peaks2 == TRUE)
  peaks3_common_peaks1_id=which(peaks3$common_peaks1 == TRUE & peaks3$common_peaks2 == FALSE)
  peaks3_common_peaks2_id=which(peaks3$common_peaks1 == FALSE & peaks3$common_peaks2 == TRUE)
  peaks3_unique_id=which(peaks3$common_peaks1 == FALSE & peaks3$common_peaks2 == FALSE)
  
  constant_peaks=sort(c(peaks1[id_common_peaks_1], peaks2[id_common_peaks_2], peaks3[id_common_peaks_3]))
  constant_peaks=IRanges::reduce(constant_peaks)
  
  possible_common_peaks_1_2=IRanges::reduce(c(peaks1[peaks1_common_peaks2_id], peaks2[peaks2_common_peaks1_id]))
  possible_common_peaks_1_3=IRanges::reduce(c(peaks1[peaks1_common_peaks3_id], peaks3[peaks3_common_peaks1_id]))
  possible_common_peaks_2_3=IRanges::reduce(c(peaks2[peaks2_common_peaks3_id], peaks3[peaks3_common_peaks2_id]))
  
  possible_unique_peaks1=IRanges::reduce(peaks1[peaks1_unique_id])
  possible_unique_peaks2=IRanges::reduce(peaks2[peaks2_unique_id])
  possible_unique_peaks3=IRanges::reduce(peaks3[peaks3_unique_id])
  
  unique_common_peaks_1_2=IRanges::reduce(subsetByOverlaps(possible_common_peaks_1_2, constant_peaks, invert=TRUE))
  unique_common_peaks_1_3=IRanges::reduce(subsetByOverlaps(possible_common_peaks_1_3, constant_peaks, invert=TRUE))
  unique_common_peaks_2_3=IRanges::reduce(subsetByOverlaps(possible_common_peaks_2_3, constant_peaks, invert=TRUE))
  
  unique_peaks_1=IRanges::reduce(subsetByOverlaps(possible_unique_peaks1, c(unique_common_peaks_1_2, unique_common_peaks_1_3, constant_peaks), invert=TRUE))
  unique_peaks_2=IRanges::reduce(subsetByOverlaps(possible_unique_peaks2, c(unique_common_peaks_1_2, unique_common_peaks_2_3, constant_peaks), invert=TRUE))
  unique_peaks_3=IRanges::reduce(subsetByOverlaps(possible_unique_peaks3, c(unique_common_peaks_1_3, unique_common_peaks_2_3, constant_peaks), invert=TRUE))
  
  if (rel_param == "trios"){
    true_common_peaks=GRangesList(constant_peaks=constant_peaks, 
                                  common_38_39=unique_common_peaks_1_2,
                                  common_38_40=unique_common_peaks_1_3,
                                  common_39_40=unique_common_peaks_2_3,
                                  unique_38=unique_peaks_1,
                                  unique_39=unique_peaks_2,
                                  unique_40=unique_peaks_3)
  }
  else if (rel_param == "experiments"){
    true_common_peaks=GRangesList(all_technologies=constant_peaks, 
                                  common_ChIAPET_HiChIP=unique_common_peaks_1_2,
                                  common_ChIAPET_ChIPSeq=unique_common_peaks_1_3,
                                  common_HiChIP_ChIPSeq=unique_common_peaks_2_3,
                                  unique_ChIAPET=unique_peaks_1,
                                  unique_HiChIP=unique_peaks_2,
                                  unique_ChIPSeq=unique_peaks_3)
  }
  else {
    true_common_peaks=GRangesList(constant_peaks=constant_peaks, 
                                  unique_common_peaks_1_2=unique_common_peaks_1_2,
                                  unique_common_peaks_1_3=unique_common_peaks_1_3,
                                  unique_common_peaks_2_3=unique_common_peaks_2_3,
                                  unique_peaks_1=unique_peaks_1,
                                  unique_peaks_2=unique_peaks_2,
                                  unique_peaks_3=unique_peaks_3)
  }
  return(true_common_peaks)
}

get_percentg_common_peaks=function(grangesList){
  
  nrows_grangeslist=elementNROWS(grangesList)
  n_all_peaks=sum(nrows_grangeslist)
  peaks_experiments=list(chiapet=sum(nrows_grangeslist[c(1,2,3,5)]),
                         hichip=sum(nrows_grangeslist[c(1,2,4,6)]),
                         chipseq=sum(nrows_grangeslist[c(1,3,4,7)]))
  
  no_peaks=c(n_all_peaks, 
             nrows_grangeslist[c(1,3,4,2,2,5,6,7)])
  
  reference=c(n_all_peaks/n_all_peaks, 
              nrows_grangeslist[1]/n_all_peaks,
              sum(nrows_grangeslist[c(1,3)])/peaks_experiments$chiapet, 
              sum(nrows_grangeslist[c(1,4)])/peaks_experiments$hichip, 
              sum(nrows_grangeslist[c(1,2)])/peaks_experiments$chiapet,
              sum(nrows_grangeslist[c(1,2)])/peaks_experiments$hichip,
              nrows_grangeslist[5]/peaks_experiments$chiapet, 
              nrows_grangeslist[6]/peaks_experiments$hichip, 
              nrows_grangeslist[7]/peaks_experiments$chipseq)
  
  
  percent_common_peaks=data.frame(no_peaks, 
                                  percentage=round(reference*100, 2))
  
  
  rownames(percent_common_peaks) = c("no_peaks", 
                                     "3 technologies", "chiapet_chipseq", "hichip_chipseq", 
                                     "chiapet_hichip","hichip_chiapet", 
                                     "unique_chiapet", "unique_hichip", "unique_chipseq")
  return(percent_common_peaks)
}

# CTCF motif

add_motif_info_to_peaks=function(peaks, CTCF_data){
  
  peaks.ctcfoverlap = findOverlaps(peaks, CTCF_data)
  
  numbers_CTCF=c()
  max_motif_CTCF=c()
  average_motif_CTCF=c()
  v=1
  
  for (p in 1:NROW(peaks)){
    id_motifs_CTCF_per_peak=c()
    while (v <= NROW(peaks.ctcfoverlap) && queryHits(peaks.ctcfoverlap)[v] == p){
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

generate_manorm_plot=function(cells, norm, raw=NULL){
  if(!is.null(raw)){
    param=list(raw[[paste0(cells[1], ".read_cnt")]], 
               raw[[paste0(cells[2], ".read_cnt")]], 
               norm[[paste0(cells[1], ".occupancy")]], 
               norm[[paste0(cells[2], ".occupancy")]])
    main="Before normalization"
  } else {
    param=list(norm[[paste0(cells[1], ".read_cnt")]], 
               norm[[paste0(cells[2], ".read_cnt")]], 
               norm[[paste0(cells[1], ".occupancy")]], 
               norm[[paste0(cells[2], ".occupancy")]])
    main="After normalization"
  }
  MAplot(param[[2]], param[[1]], param[[4]], param[[3]], 
         ylim = set_ylim, xlim=set_xlim,
         main=main, cex=1.5, cex.main=2.3, 
         line=0.2,
         cex.axis=2, cex.lab=2,
         args.legend = list(
           x="topright",
           inset = c(0, 0),
           legend=c("common", paste0(cells[1]," specific"), paste0(cells[2]," specific"), "others"),
           bty = "n", xpd = TRUE,
           y.intersp=0.8,
           x.intersp=1,
           cex=2))
  abline(h = 0, lwd = 2, lty = 5)
  
}

generate_manorm_plot_trios=function(cells, norm, type){
  if(type=="before"){main="Before normalization"} 
  else {main="After normalization"}
  
  MAplot(norm[[cells[2]]], norm[[cells[1]]], ylim = set_ylim, xlim=set_xlim,
         main=main, cex=1.5, cex.main=2.3, line=0.2, 
         cex.axis=2, cex.lab=2,
         # args.legend = list(
         #   x="topright",
         #   # inset = c(0, 0),
         #   legend=legend,
         #   bty = "o", xpd = TRUE,
         #   y.intersp=0.8,
         #   # x.intersp=1,
         #   cex=1.5)
  )
  abline(h = 0, lwd = 2, lty = 5)
}

create_table_common_peaks_with_info=function(norm, differential, conds, cell_lines){
  common=data.frame(chrom=norm$chrom,
                    start=norm$start,
                    end=norm$end,
                    conds[[cell_lines[1]]]$occupancy,
                    conds[[cell_lines[2]]]$occupancy,
                    differential[[paste0(cell_lines[1],"_",cell_lines[2])]],
                    conds[[cell_lines[1]]]$norm.signal,
                    conds[[cell_lines[2]]]$norm.signal,
                    norm$index)
  colnames(common)[[4]]=paste(cell_lines[1],'occupancy', sep='.')
  colnames(common)[[5]]=paste(cell_lines[2],'occupancy', sep='.')
  
  common$isPadj=ifelse(common$padj<0.001, TRUE, FALSE)
  common$isPval=ifelse(common$pval<p_value, TRUE, FALSE)
  common$common_padj=ifelse(common$padj>=0.001, TRUE, FALSE)
  common$common_pval=ifelse(common$pval>=p_value, TRUE, FALSE)
  common$predicted_rev=ifelse(common[[paste0(cell_lines[1],".occupancy")]]==TRUE & 
                                common[[paste0(cell_lines[2],".occupancy")]]==TRUE , TRUE, FALSE)
  return(common)
}

create_fp_tables=function(common,cell_lines){
  # Devide the picks into different groups
  cell1=paste0(str_sub(cell_lines[1],1,2),str_sub(cell_lines[1],start=6))
  cell2=paste0(str_sub(cell_lines[2],1,2),str_sub(cell_lines[2],start=6))
  
  cell1.occupancy=paste0(cell_lines[1],".occupancy")
  cell2.occupancy=paste0(cell_lines[2],".occupancy")
  cell1.mean=paste0(cell_lines[1],".mean")
  cell2.mean=paste0(cell_lines[2],".mean")
  
  tf_table=list(common[common$predicted_rev==TRUE,],
                common[common$predicted_rev==FALSE,],
                common[common[[cell1.occupancy]]==FALSE & 
                         common[[cell2.occupancy]]==FALSE ,],
                common[common[[cell1.occupancy]]==TRUE & 
                         common[[cell2.occupancy]]==FALSE ,],
                common[common[[cell1.occupancy]]==FALSE & 
                         common[[cell2.occupancy]]==TRUE ,])
  names(tf_table)=c(paste0(cell1,"T_",cell2,"T"),
                    "false_peaks",
                    paste0(cell1,"F_",cell2,"F"),
                    paste0(cell1,"T_",cell2,"F"),
                    paste0(cell1,"F_",cell2,"T"))
  
  return(tf_table)
}

save_fp_graphs_with_mean_values=function(tf_table,cell_lines, experiment, prefix, padj_or_pval="isPadj"){
  png(paste0("plots/mean_values/",experiment, "_", cell1,"_",cell2,"_diff_test_",prefix,".png"), res = 300, width = 3500, height = 2200)
  plot.new()
  p1=ggplot(tf_table[[1]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("Common peaks in 2 cell lines")+labs(x=cell1.mean, y=cell2.mean)
  p2=ggplot(tf_table[[2]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("Uncommon peaks in 2 cell lines")+labs(x=cell1.mean, y=cell2.mean)
  p3=ggplot(tf_table[[3]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle('No peaks\nin both cell lines')+theme(legend.position = "none")+labs(x=cell1.mean, y=cell2.mean)
  p4=ggplot(tf_table[[4]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle(paste0('Peaks only\nin ', cell_lines[1]))+theme(legend.position = "none")+labs(x=cell1.mean, y=cell2.mean)
  p5=ggplot(tf_table[[5]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle(paste0('Peaks only\nin ',cell_lines[2]))+labs(x=cell1.mean, y=cell2.mean)
  p6=ggplot(common, aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("All peaks")+geom_smooth(method="lm", formula=y~x)+labs(x=cell1.mean, y=cell2.mean)
  print(p6+(p1 | p2) / (p3 | p4 | p5) + plot_layout(widths=c(1, 2)))
  dev.off()
}

add_greater_and_lower_mean_ranges_tables=function(common_manorm, mean_cutoff, cell_lines){
  comp_cell_lines=paste0(paste0(str_sub(cell_lines[1],1,2),str_sub(cell_lines[1],start=6)),"_",paste0(str_sub(cell_lines[2],1,2),str_sub(cell_lines[2],start=6)))
  
  i=length(common_manorm)
  
  common_manorm[[i+1]]=common_manorm[[comp_cell_lines]][common_manorm[[comp_cell_lines]][[paste0(cell_lines[1],".mean")]]>=mean_cutoff&                                                   
                                                          common_manorm[[comp_cell_lines]][[paste0(cell_lines[2],".mean")]]>=mean_cutoff,]
  
  common_manorm[[i+2]]=common_manorm[[comp_cell_lines]][common_manorm[[comp_cell_lines]][[paste0(cell_lines[1],".mean")]]<mean_cutoff&                                                                        
                                                          common_manorm[[comp_cell_lines]][[paste0(cell_lines[2],".mean")]]<mean_cutoff,]
  
  names(common_manorm)[(i+1):(i+2)]=c(paste0(comp_cell_lines, "_mean_gr10"),
                                      paste0(comp_cell_lines, "_mean_lw10"))
  
  return(common_manorm)
}


create_short_name_for_col=function(trios){
  names=c(
    paste0(paste0(str_sub(trios[1],1,2),str_sub(trios[1],start=6)),"_",paste0(str_sub(trios[2],1,2),str_sub(trios[2],start=6))),
    paste0(paste0(str_sub(trios[1],1,2),str_sub(trios[1],start=6)),"_",paste0(str_sub(trios[3],1,2),str_sub(trios[3],start=6))),
    paste0(paste0(str_sub(trios[2],1,2),str_sub(trios[2],start=6)),"_",paste0(str_sub(trios[3],1,2),str_sub(trios[3],start=6))))
  return(names)
}


#For Venn Diagrams

color=function(word){
  if (word == "GM19238") {color="deeppink"}
  else if (word == "GM19239") {color="deepskyblue"}
  else if (word == "GM19240") {color="orange"}
  else if (word == "chiapet") {color="blue"}
  else if (word == "hichip") {color="lawngreen"}
  else if (word == "chipseq") {color="red"}
  
  else if(!is.na(str_extract(word, "rep1"))){
    if (experiment=="chiapet"){color="cornflowerblue"}
    else if (experiment=="hichip") {color="greenyellow"}
  }
  else if(!is.na(str_extract(word, "rep2"))){
    if (experiment=="chiapet"){color="darkblue"}
    else if (experiment=="hichip") {color="green4"}
  }
  return(color)
}

get_color_all_overlaps=function(word, cell_line){
  if (word=="all_trios") {color="#ccd47c"}
  else {
    if (cell_line == "GM19238") {
      if (word=="unique") {color="#ff91cd"}
      else if (word=="GM19239") {color="#98aaea"}
      else if (word=="GM19240") {color="#fac868"}
    }
    if (cell_line == "GM19239") {
      if (word=="unique") {color="deeppink"}
      else if (word=="GM19238") {color="#98aaea"}
      else if (word=="GM19240") {color="#c6ed8c"}
    }
    if (cell_line == "GM19240") {
      if (word=="unique") {color="yellow"}
      else if (word=="GM19238") {color="#fac868"}
      else if (word=="GM19239") {color="#c6ed8c"}
    }
  }
  return(color)
}

create_pie_plot_my=function(slices.df, cell_line){
  df2 <- slices.df %>% 
    mutate(perc = quantity/ sum(quantity)) %>% 
    mutate(perc = scales::percent(perc)) %>% 
    arrange(desc(labels)) %>% ## arrange in the order of the legend
    mutate(pos = cumsum(quantity) - quantity/1.35) ### calculate where to place the text labels  
  
  
  ggplot(slices.df, aes(x = "" , y = quantity, fill = labels)) +
    geom_col(width = 0.8, color = 1) +
    coord_polar(theta = "y", start=5)  +
    geom_label_repel(data = df2,
                     aes(y = pos, label = paste0(quantity, "\n", round(quantity/sum(quantity)*100, 2), "%")),
                     size = 4.5, nudge_x = 1, show.legend = FALSE) +
    scale_fill_manual(values=c(get_color_all_overlaps(slices.df$labels[[1]], cell_line),
                               get_color_all_overlaps(slices.df$labels[[2]], cell_line),
                               get_color_all_overlaps(slices.df$labels[[3]], cell_line),
                               get_color_all_overlaps(slices.df$labels[[4]], cell_line))) +
    theme_void()
}  

prop_name=function(technology){
  if (technology == "chiapet") {name="ChIA-PET"}
  else if (technology == "hichip") {name="Hi-ChIP"}
  else if (technology == "chipseq") {name="ChIP-Seq"}
  
  else if(!is.na(str_extract(technology, "rep1"))){
    if (experiment=="chiapet"){name="ChIA-PET_rep1"}
    else if (experiment=="hichip") {name="Hi-ChIP_rep2"}}
  else if(!is.na(str_extract(technology, "rep2"))){
    if (experiment=="chiapet"){name="ChIA-PET_rep2"}
    else if (experiment=="hichip") {name="Hi-ChIP_rep2"}}
  
  else {name=technology}
  
  return(name)
}


draw_2_venn_diagram_from_raw=function(area1, area2, cross.area, region=12){
  
  if (area1 < area2){
    rotation_degree=180
    cat_pos=c(60,300)
    ext_pos=c(140, -140)
  }
  else {
    rotation_degree=0
    cat_pos=c(240,120)
    ext_pos=c(-30, 30)
  }
  
  if (region==12) {exp_1=experiment_1 
  exp_2=experiment_2}
  else if (region==13){exp_1=experiment_1 
  exp_2=experiment_3}
  else if (region==23){exp_1=experiment_2
  exp_2=experiment_3}
  
  if (experiment_1==experiment_2){
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
                                              scaled = T,
                                              category=c(run_1, run_2),
                                              print.mode=c('raw','percent'),
                                              fill=c(color(run_1), color(run_2)),
                                              rotation.degree=rotation_degree,
                                              lwd = 3,
                                              cat.cex = 2,
                                              cat.pos = cat_pos,
                                              cat.dist = c(0.11, 0.11),
                                              ext.pos = ext_pos,
                                              ext.percent=.25,
                                              margin = 0.11,
                                              cex = 2.5)
  } else {
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
                                              scaled = T,
                                              category=c(prop_name(exp_1), prop_name(exp_2)),
                                              print.mode=c('raw'),
                                              fill=c(color(exp_1), color(exp_2)),
                                              rotation.degree=rotation_degree,
                                              lwd = 3,
                                              cat.cex = 2,
                                              cat.pos = cat_pos,
                                              cat.dist = c(0.11, 0.11),
                                              ext.pos = ext_pos,
                                              ext.percent=.25,
                                              margin = 0.11,
                                              cex = 2.5)
    
  }
}


draw_2_venn_diagram=function(common_peaks, region=12){
  peaks_quantity=elementNROWS(common_peaks)
  
  if (region==12) {
    cross.area=peaks_quantity[[2]]+peaks_quantity[["constant_peaks"]]
    area1=cross.area+peaks_quantity[[3]]+peaks_quantity[[5]]
    area2=cross.area+peaks_quantity[[4]]+peaks_quantity[[6]]
  }
  else if (region==13){
    cross.area=peaks_quantity[3]+peaks_quantity[["constant_peaks"]]
    area1=cross.area+peaks_quantity[2]+peaks_quantity[5]
    area2=cross.area+peaks_quantity[4]+peaks_quantity[7]
  }
  else if (region==23){
    cross.area=peaks_quantity[4]+peaks_quantity[["constant_peaks"]]
    area1=cross.area+peaks_quantity[2]+peaks_quantity[6]
    area2=cross.area+peaks_quantity[3]+peaks_quantity[7]
  }
  
  if (area1 < area2){
    rotation_degree=180
    cat_pos=c(60,300)
    ext_pos=c(140, -140)
  }
  else {
    rotation_degree=0
    cat_pos=c(240,120)
    ext_pos=c(-30, 30)
  }
  
  if (region==12) {exp_1=experiment_1 
  exp_2=experiment_2}
  else if (region==13){exp_1=experiment_1 
  exp_2=experiment_3}
  else if (region==23){exp_1=experiment_2
  exp_2=experiment_3}
  
  if (experiment_1==experiment_2){
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
                                              scaled = T,
                                              category=c(run_1, run_2),
                                              print.mode=c('raw','percent'),
                                              fill=c(color(run_1), color(run_2)),
                                              rotation.degree=rotation_degree,
                                              lwd = 3,
                                              cat.cex = 2,
                                              cat.pos = cat_pos,
                                              cat.dist = c(0.11, 0.11),
                                              ext.pos = ext_pos,
                                              ext.percent=.25,
                                              margin = 0.11,
                                              cex = 2.5)
    # grid.arrange(gTree(children = venndiag),
    #              top = textGrob( paste0("Peaks comparison between the replicates for ", prop_name(experiment_1)), gp=gpar(cex=4)))
  } else {
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
                                              scaled = T,
                                              category=c(prop_name(exp_1), prop_name(exp_2)),
                                              print.mode=c('raw'),
                                              fill=c(color(exp_1), color(exp_2)),
                                              rotation.degree=rotation_degree,
                                              lwd = 3,
                                              cat.cex = 2,
                                              cat.pos = cat_pos,
                                              cat.dist = c(0.11, 0.11),
                                              ext.pos = ext_pos,
                                              ext.percent=.25,
                                              margin = 0.11,
                                              cex = 2.5)
    
    # venndiag.ga=grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between technologies for ", cell1))
  }
}

draw_3_venn_diagram=function(common_peaks){
  
  areas=define_subsets_for_3_venn_from_common(common_peaks)
  
  # Define the required overlaps for Venn Diagram
  peaks_quantity=elementNROWS(common_peaks)
  n123=peaks_quantity[["constant_peaks"]]
  n12=peaks_quantity[[2]]+n123
  n13=peaks_quantity[[3]]+n123
  n23=peaks_quantity[[4]]+n123
  
  # Calculate the peaks quantity from the overlaps
  area1=n12+peaks_quantity[[3]]+peaks_quantity[[5]]
  area2=n12+peaks_quantity[[4]]+peaks_quantity[[6]]
  area3=n13+peaks_quantity[[4]]+peaks_quantity[[7]]
  
  #Draw plots
  if (experiment_1 == experiment_2 & experiment_1 == experiment_3){
    venndiag = VennDiagram::draw.triple.venn(areas$area1, areas$area2, areas$area3, 
                                             areas$n12, areas$n23, areas$n13, areas$n123, 
                                             fill=c(color(trios[1]), color(trios[2]), color(trios[3])),
                                             alpha=c(0.1, 0.1, 0.1),
                                             category=c(trios[1], trios[2], trios[3]),
                                             print.mode=c('raw','percent'),
                                             lwd = 4,
                                             cat.dist = c(0.06, 0.06, 0.03),
                                             cat.cex = 3, 
                                             margin = 0.12,
                                             cex = 3)
    # grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between tiosFamily for ", prop_name(experiment_1)))
  } else {
    venndaig=VennDiagram::draw.triple.venn(areas$area1, areas$area2, areas$area3, 
                                           areas$n12, areas$n23, areas$n13, areas$n123,
                                           fill=c(color(experiment_1), color(experiment_2), color(experiment_3)),
                                           category=c(prop_name(experiment_1), prop_name(experiment_2), prop_name(experiment_3)),
                                           print.mode=c('raw', 'percent'),
                                           lwd = 4,
                                           cat.dist = c(0.06, 0.06, 0.03),
                                           cat.cex = 3, 
                                           margin = 0.12,
                                           cex = 3)
    # venndiag.ga=grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between technologies for ", cell_line_1))
  }
}

define_subsets_for_3_venn_from_common=function(common_peaks){
  
  # Define the required overlaps for Venn Diagram
  peaks_quantity=elementNROWS(common_peaks)
  n123=peaks_quantity[["all_technologies"]]
  n12=peaks_quantity[[2]]+n123
  n13=peaks_quantity[[3]]+n123
  n23=peaks_quantity[[4]]+n123
  
  # Calculate the peaks quantity from the overlaps
  area1=n12+peaks_quantity[[3]]+peaks_quantity[[5]]
  area2=n12+peaks_quantity[[4]]+peaks_quantity[[6]]
  area3=n13+peaks_quantity[[4]]+peaks_quantity[[7]]
  
  areas=list(n123=n123,
             n12=n12,
             n13=n13,
             n23=n23,
             area1=area1,
             area2=area2,
             area3=area3)
  return(areas)
}

clean_all_loops = function(bedp, peaks, score.threshold, CTCF.motifs, maxwidth=1000000) {
  # convert bedp to grange
  anchor1.grange = GRanges(seqnames = bedp$V1, IRanges(start = bedp$V2, end = bedp$V3), strand = '*')
  anchor2.grange = GRanges(seqnames = bedp$V4, IRanges(start = bedp$V5, end = bedp$V6), strand = '*')
  scores= bedp$V7
  
  # Convert peaks to grange
  peaks.grange = GRanges(seqnames = peaks$V1, IRanges(start = peaks$V2, end = peaks$V3), strand = '*', minus_log10pval=peaks$V7)
  print(NROW(peaks.grange))
  
  "Next filters will be applied:
  0) They should have ctcf motif
  1) Support in two of the anchors
  2) min 3 PETS
  3) motifs are from the ENCODE
  " 
  
  #0 peaks with motif
  peak.with.motif = findOverlaps(peaks.grange, CTCF.motifs)
  peak.with.motif = peaks.grange[unique(queryHits(peak.with.motif)), ]
  print(NROW(peak.with.motif))
  
  
  #1 Loop with 2 anchors
  anchor1.ovelap.peaks = findOverlaps(anchor1.grange, peak.with.motif)
  anchor2.ovelap.peaks = findOverlaps(anchor2.grange, peak.with.motif)
  
  both.overlap = unique(queryHits(anchor1.ovelap.peaks)[queryHits(anchor1.ovelap.peaks) %in% queryHits(anchor2.ovelap.peaks)])
  
  anchor1.peaksupport.withmotif = anchor1.grange[both.overlap, ]
  anchor2.peaksupport.withmotif = anchor2.grange[both.overlap, ]
  scores.peaksupport.withmotif = scores[both.overlap]
  
  #3 Peaks support + CTCf motif (if there are not already just comment)
  #anchor1.motif = findOverlaps(anchor1.peaksupport, CTCF.motifs)
  #anchor2.motif = findOverlaps(anchor2.peaksupport, CTCF.motifs)
  
  #both.motif = unique(queryHits(anchor1.motif)[queryHits(anchor1.motif) %in% queryHits(anchor2.motif)])
  
  #anchor1.peaksupport.wihhmotif = anchor1.peaksupport[both.motif, ]
  #anchor2.peaksupport.withmotif = anchor2.peaksupport[both.motif, ]
  #scores.peaksupport.withmotif = scores.peaksupport[both.motif]
  
  
  #2 Remove peaks with low score 
  above.threshold<- which(scores.peaksupport.withmotif >= score.threshold)
  
  anchor1.peaksupport.wihhmotif.nosingle = anchor1.peaksupport.withmotif[above.threshold, ]
  anchor2.peaksupport.wihhmotif.nosingle = anchor2.peaksupport.withmotif[above.threshold, ]
  scores.peaksupport.wihhmotif.nosingle = scores.peaksupport.withmotif[above.threshold]
  
  
  # Filter on the length
  lengths = end(anchor2.peaksupport.wihhmotif.nosingle)-start(anchor1.peaksupport.wihhmotif.nosingle)
  to_keep = which(lengths < maxwidth)
  
  anchor1.peaksupport.wihhmotif.nosingle.maxwidth = anchor1.peaksupport.wihhmotif.nosingle[to_keep, ]
  anchor2.peaksupport.wihhmotif.nosingle.maxwidth = anchor2.peaksupport.wihhmotif.nosingle[to_keep, ]
  scores.peaksupport.wihhmotif.nosingle.maxwidth = scores.peaksupport.wihhmotif.nosingle[to_keep]
  
  return(list(anchor1.peaksupport.wihhmotif.nosingle.maxwidth, 
              anchor2.peaksupport.wihhmotif.nosingle.maxwidth,
              scores.peaksupport.wihhmotif.nosingle.maxwidth ))
}


common_loops_new = function(lista1, lista2, max_gap = -1) {
  "
  Function search for the ids between the samples. 
  Input:left anchors, right anchors and vector scale
  Output: 1) ccommon first samples loops with second
          2) ccommon seond samples loops with first
          3) relation of the first group with the strength of the second
          4) relation of the second group with the strength of the first
    
  "
  anchor1.a <- lista1[[1]]
  anchor2.a <- lista1[[2]]
  
  anchor1.b <- lista2[[1]]
  anchor2.b <- lista2[[2]]
  
  anchor1.overlap <- findOverlaps(anchor1.a, anchor1.b, maxgap = max_gap)
  anchor2.overlap <- findOverlaps(anchor2.a, anchor2.b, maxgap = max_gap)
  
  common_loops_sample1 <- c()
  common_loops_sample2 <- c()
  
  list_mapping_sample1 <- list()
  list_mapping_sample2 <- list()
  
  # Now there will be no loop without it
  i <- 1
  for(id in unique(queryHits(anchor1.overlap))) {
    # Check if the loops 1 are mapped to loops 2 (right and left anchors)
    loops_sample2_anchor_1 <- subjectHits(anchor1.overlap[queryHits(anchor1.overlap) == id])
    loops_sample2_anchor_2 <- subjectHits(anchor2.overlap[queryHits(anchor2.overlap) == id])                                  
    
    both_mapping <- loops_sample2_anchor_1[loops_sample2_anchor_1 %in% loops_sample2_anchor_2 ]
    if(length(both_mapping) > 0) {
      # Check if loop with id from sample 1 one is present in sample 2
      common_loops_sample1 <- c(common_loops_sample1, id)
      list_mapping_sample1[[i]] <- both_mapping
      i <-  i +1
    } 
    
  } #End of the iteraction on loops 1
  i <- 1
  for(id in unique(subjectHits(anchor1.overlap))) {
    # Check if loops ids from sample 1 map to sample 2 (right and left anchors)
    loopy_sample1_anchor_1 <- queryHits(anchor1.overlap[subjectHits(anchor1.overlap) == id])
    loopy_sample1_anchor_2 <- queryHits(anchor2.overlap[subjectHits(anchor2.overlap) == id])                                  
    
    both_mapping <- loopy_sample1_anchor_1[loopy_sample1_anchor_1 %in% loopy_sample1_anchor_2 ]
    if(length(both_mapping) > 0) {
      # check whether loop with the same id from sample 1 have common with sample 2
      common_loops_sample2 <- c(common_loops_sample2, id)
      list_mapping_sample2[[i]] <- both_mapping
      i <- i + 1
    }
    
  } # End of the iterration
  
  names(list_mapping_sample1) <- common_loops_sample1
  names(list_mapping_sample2) <- common_loops_sample2
  return(list(common_loops_sample1, common_loops_sample2, list_mapping_sample1, 
              list_mapping_sample2))
}
