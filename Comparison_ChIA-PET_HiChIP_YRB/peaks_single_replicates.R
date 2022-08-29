
# Upload peaks for single replicates 

```{r}
experiment<<-"chiapet"
cell_lines=c("GM19238_rep1", "GM19238_rep2", "GM19239_rep1", "GM19239_rep2","GM19240_rep1","GM19240_rep2")

rep2_peaks_raw_ChIAPET=GRangesList("GM19238_rep1"=GRanges(), "GM19238_rep2"=GRanges(),
                                  "GM19239_rep1"=GRanges(),"GM19239_rep2"=GRanges(),
                                  "GM19240_rep1"=GRanges(),"GM19240_rep2"=GRanges())

rep2_peaks_wout_sex_chr_ChIAPET=rep2_peaks_with_motif_ChIAPET=rep2_peaks_reduced_with_motif_ChIAPET=rep2_peaks_raw_ChIAPET

#Create GRangesLists that will conclude all the trios data

for (cell_line in cell_lines){
  
  rep2_peaks_raw_ChIAPET[[cell_line]]=read_peaks(get_peaks_path(experiment, cell_line))
  
  rep2_peaks_with_motif_ChIAPET[[cell_line]]=subsetByOverlaps(rep2_peaks_raw_ChIAPET[[cell_line]],
                                                                              CTCF.hg38)
  
  rep2_peaks_reduced_with_motif_ChIAPET[[cell_line]]=IRanges::reduce(rep2_peaks_with_motif_ChIAPET[[cell_line]]+p_exp)
  
  rep2_peaks_wout_sex_chr_ChIAPET[[cell_line]]=subset(
    rep2_peaks_with_motif_ChIAPET[[cell_line]],
    !(rep2_peaks_with_motif_ChIAPET[[cell_line]]@seqnames %in% c('chrY', 'chrX', 'chrM')))
}

experiment<<-"hichip"
cell_lines=c("GM19238_rep1", "GM19238_rep2", "GM19239_rep1", "GM19239_rep2","GM19240_rep1","GM19240_rep2")

rep2_peaks_raw_HiChIP=GRangesList("GM19238_rep1"=GRanges(), "GM19238_rep2"=GRanges(),
                                  "GM19239_rep1"=GRanges(),"GM19239_rep2"=GRanges(),
                                  "GM19240_rep1"=GRanges(),"GM19240_rep2"=GRanges())

rep2_peaks_wout_sex_chr_HiChIP=rep2_peaks_with_motif_HiChIP=rep2_peaks_reduced_with_motif_HiChIP=rep2_peaks_raw_HiChIP

#Create GRangesLists that will conclude all the trios data

for (cell_line in cell_lines){
  
  rep2_peaks_raw_HiChIP[[cell_line]]=read_peaks(get_peaks_path(experiment, cell_line))
  
  rep2_peaks_with_motif_HiChIP[[cell_line]]=rep2_peaks_raw_HiChIP[[cell_line]][unique(queryHits(findOverlaps(rep2_peaks_raw_HiChIP[[cell_line]], CTCF.hg38)))]
  
  rep2_peaks_reduced_with_motif_HiChIP[[cell_line]]=IRanges::reduce(rep2_peaks_with_motif_HiChIP[[cell_line]]+p_exp)
  
  rep2_peaks_wout_sex_chr_HiChIP[[cell_line]]=subset(
    rep2_peaks_with_motif_HiChIP[[cell_line]],
    !(rep2_peaks_with_motif_HiChIP[[cell_line]]@seqnames %in% c('chrY', 'chrX', 'chrM')))
}
 
```

# [ChIA-PET single replicates] Comparison between the replicates (raw and peaks with motif)

```{r}
experiment<<-"chiapet"
cell_lines=c("GM19238", "GM19239", "GM19240")

experiment_2<<-experiment_1<<-experiment
for (cell_line in cell_lines){
  # Specify the runs' name
  run_1<<-paste0(cell_line, "_rep1")
  run_2<<-paste0(cell_line, "_rep2")
  
  # Draw the plots
  png(paste0('plots/venn_diagrams/Common_raw_peaks_ChIA-PET_', cell_line, '_replicates.png'), res = 300, width = 2400, height = 2400)
  plot.new()
  draw_2_venn_diagram_from_raw(NROW(rep2_peaks_raw_ChIAPET[[run_1]]), 
                      NROW(rep2_peaks_raw_ChIAPET[[run_2]]), 
                      NROW(subsetByOverlaps(rep2_peaks_raw_ChIAPET[[run_1]], 
                                            rep2_peaks_raw_ChIAPET[[run_2]])))
  dev.off()
}


experiment_2<<-experiment_1<<-experiment
for (cell_line in cell_lines){
  # Specify the runs' name
  run_1<<-paste0(cell_line, "_rep1")
  run_2<<-paste0(cell_line, "_rep2")
  
  # Draw the plots
  png(paste0('plots/venn_diagrams/Common_peaks_with_motif_ChIA-PET_', cell_line, '_replicates.png'), res = 300, width = 2400, height = 2400)
  plot.new()
  draw_2_venn_diagram_from_raw(NROW(rep2_peaks_with_motif_ChIAPET[[run_1]]), 
                      NROW(rep2_peaks_with_motif_ChIAPET[[run_2]]), 
                      NROW(subsetByOverlaps(rep2_peaks_with_motif_ChIAPET[[run_1]], 
                                            rep2_peaks_with_motif_ChIAPET[[run_2]])))
  dev.off()
}
```

# [HiChIP single replicates] Comparison between the replicates

```{r}
experiment<<-"hichip"
cell_lines=c("GM19238", "GM19239", "GM19240")
common_peaks_HiChIP_replicates_raw=common_peaks_HiChIP_replicates_motif=
  GRangesList("GM19238"=GRanges(), "GM19239"=GRanges(), "GM19240"=GRanges())

experiment_2<<-experiment_1<<-experiment

for (cell_line in cell_lines){
  # Specify the runs' name
  run_1<<-paste0(cell_line, "_rep1")
  run_2<<-paste0(cell_line, "_rep2")

  if (NROW(rep2_peaks_raw_HiChIP[[run_1]]) > NROW(rep2_peaks_raw_HiChIP[[run_2]])){
    common_peaks_HiChIP_replicates_raw[[cell_line]] <- subsetByOverlaps(rep2_peaks_raw_HiChIP[[run_2]], 
                                                                    rep2_peaks_raw_HiChIP[[run_1]])
  } else{
  common_peaks_HiChIP_replicates_raw[[cell_line]] <- subsetByOverlaps(rep2_peaks_raw_HiChIP[[run_1]], 
                                                                  rep2_peaks_raw_HiChIP[[run_2]])
  }
  # Draw the plots
  png(paste0('plots/venn_diagrams/Common_raw_peaks_HiChIP_trial_', cell_line, '_replicates.png'), res = 300, width = 2400, height = 2400)
  plot.new()
  draw_2_venn_diagram_from_raw(NROW(rep2_peaks_raw_HiChIP[[run_1]]), 
                      NROW(rep2_peaks_raw_HiChIP[[run_2]]), 
                      NROW(common_peaks_HiChIP_replicates_raw[[cell_line]]))
  dev.off()
}

for (cell_line in cell_lines){
  # Specify the runs' name
  run_1<<-paste0(cell_line, "_rep1")
  run_2<<-paste0(cell_line, "_rep2")

  if (NROW(rep2_peaks_raw_HiChIP[[run_1]]) > NROW(rep2_peaks_raw_HiChIP[[run_2]])){
    common_peaks_HiChIP_replicates_motif[[cell_line]] <- subsetByOverlaps(rep2_peaks_with_motif_HiChIP[[run_2]], 
                                                                    rep2_peaks_with_motif_HiChIP[[run_1]])
  } else{
  common_peaks_HiChIP_replicates_motif[[cell_line]] <- subsetByOverlaps(rep2_peaks_with_motif_HiChIP[[run_1]], 
                                                                  rep2_peaks_with_motif_HiChIP[[run_2]])
  }
  # Draw the plots
  png(paste0('plots/venn_diagrams/Common_peaks_with_motif_HiChIP_trial_', cell_line, '_replicates.png'), res = 300, width = 2400, height = 2400)
  plot.new()
  draw_2_venn_diagram_from_raw(NROW(rep2_peaks_with_motif_HiChIP[[run_1]]), 
                      NROW(rep2_peaks_with_motif_HiChIP[[run_2]]), 
                      NROW(common_peaks_HiChIP_replicates_motif[[cell_line]]))
  dev.off()
}

```
