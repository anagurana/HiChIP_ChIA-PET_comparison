library(MAnorm2)
library(caret)

### Function that generates MAnorm single plot and scales it properly

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
  MAplot(param[[2]], param[[1]], param[[4]], param[[3]], ylim = set_ylim,
         main=main, cex=1.5, cex.main=2.3, 
         line=0.2,
         cex.axis=2, cex.lab=2,
         args.legend = list(
           x="topright",
           inset = c(0, 0),
           legend=c("common", paste0(cells[1]," specific"), paste0(cells[2]," specific"), "others"),
           bty = "o", xpd = TRUE,
           y.intersp=0.8,
           x.intersp=1,
           cex=1.5))
  abline(h = 0, lwd = 2, lty = 5)
  
}

### Function that combine data from different tables and create 
### single data.frame with all important info: chrom start, end, mean values, pval and ect.
### !important - p_value is being defined in the script as global variable

create_table_common_peaks_with_info=function(norm, differential, conds, cell_lines){
  common=data.frame(chrom=norm$chrom,
                    start=norm$start,
                    end=norm$end,
                    conds[[cell_lines[1]]]$occupancy,
                    conds[[cell_lines[2]]]$occupancy,
                    differential[[paste0(cell_lines[1],"_",cell_lines[2])]],
                    conds[[cell_lines[1]]]$norm.signal,
                    conds[[cell_lines[2]]]$norm.signal)
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

### Function that creates short names for the columns in shorter format like "GM38_GM39"

create_short_name_for_col=function(trios){
  names=c(
    paste0(paste0(str_sub(trios[1],1,2),str_sub(trios[1],start=6)),"_",paste0(str_sub(trios[2],1,2),str_sub(trios[2],start=6))),
    paste0(paste0(str_sub(trios[1],1,2),str_sub(trios[1],start=6)),"_",paste0(str_sub(trios[3],1,2),str_sub(trios[3],start=6))),
    paste0(paste0(str_sub(trios[2],1,2),str_sub(trios[2],start=6)),"_",paste0(str_sub(trios[3],1,2),str_sub(trios[3],start=6))))
  return(names)
}

### Function that creates separate tables with mean value cut-off defined as the parameter

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

### Function that creates separate tables with peaks that are present/not present in each cell line
### As the result there will be tables like cell1T_cell2T, cell1T_cell2F,cell1F_cell2T,cell1F_cell2F

create_fp_tables_and_save_graphs=function(common,cell_lines, experiment, prefix, padj_or_pval="isPadj"){
  # Devide picks into different groups
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
  
  png(paste0("plots/mean_values/",experiment, "_MAnorm_2rep_", cell1,"_",cell2,"_diff_test_",prefix,".png"), res = 300, width = 3500, height = 2200)
  plot.new()
  p1=ggplot(tf_table[[1]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("Common peaks in 2 cell lines")+labs(x=cell1.mean, y=cell2.mean)
  p2=ggplot(tf_table[[2]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("Uncommon peaks in 2 cell lines")+labs(x=cell1.mean, y=cell2.mean)
  p3=ggplot(tf_table[[3]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle('No peaks\nin both cell lines')+theme(legend.position = "none")+labs(x=cell1.mean, y=cell2.mean)
  p4=ggplot(tf_table[[4]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle(paste0('Peaks only\nin ', cell_lines[1]))+theme(legend.position = "none")+labs(x=cell1.mean, y=cell2.mean)
  p5=ggplot(tf_table[[5]], aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle(paste0('Peaks only\nin ',cell_lines[2]))+labs(x=cell1.mean, y=cell2.mean)
  p6=ggplot(common, aes(x=get(cell1.mean), y=get(cell2.mean)))+geom_point(aes(col=get(padj_or_pval)))+ggtitle("All peaks")+geom_smooth(method="lm", formula=y~x)+labs(x=cell1.mean, y=cell2.mean)
  print(p6+(p1 | p2) / (p3 | p4 | p5) + plot_layout(widths=c(1, 2)))
  dev.off()
  
  return(tf_table)
}

