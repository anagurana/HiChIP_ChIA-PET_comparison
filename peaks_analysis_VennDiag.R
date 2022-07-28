## Libraries installations, commented due to overload in the console

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install()
# BiocManager::install("ChIPseeker")
# BiocManager::install("GenomicRanges")
# install.packages('VennDiagram')
# install.packages("gridExtra")       # Install gridExtra package
# suppressMessages(library(ChIPseeker))
# suppressMessages(library(GenomicRanges))
# suppressMessages(library(VennDiagram))
# suppressMessages(library(gridExtra))

source("/Users/aleksander.ostruk/Documents/asia_master/functions_ChIA-PET_HiChIP_YRB.R")

## Definition of the color for the each technology and cell lines
color=function(word){
  if (word == "GM19238") {color="deeppink"}
  # else if (word == "GM19240_rep1") {color="greenyellow"}
  # else if (word == "GM19240_rep2") {color="green4"}
  else if (word == "GM19240_rep1") {color="cornflowerblue"}
  else if (word == "GM19240_rep2") {color="darkblue"}
  else if (word == "GM19239") {color="deepskyblue"}
  else if (word == "GM19240") {color="yellow"}
  else if (word == "chiapet") {color="blue"}
  else if (word == "hichip") {color="lawngreen"}
  else if (word == "chipseq") {color="red"}
  return(color)
}

## Converting the names of the technologies to the proper formats
prop_name=function(technology){
  if (technology == "chiapet") {name="ChIA-PET"}
  else if (technology == "hichip") {name="HiChIP"}
  else if (technology == "chipseq") {name="ChIP-Seq"}
  return(name)
}

draw_2_venn_diagram <- function(peaks1, peaks2) {

  common_peaks=find_main_overlaps(peaks1, peaks2)
  
  area1 = NROW(peaks1)
  area2 = NROW(peaks2)
  cross.area = NROW(common_peaks)
  
  graphics.off()
  
  # FOR hichip 38 and 40
  
  if (experiment_1==experiment_2){
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
                                              scaled = T,
                                              fill=c(color(cell1), color(cell2)),
                                              category=c(cell1, cell2),
                                              print.mode=c('raw','percent'),
                                              reversed=TRUE, cat.dist=-0.45, cat.cex=3.5, cat.pos =  c(-30, 30),
                                              cex=4, margin=0.05,
                                              ext.percent=.2)
    grid.arrange(gTree(children = venndiag),
                 top = textGrob( paste0("Peaks comparison between the replicates for ", prop_name(experiment_1)), gp=gpar(cex=4)))
    } else {
    venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area, 
                                              scaled = T,
                                              fill=c(color(experiment_1), color(experiment_2)),
                                              category=c(prop_name(experiment_1), prop_name(experiment_2)),
                                              print.mode=c('raw'), 
                                              cat.cex=1.3, cat.dist=0.05, cat.pos=c(-120, 120), cex=1.5, 
                                              ext.percent=0.3, margin=0.05)
    venndiag.ga=grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between technologies for ", cell1))
    }
  if (save == TRUE){
    print("save")
  }
    #for hichip 39
    # if (exp1==exp2){
    #   venndiag=VennDiagram::draw.pairwise.venn( area1, area2, cross.area,
    #                                             scaled = T,
    #                                             fill=c(color(cell1), color(cell2)),
    #                                             category=c(cell1, cell2),
    #                                             print.mode=c('raw','percent'),
    #                                             inverted=FALSE, cat.dist=0.05, cat.cex=3.5, cat.pos =  c(220, -220),
    #                                             cex=4, margin=0.05,
    #                                             ext.percent=.2)
    #   grid.arrange(gTree(children = venndiag),
    #                top = textGrob( paste0("Peaks comparison between the replicates for ", prop_name(exp1)), gp=gpar(cex=4)))
    
}

draw_3_venn_diagram<- function(peaks1, peaks2, peaks3){
  
  common_peaks=find_main_overlaps(peaks1, peaks2, peaks3)
  
  #Define the required areas for Venn Diagram
  area1=NROW(peaks1)
  area2=NROW(peaks2)
  area3=NROW(peaks3)
  
  n12=NROW(common_peaks$com12)
  n23=NROW(common_peaks$com23)
  n13=NROW(common_peaks$com13)
  n123=NROW(common_peaks$com123)
  
  #Draw plots
  graphics.off()
  if (experiment_1 == experiment_2 & experiment_1 == experiment_3){
    venndiag = draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, 
                                fill=c(color(cell1), color(cell2), color(cell3)),
                                category=c(cell1, cell2, cell3),
                                print.mode=c('raw','percent'))
    grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between cell lines for ", prop_name(experiment_1)))
  } else {
    venndiag = draw.triple.venn(area1, area2, area3, n12, n23, n13, n123, 
                                fill=c(color(experiment_1), color(experiment_2), color(experiment_3)),
                                category=c(prop_name(experiment_1), prop_name(experiment_2), prop_name(experiment_3)),
                                print.mode=c('raw', 'percent'))
    venndiag.ga=grid.arrange(gTree(children = venndiag), top = paste0("Peaks comparison between technologies for ", cell1))
  }
  return(venndiag.ga)
}

# Path is the place, which contains directory(s) with experiment(s) type,
# which is specified with "experiment" parameter
# Inside of the "experiment" folder, cell lines folders are expected,
# which could be specifies with "cell_line" parameter
abs_path<<-"/Users/aleksander.ostruk/Documents/asia_master"

## Data for the Venn Diagrams, if you want to use only 2 Venn diagram,
## please specify only 1 and 2 parameters

"
IMPORTANT - the pipeline was build in the way, that peaks in the 3 file are
expected to have the highest number of peaks, adn rest in order 3 > 1 > 2
"

experiment_1<<-"chiapet"
experiment_<<-"hichip"

experiment_3<<-"chipseq"

cell_1<<-"GM19238"
cell_2<<-"GM19238"

cell_3<<-"GM19238"

## If you have your directories sorted like was descibed in the introduction,
## please use the function get_peaks_path, otherwise please specify the pathes manually
peak_file_1=get_peaks_path(experiment_1, cell_1)
peak_file_2=get_peaks_path(experiment_2, cell_2)
peak_file_3=get_peaks_path(experiment_3, cell_3)

# Filter you peaks
peaks_1=read_and_filter_peaks(peak_file_1)
peaks_2=read_and_filter_peaks(peak_file_2)
peaks_3=read_and_filter_peaks(peak_file_3)

venn_diagram=draw_3_venn_diagram(peaks_1, peaks_2, peaks_3)
# venn_diagram=draw_2_venn_diagram(peaks_1, peaks_2, save=FALSE)








