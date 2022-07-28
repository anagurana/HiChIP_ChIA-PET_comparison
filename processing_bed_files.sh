
samtools view GM1923_merged_replicates.for.BROWSER.bam | sam2bed -o GM19238_merged_replicates.for.BROWSER.bed

profile_bins --peaks=/home2/sfglab/agurova/data/chiapet/GM19238/peaks_ChIAPET_GM19238.bed,/home2/sfglab/agurova/data/chiapet/GM19239/peaks_ChIAPET_GM19239.bed,/home2/sfglab/agurova/data/chiapet/GM19240/peaks_ChIAPET_GM19240.bed  \
--reads=/home2/sfglab/agurova/data/chiapet/GM19238/GM19238_merged_replicates.for.BROWSER.bed,/home2/sfglab/agurova/data/chiapet/GM19239/GM19239_merged_replicates.for.BROWSER.bed,/home2/sfglab/agurova/data/chiapet/GM19240/GM19240_merged_replicates.for.BROWSER.bed \
--labs=GM19238,GM19239,GM19240 \
-n reference_bins \
--min-peak-gap=1000 \
#--paired \
#--keep-dup=1 \
--bins=/home2/sfglab/agurova/data/chiapet/bins_ChIAPET.bed


#--min-peak-gap=1000 \# (default 150) \