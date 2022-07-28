#!/bin/bash

chr='"chr1"'
petN=3
technology="chiapet"

cell_line="GM19238_merged_replicates"
path='/Users/aleksander.ostruk/Documents/asia_master/'

cat $path/$technology/$cell_line/$cell_line.e500.clusters.cis.BE3.peak_annot.E2 |

awk '$1 == '"$chr" |
awk '$4 == '"$chr" |
awk '($8 = ($2 + $3) /2)' |
awk '($9 = ($5 + $6) /2)' |
awk '($10 = int($9 - $8))' |
awk '{print $10, $7}' |
awk '$2>='"$petN" |

# specify the intervals
#sort -g |
#minbp=7000
#times=2

#awk '{if ($1 >= "7000" && $1 <= "14000"}' | head


# create a file
sort -g > $path/$technology/$cell_line/length_analysis/len_freq_"$cell_line"_"$(echo $chr | sed 's/"//g')"_pet$petN.csv

