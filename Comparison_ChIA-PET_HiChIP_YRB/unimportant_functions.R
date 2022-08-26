# Tu będą się znajdywać funkcje, które raczej już nie są potrzebne

# Raczej ta funcka nie jest już potrzebna, została zastąpiona przez true_overlaps
find_main_overlaps_between_3=function(peaks1, peaks2, peaks3){
  nrow_1=NROW(peaks1)
  nrow_2=NROW(peaks2)
  nrow_3=NROW(peaks3)
  
  if (nrow_1<nrow_2){common_12=subsetByOverlaps(peaks1, peaks2)} 
  else {common_12=subsetByOverlaps(peaks2, peaks1)}
  if (nrow_2<nrow_3){common_23=subsetByOverlaps(peaks2, peaks3)}
  else {common_23=subsetByOverlaps(peaks3, peaks2)}
  if (nrow_1<nrow_3){common_13=subsetByOverlaps(peaks1, peaks3)}
  else {common_13=subsetByOverlaps(peaks3, peaks1)}
  
  if (nrow_3>nrow_1 & nrow_3>nrow_2){common_peaks=subsetByOverlaps(common_12, peaks3)}
  else if (nrow_2>nrow_1 & nrow_2>nrow_3){common_peaks=subsetByOverlaps(common_13, peaks2)}
  else if (nrow_1>nrow_2 & nrow_1>nrow_3){common_peaks=subsetByOverlaps(common_23, peaks1)}
  
  common_overlaps=GRangesList("com12"=common_12,"com23"=common_23, "com13" = common_13, "com123"= common_peaks)
  return(common_overlaps)
}