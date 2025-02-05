library(data.table)
library(Rcpp)
# load the reference peaks
ref <- fread('ref2d.csv')
sourceCpp('peakalign.cpp')
sourceCpp("one_over_k0_to_ccs.cpp")
lib_path <- "libtimsdata.so"
# load the raw data
library(opentimsr)
path = 'PATH_TO_RAW.d'
setup_bruker_so("libtimsdata.so")
all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
D = OpenTIMS(path)
xy <- table2df(D, "MaldiFrameInfo")
idx <- which(xy$MaldiFrameInfo$XIndexPos>2900&xy$MaldiFrameInfo$XIndexPos<3000)
# set a batch size to avoid memory issues by process the data in batches
total_frames <- dim(xy$MaldiFrameInfo[-idx,])[1]
frame <- xy$MaldiFrameInfo$Frame[-idx]
coord <- xy$MaldiFrameInfo[-idx,c('Frame','XIndexPos','YIndexPos')]
# this setting will make sure the data.frame's row number is within the limits
batch <- floor(2e+09/nrow(ref))
batchs <- ceiling(total_frames/batch)
result <- list()
for(i in 1:batchs){
  upper <- min(total_frames,i*batch)
  idx<- c(((i-1)*batch+1):upper)
  data <- query(D, frames = frame[idx], columns = all_columns)
  data <- as.data.table(data)
  data[,ccs:=one_over_k0_to_ccs_parallel(data$inv_ion_mobility, rep(1,length(data$inv_ion_mobility)), data$mz, lib_path)]
  resultt <- findpeakalign(data$mz,data$ccs,data$intensity,data$frame,ref$mz,ref$ccs, ppm = 20, ccs_shift = 0.03)
  resultt <- as.data.table(resultt)
  # perform TIC normalization
  resultt[, TIC := sum(intensity, na.rm = TRUE), by = frame]
  resultt[, normalized_intensity := as.numeric(intensity / TIC * .N), by = frame]
  resultt[, TIC := NULL]
  # perform RMS normalization
  # resultt[, rms := sqrt(mean(sum(intensity^2, na.rm = TRUE))), by = frame]
  # resultt[, normalized_intensity := as.numeric(intensity / rms), by = frame]
  # resultt[, rms := NULL]
  result[[i]] <- dcast(
    resultt, 
    frame ~ mz_ccs, 
    value.var = "normalized_intensity", 
    fun.aggregate = sum,
    fill = 0)
}
result <- rbindlist(result,fill=T)
for (col in names(result)) {
  set(result, i = which(is.na(result[[col]])), j = col, value = 0)
}
# remove the peaks with zero_proportion larger than 95%
zero_proportion <- result[, lapply(.SD, function(x) mean(x == 0)), .SDcols = 2:ncol(result)]
columns_to_keep <- names(zero_proportion)[zero_proportion <= 0.95]
cn <- columns_to_keep[!is.na(columns_to_keep)]
result_filtered <- result[, ..cn, with = FALSE]
result_filtered[,frame:=result$frame]
# save the data with pixel location information
coord[,location:=paste0(coord$XIndexPos,'_',coord$YIndexPos)]
result_filtered[,location:=coord$location[match(result_filtered$frame,coord$Frame)]]
setcolorder(result_filtered, c("location", setdiff(names(result_filtered), "location")))
result_filtered[, frame := NULL]
fwrite(result_filtered,'ticmzccs.csv')

proc.time() - ptm