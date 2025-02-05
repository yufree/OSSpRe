library(opentimsr)
library(data.table)
library(Rcpp)

# set the path to the .d folder
path = 'PATH_TO_RAW.d'
# accept Bruker's lisence
accept_Bruker_EULA_and_on_Windows_or_Linux = TRUE
# download dll/so file and set the column to be collected
if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
  folder_to_stode_priopriatary_code = "~"
  path_to_bruker_dll = download_bruker_proprietary_code(folder_to_stode_priopriatary_code)
  setup_bruker_so(path_to_bruker_dll)
  all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
}
# Build the connection with raw data
D = OpenTIMS(path)
# extract frame information or define a range to subset the data
xy <- table2df(D, "MaldiFrameInfo")
# subset coords idx for region of interests
idx <- which(xy$MaldiFrameInfo$XIndexPos>2900&xy$MaldiFrameInfo$XIndexPos<3000)
# set a batch size to avoid memory issues by process the data in batches
total_frames <- dim(xy$MaldiFrameInfo[-idx,])[1]
frame <- xy$MaldiFrameInfo$Frame[-idx]
batch_size <- 100000
total_batches <- ceiling(total_frames / batch_size)
# save the data as h5 files with memory release
dtx_list <- list()
for (i in 1:total_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, total_frames)
  idx <- frame[start_idx:end_idx]

  data <- query(D, frames = idx, columns = all_columns)
  data <- as.data.table(data)
  data[, mz := round(mz, 4)]
  data[, inv_ion_mobility := round(inv_ion_mobility, 3)]

  dt <- data[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
  dtx_list[[i]] <- dt
  rm(data)
  gc()
}

dtx <- rbindlist(dtx_list)
dtx <- dtx[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]

# add ccs value for super pixel peaks
sourceCpp("one_over_k0_to_ccs.cpp")
lib_path <- "libtimsdata.so"
result <- one_over_k0_to_ccs_parallel(dtx$inv_ion_mobility, rep(1,length(dtx$inv_ion_mobility)), dtx$mz, lib_path)
dtx[,ccs:=result]
# save the location of pixels
coord <- xy$MaldiFrameInfo[-idx,c('Frame','XIndexPos','YIndexPos')]

sourceCpp("peak_finder.cpp")
find_2d_peaks <- function(mz, ccs, intensity, 
                          mz_ppm = 20, 
                          ccs_tolerance = 0.03,
                          snr = 3.0,
                          mz_bins = 1000,
                          ccs_bins = 50) {
  
  if (!is.numeric(mz) || !is.numeric(ccs) || !is.numeric(intensity)) {
    stop("All inputs must be numeric vectors")
  }
  if (length(mz) != length(ccs) || length(mz) != length(intensity)) {
    stop("All input vectors must have the same length")
  }
  
  find_2d_peaks_spatial_parallel_openmp(mz, ccs, intensity, 
                                        mz_ppm, ccs_tolerance, snr,
                                        mz_bins, ccs_bins)
}

mz_ppm <- 20
ccs_tolerance <- 0.03
snr <- 3
mz_bins <- 1000
ccs_bins <- 50
result <- find_2d_peaks(dtx$mz,dtx$ccs,dtx$intensity,mz_ppm,ccs_tolerance,snr,mz_bins,ccs_bins)
fwrite(result,'ref2d.csv')
