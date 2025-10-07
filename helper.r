library(opentimsr)
library(data.table)
library(Rcpp)
sourceCpp("one_over_k0_to_ccs.cpp")
sourceCpp('peakalign.cpp')
sourceCpp("one_over_k0_to_ccs.cpp")
sourceCpp("peak_finder.cpp")

getrefpeak <- function(path,accept_Bruker_EULA_and_on_Windows_or_Linux,libpath, batch_size, outref, outcoord, xrange,yrange){
  # download dll/so file and set the column to be collected
  if (accept_Bruker_EULA_and_on_Windows_or_Linux) {
    folder_to_stode_priopriatary_code = libpath
    path_to_bruker_dll = download_bruker_proprietary_code(folder_to_stode_priopriatary_code)
    setup_bruker_so(path_to_bruker_dll)
    all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  }
  # Build the connection with raw data
  D = OpenTIMS(path)
  # extract frame information or define a range to subset the data
  xy <- table2df(D, "MaldiFrameInfo")
  # subset coords idx for region of interests
  idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- dim(xy$MaldiFrameInfo[idx,])[1]
  frame <- xy$MaldiFrameInfo$Frame[idx]
  
  # Process data in batches
  
  total_batches <- ceiling(total_frames / batch_size)
  
  dtx_list <- list()
  for (i in 1:total_batches) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, total_frames)
    idx <- start_idx:end_idx
    
    data <- query(D, frames = frame[idx], columns = all_columns)
    data <- as.data.table(data)
    # round up to make sure binning is working and fit the mass/ion mobility accuracy
    data[, mz := round(mz, 4)]
    data[, inv_ion_mobility := round(inv_ion_mobility, 3)]
    
    dt <- data[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
    dtx_list[[i]] <- dt  
    rm(data)
    gc()
  }
  # reshape the data by mz and ion mobility
  dtx <- rbindlist(dtx_list)
  dtx <- dtx[, .(intensity = sum(intensity)), by = .(mz, inv_ion_mobility)]
  
  # add ccs value for super pixel peaks
  # change to libtimsdata.dll for windows
  lib_path <- paste0(libpath,"libtimsdata.so")
  result <- one_over_k0_to_ccs_parallel(dtx$inv_ion_mobility, rep(1,length(dtx$inv_ion_mobility)), dtx$mz, lib_path)
  dtx[,ccs:=result]
  # save super pixel peaks
  fwrite(dtx,outref)
  # save the location of pixels
  coord <- xy$MaldiFrameInfo[-idx,c('Frame','XIndexPos','YIndexPos')]
  fwrite(coord,outcoord)
}
# This section define the function for peak picking with rcpp support.
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
# function to generate quantative peaks list
getqlist <- function(refpath,lib_path,path, method, zero_proportion_cutoff, coordpath, normpath){
  ref <- fread(refpath)
  setup_bruker_so(lib_path)
  all_columns = c('frame', 'intensity', 'mz', 'inv_ion_mobility')
  D = OpenTIMS(path)
  xy <- table2df(D, "MaldiFrameInfo")
  idx <- which(xy$MaldiFrameInfo$XIndexPos>xrange[1]&xy$MaldiFrameInfo$XIndexPos<xrange[2]&xy$MaldiFrameInfo$YIndexPos>yrange[1]&xy$MaldiFrameInfo$YIndexPos<yrange[2])
  # set a batch size to avoid memory issues by process the data in batches
  total_frames <- dim(xy$MaldiFrameInfo[idx,])[1]
  frame <- xy$MaldiFrameInfo$Frame[idx]
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
    if(method == 'tic'){
      # perform TIC normalization
      resultt[, TIC := sum(intensity, na.rm = TRUE), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / TIC * .N), by = frame]
      resultt[, TIC := NULL]
    }else{
      # perform RMS normalization
      resultt[, rms := sqrt(mean(sum(intensity^2, na.rm = TRUE))), by = frame]
      resultt[, normalized_intensity := as.numeric(intensity / rms), by = frame]
      resultt[, rms := NULL]
    }
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
  columns_to_keep <- names(zero_proportion)[zero_proportion <= zero_proportion_cutoff]
  cn <- columns_to_keep[!is.na(columns_to_keep)]
  result_filtered <- result[, ..cn, with = FALSE]
  result_filtered[,frame:=result$frame]
  # save the data with pixel location information
  coord <- fread(coordpath)
  coord[,location:=paste0(coord$XIndexPos,'_',coord$YIndexPos)]
  result_filtered[,location:=coord$location[match(result_filtered$frame,coord$Frame)]]
  setcolorder(result_filtered, c("location", setdiff(names(result_filtered), "location")))
  result_filtered[, frame := NULL]
  fwrite(result_filtered,normpath)
}
getanno <- function(database,mode,peakpath,annofile){
  lipid <- fread(database)
  # set mode
  if(mode == 'pos'){
    lipid <- lipid[adducts%in%c('[M+H]','[M+Na]','[M+NH4]','[M+H-H2O]')]
  }else{
    lipid <- lipid[adducts%in%c('[M-H]','[M-H+FA]','[M+Cl]')]
  }
  
  # load ref peaks
  ref <- fread(peakpath)
  mz <- sapply(strsplit(colnames(ref)[-1],'\\_'),function(x) round(as.numeric(x[1]),4))
  im <- sapply(strsplit(colnames(ref)[-1],'\\_'),function(x) as.numeric(x[2]))
  # align
  align <- enviGCMS::getalign(mz,lipid$mz,im,lipid$ccs,ppm=20,deltart = 5)
  anno <- cbind.data.frame(mz=mz[align$xid],im=im[align$xid],db_mz=align$mz2,db_im=align$rt2)
  lipidanno <- merge(anno,lipid,by.x = c('db_mz','db_im'),by.y = c('mz','ccs'))
  lipidanno <- lipidanno[!duplicated(lipidanno),]
  # save annotation result as csv file
  fwrite(lipidanno, annofile)
}
