library(opentimsr)
library(data.table)
library(Rcpp)
library(enviGCMS)
library(irlba)
library(ClusterR)
library(reticulate)
library(umap)
library(dbscan)
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
getqlist <- function(refpath,lib_path,path, method, zero_proportion_cutoff, coordpath, normpath, xrange, yrange){
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

plot_peak_stats <- function(peak_file) {
  library(data.table)
  library(ggplot2)

  dt <- fread(peak_file, header = TRUE)
  dt[, np := rowSums(.SD != 0)]
  
  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))

  df <- cbind.data.frame(x = x, y = y, np = dt$np)

  p1 <- ggplot(df, aes(np)) +
    ggtitle('Peak number for each pixel') + xlab('Peak number') +
    geom_histogram(binwidth = 1) + theme_bw()

  p2 <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = np), size = 0.001) +
    scale_color_gradient(low = "white", high = "red") + theme_void()

  print(p1)
  print(p2)
}

save_ion_images <- function(peak_file, mz_values) {
  library(data.table)

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))

  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))

  dt_filtered <- as.data.frame(dt_values[, .SD, .SDcols = mz %in% mz_values])

  width <- max(x) - min(x) + 1
  height <- max(y) - min(y) + 1

  for (i in 1:ncol(dt_filtered)) {
    dfx <- dt_filtered[, i]
    norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx))
    png(paste0(colnames(dt_filtered)[i], '.png'), width = width, height = height)
    plot.new()
    par(mar = c(0, 0, 0, 0))
    plot.window(xlim = c(0, width), ylim = c(0, height), xaxs = "i", yaxs = "i", asp = NA)
    points(x - min(x) + 1, y - min(y) + 1, pch = 16, col = grDevices::gray(1 - norm), cex = 0.3)
    dev.off()
  }
}

perform_pca_segmentation <- function(peak_file, output_file, n_components = 20, n_clusters = 5) {

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))

  mat <- t(dt_values)
  mat_centered <- scale(t(mat), center = TRUE, scale = FALSE)
  svd_result <- irlba(t(mat_centered), nv = n_components)
  pca_scores <- t(mat_centered) %*% svd_result$v

  km <- KMeans_arma(as.matrix(pca_scores), clusters = n_clusters, n_iter = 10, seed_mode = "random_subset", verbose = TRUE, CENTROIDS = NULL)
  pr <- predict_KMeans(as.matrix(pca_scores), km)

  plot(x, y, col = pr, cex = 0.1, pch = 19)
  legend('topright', legend = unique(pr), col = unique(pr), pch = 19, cex = 1)

  seg <- cbind.data.frame(x = x, y = y, pca = pr)
  fwrite(seg, output_file)
}

perform_umap_segmentation <- function(peak_file, output_file, n_threads = 50, eps = 0.2, minPts = 20) {

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  x <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[1]))
  y <- sapply(strsplit(dt$location, '_'), function(x) as.numeric(x[2]))

  mat <- t(dt_values)
  Sys.setenv(NUMBA_NUM_THREADS = n_threads)
  py_run_string(paste0("import os; print(os.environ['NUMBA_NUM_THREADS'])"))

  viz <- umap::umap(t(mat), method = 'umap-learn', metric = 'cosine')
  dbscan_result <- dbscan(viz$layout, eps = eps, minPts = minPts)
  
  plot(x = x, y = y, col = dbscan_result$cluster+1,xlab='',ylab = '',main='',xaxt = "n", yaxt = "n",bty = "n",cex=0.1)
  legend('bottomright',legend = unique(dbscan_result$cluster+1),col = unique(dbscan_result$cluster+1), pch=19,bty = "n")

  seg <- cbind.data.frame(x = x, y = y, umap = dbscan_result$cluster+1)
  fwrite(seg, output_file)
}

cluster_ions <- function(peak_file, output_file, hclust_cutoff = 0.6, min_cluster_size = 10) {
  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))

  Matrix <- t(dt_values)
  row_norms <- sqrt(rowSums(Matrix^2))
  sim <- sweep(Matrix, 1, row_norms, "/")
  sim <- sim %*% t(sim)
  D_sim <- as.dist(1 - sim)

  t <- hclust(D_sim)
  s <- cutree(t, h = hclust_cutoff)

  name <- as.numeric(names(table(s)[table(s) > min_cluster_size]))
  
  Matrix <- t(dt_values)
  split_matrices <- lapply(name, function(category) {
    rows <- which(s == as.numeric(category))
    subset_matrix <- Matrix[rows, , drop = FALSE]
    return(subset_matrix)
  })
  
  summed_matrices <- lapply(split_matrices, function(subset_matrix) {
    xxx <- apply(subset_matrix, 1, scale)
    x <- rowSums(xxx) / ncol(xxx)
    return(x)
  })
  
  result_matrix <- do.call(cbind, summed_matrices)
  
  clpan <- cbind.data.frame(x = x, y = y, result_matrix)
  
  dir.create('cluster')
  for (i in c(1:length(name))) {
    dfx <- result_matrix[, i]
    norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx))
    color_palette <- colorRamp(c("yellow", "red"))
    color_sequence <- rgb(color_palette(norm) / 255, alpha = 1)
    xlim <- max(x) - min(x) + 1
    ylim <- max(y) - min(y) + 1
    
    png(paste0('cluster/cluster', name[i], '.png'), width = xlim, height = ylim)
    plot.new()
    par(mar = c(0, 0, 0, 0))
    plot.window(xlim = c(0, xlim), ylim = c(0, ylim), xaxs = "i", yaxs = "i", asp = NA)
    points(x - min(x) + 1, y - min(y) + 1, pch = 16, col = color_sequence, cex = 0.3)
    dev.off()
  }
  ioncluster <- cbind.data.frame(mz, im, class = s)
  fwrite(ioncluster, output_file)
}

cluster_roi_ions <- function(peak_file, segmentation_file, output_file,roi_cluster = 2, hclust_cutoff = 0.6, min_cluster_size = 10) {
  library(data.table)

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))

  seg <- fread(segmentation_file, header = TRUE)

  subdt <- dt_values[seg$umap == roi_cluster, ]
  x <- seg$x
  y <- seg$y

  Matrix <- t(subdt)

  n <- ncol(Matrix)
  Matrix_scaled <- scale(Matrix)
  sim <- (Matrix_scaled %*% t(Matrix_scaled)) / (n - 1)
  D_sim <- as.dist(1 - sim)

  t <- hclust(D_sim)
  s <- cutree(t, h = hclust_cutoff)

  name <- as.numeric(names(table(s)[table(s) > min_cluster_size]))

  Matrix <- t(dt_values)
  split_matrices <- lapply(name, function(category) {
    rows <- which(s == as.numeric(category))
    subset_matrix <- Matrix[rows, , drop = FALSE]
    return(subset_matrix)
  })

  summed_matrices <- lapply(split_matrices, function(subset_matrix) {
    xxx <- apply(subset_matrix, 1, scale)
    x <- rowSums(xxx) / ncol(xxx)
    return(x)
  })

  result_matrix <- do.call(cbind, summed_matrices)

  clpan <- cbind.data.frame(x = x, y = y, result_matrix)

  dir.create('cluster')
  for (i in c(1:length(name))) {
    dfx <- result_matrix[, i]
    norm <- (dfx - min(dfx)) / (max(dfx) - min(dfx))
    color_palette <- colorRamp(c("yellow", "red"))
    color_sequence <- rgb(color_palette(norm) / 255, alpha = 1)
    xlim <- max(x) - min(x) + 1
    ylim <- max(y) - min(y) + 1

    png(paste0('cluster/cluster', name[i], '.png'), width = xlim, height = ylim)
    plot.new()
    par(mar = c(0, 0, 0, 0))
    plot.window(xlim = c(0, xlim), ylim = c(0, ylim), xaxs = "i", yaxs = "i", asp = NA)
    points(x - min(x) + 1, y - min(y) + 1, pch = 16, col = color_sequence, cex = 0.3)
    dev.off()
  }
  ioncluster <- cbind.data.frame(mz, im, class = s)
  fwrite(ioncluster, output_file)
}

perform_reactomics_analysis <- function(peak_file, segmentation_file, islet_file, output_file) {
  library(pmd)
  library(data.table)

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))

  seg <- fread(segmentation_file, header = TRUE)
  isletdf <- fread(islet_file, header = TRUE)

  df <- getpmddf(mz, pmd = c(14.015, 14.016, 15.995, 15.996, 2.015, 2.016, 18.010, 18.011), digits = 3, group = ifelse(colnames(dt)[-1] %in% isletdf$`islet`, 1, 0))
  xy <- data.frame(x = seg$x, y = seg$y)
  xy$factor <- seg$umap

  df$diff3 <- round(df$diff, 3)
  df3 <- df[df$diff3 == 14.015 | df$diff3 == 14.016, ]
  df4 <- df[df$diff3 == 15.995 | df$diff3 == 15.996, ]
  df5 <- df[df$diff3 == 2.015 | df$diff3 == 2.016, ]
  df6 <- df[df$diff3 == 18.010 | df$diff3 == 18.011, ]

  dfall <- apply(dt_values, 1, function(x) x / max(x))
  dfall[is.na(dfall)] <- 0
  dfall <- t(dfall)

  dfms1 <- dfall[, mz %in% df3$ms1]
  dfms2 <- dfall[, mz %in% df3$ms2]
  xy$pmd14h <- apply(dfms1, 1, sum)
  xy$pmd14l <- apply(dfms2, 1, sum)

  dfms1 <- dfall[, mz %in% df4$ms1]
  dfms2 <- dfall[, mz %in% df4$ms2]
  xy$pmd16h <- apply(dfms1, 1, sum)
  xy$pmd16l <- apply(dfms2, 1, sum)

  dfms1 <- dfall[, mz %in% df5$ms1]
  dfms2 <- dfall[, mz %in% df5$ms2]
  xy$pmd2h <- apply(dfms1, 1, sum)
  xy$pmd2l <- apply(dfms2, 1, sum)

  dfms1 <- dfall[, mz %in% df6$ms1]
  dfms2 <- dfall[, mz %in% df6$ms2]
  xy$pmd18h <- apply(dfms1, 1, sum)
  xy$pmd18l <- apply(dfms2, 1, sum)

  xy$pmd2 <- xy$pmd2h + xy$pmd2l
  xy$pmd14 <- xy$pmd14h + xy$pmd14l
  xy$pmd16 <- xy$pmd16h + xy$pmd16l
  xy$pmd18 <- xy$pmd18h + xy$pmd18l

  fwrite(xy, output_file)
}

generate_molecular_network <- function(peak_file, ion_cluster_file, annotation_file, output_file) {
  library(pmd)
  library(data.table)
  library(igraph)
  library(ggraph)
  library(tidygraph)

  dt <- fread(peak_file, header = TRUE)
  dt_values <- dt[, -1, with = FALSE]

  mz <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
  im <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) as.numeric(x[2]))

  ioncluster <- fread(ion_cluster_file)
  anno <- fread(annotation_file)

  data("keggrall")
  hfpmd <- as.numeric(names(table(keggrall$pmd)[table(keggrall$pmd) > mean(table(keggrall$pmd))]))
  df <- getpmddf(mz, group = ioncluster$class, pmd = hfpmd, digits = 3)

  df$anno1 <- anno$name[match(df$ms1, anno$mz)]
  df$anno2 <- anno$name[match(df$ms2, anno$mz)]

  dfanno <- df[complete.cases(df), ]

  net <- graph_from_data_frame(dfanno)
  graph <- as_tbl_graph(net)
  p <- ggraph(graph, layout = 'fr') +
    geom_edge_link(aes(color = net)) + theme_void()

  print(p)

  fwrite(dfanno, output_file)
}