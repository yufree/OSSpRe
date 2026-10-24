setwd("/content/OSSpRe")
if (file.exists("renv/activate.R")) source("renv/activate.R")
if (file.exists("/content/config.R")) source("/content/config.R")
suppressPackageStartupMessages(source("helper.r"))

## ---- data summary ----
cat("Summarising ticmzccs.csv...\n")
getsummary("ticmzccs.csv")

## ---- peak statistics: two ggplots -> two numbered files ----
cat("Plotting peak statistics...\n")
png("peakstats_%d.png", width = 900, height = 600, res = 110)
plot_peak_stats("ticmzccs.csv")
dev.off()

## ---- ion images for the 7 most abundant reference ions ----
cat("Extracting top 7 abundant ions for imaging...\n")
dt <- fread("ticmzccs.csv")
mz_col <- sapply(strsplit(colnames(dt)[-1], '_'), function(x) round(as.numeric(x[1]), 4))
abund  <- colSums(dt[, -1, with = FALSE])
mzdemo <- mz_col[order(-abund)][1:7]
rm(dt, abund)
gc(verbose = FALSE)
save_ion_images("ticmzccs.csv", mzdemo)
cat("ion images:", paste0(mzdemo, ".png"), "\n")

## ---- PCA segmentation ----
cat("Performing PCA segmentation...\n")
png("pca_segmentation.png", width = 800, height = 800, res = 110)
perform_pca_segmentation("ticmzccs.csv", "segmentation.csv")
dev.off()

## ---- UMAP segmentation (optional: needs python umap-learn) ----
if (identical(Sys.getenv("RUN_UMAP"), "TRUE")) {
  tryCatch({
    Sys.setenv(RETICULATE_PYTHON = Sys.which("python3"))
    png("umap_segmentation.png", width = 800, height = 800, res = 110)
    perform_umap_segmentation("ticmzccs.csv", "segmentation_umap.csv")
    dev.off()
  }, error = function(e) message("UMAP segmentation skipped: ", conditionMessage(e)))
} else {
  cat("UMAP segmentation skipped (optional; downstream steps use PCA segmentation.csv)\n")
}

## ---- ion clustering across the whole sample ----
cat("Clustering ions across whole sample (generating ioncluster.csv)...\n")
hclust_cut <- if (exists("HCLUST_CUTOFF")) HCLUST_CUTOFF else 0.6
min_cl_sz  <- if (exists("MIN_CLUSTER_SIZE")) MIN_CLUSTER_SIZE else 10
perform_cluster_ions(peak_file = "ticmzccs.csv", output_file = "ioncluster.csv",
                     hclust_cutoff = hclust_cut,
                     min_cluster_size = min_cl_sz,
                     folder_name = "fullsample")
cat("ioncluster.csv written successfully.\n")

## ---- ROI-specific clustering and reactomics ----
seg      <- fread("segmentation.csv")
cls      <- as.integer(names(sort(table(seg$seg), decreasing = TRUE))[1])
locindex <- seg$seg == cls
cat("ROI segmentation class:", cls, "| pixels:", sum(locindex), "\n")

cat("Performing ROI ion clustering...\n")
perform_cluster_ions(peak_file = "ticmzccs.csv", output_file = "roi_ion_cluster.csv",
                     locindex = locindex, hclust_cutoff = 0.8,
                     min_cluster_size = min_cl_sz, folder_name = "roi")

cat("Performing whole-sample and ROI reactomics analysis...\n")
perform_reactomics_analysis("ticmzccs.csv", "reactomics.csv")
perform_reactomics_analysis("ticmzccs.csv", "roireactomics.csv", locindex = locindex)

cat("Saving PMD images...\n")
save_pmd_images("reactomics.csv", pmd_values = 2.016, "pmd2.016.png")
cat("Step 5 exploratory analysis complete! ioncluster.csv and reactomics.csv are ready.\n")
