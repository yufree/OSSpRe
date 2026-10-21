# OSSpRe

This document is the reproducible code for Open Source Spatial Reactomics (OSSpRe) workflow. This workflow is designed to process mass spectrometry imaging (MSI) data from raw data to exploratory data analysis.

# Install the packages needed for this workflow

We use renv package to restore all the R packages needed for this workflow.

``` r
install.packages('renv')
renv::restore()
```
# Tutorial

Part 1: Introduction: https://youtu.be/j26o3j2B1FI

Part 2: Installation: https://youtu.be/l_3MuLlQbbc

Part 3: Rmarkdown work-through: https://youtu.be/e_BcrRFrNCA

Part 4: Shiny application: https://youtu.be/W7lre03bUuc

# Files

- `workflow.Rmd`: The reproducible R Markdown document for the OSSpRe workflow. It contains all the steps from raw data processing to exploratory data analysis.
- `helper.r`: An R script containing helper functions used in the workflow.
- `lipidall.csv`: A lipid m/z-CCS annotation database from Metlin and Lipid Atlas.
- `metaall.csv`: A metabolomics m/z-CCS annotation database from Metlin.
- `ccs_to_one_over_k0.cpp`: A C++ script to convert CCS values to 1/k0.
- `one_over_k0_to_ccs.cpp`: A C++ script to convert 1/k0 values to CCS.
- `peak_finder.cpp`: A C++ script to find local maximum peaks in the m/z-CCS 2D plane.
- `peakalign.cpp`: A C++ script to extract peaks with a reference peak list across pixels.

# Usage

1.  **Set the path to your raw data:** In the `workflow.Rmd` file, change the `path` variable to the path of your `.d` folder.
2.  **Accept the Bruker EULA:** Set `accept_Bruker_EULA_and_on_Windows_or_Linux` to `TRUE` to accept the Bruker End User License Agreement.
3.  **Run the workflow:** The `workflow.Rmd` chunks are set to `eval = FALSE` by default (tutorial mode, since the data path is a placeholder). After setting the `path` variable to your `.d` folder, run each code chunk individually in RStudio to better understand the process, or switch the setup chunk to `eval = TRUE` and render with the "Knit" button or `rmarkdown::render('workflow.Rmd')`.

# Workflow Steps

The OSSpRe workflow consists of the following major steps:

1.  **Raw Data Processing:** This step processes the raw MSI data to generate a reference peak list.
2.  **Reference Peak Picking:** This step identifies reference peaks from the super pixel data.
3.  **Quantitative Peak List Generation:** This step extracts peak intensities for all pixels based on the reference peaks and normalizes the data.
4.  **Qualitative Peak List Generation:** This step annotates the peaks using public databases for lipids and metabolites.
5.  **Exploratory Data Analysis:** This step includes data summary, peak statistics, visualization, segmentation, ion clustering, and reactomics analysis.

# HPC Usage

For large datasets, it is recommended to run the workflow on a High-Performance Computing (HPC) cluster. The following scripts are provided for HPC usage:

-   `refpeakpicking.R`: An R script to extract reference peaks from a super pixel.
-   `peakpicking.R`: An R script to extract peak profiles with a reference peak list.
-   `msi.job`: A Slurm job script for the peak picking process.

*Please check the code and update the paths for your own data.*

# Testing

A self-contained Colab notebook `OSSpRe_colab_test.ipynb` verifies the workflow on Linux: C++ compilation, peak alignment correctness and performance, ion mobility conversion against the Bruker library, and the full `helper.r` function suite, with the R environment restored from `renv.lock`. Open it directly in Colab (takes about 40-60 minutes, mostly compiling Bioconductor packages from source):

https://colab.research.google.com/github/yufree/OSSpRe/blob/main/OSSpRe_colab_test.ipynb

# Full Workflow on Real Data (Colab)

`OSSpRe_workflow.ipynb` runs the **complete workflow on a real Bruker `.d`
dataset**: it mounts Google Drive, restores the R environment from `renv.lock`,
then walks through every `workflow.Rmd` step — reference peaks, 2D peak
picking, the quantitative peak list, lipid/metabolite annotation, the
exploratory analyses (segmentation, ion clustering, reactomics) and molecular
networking — writing and displaying the results. Edit `/content/config.R` in
the notebook to point at your own `.d` and ROI:

https://colab.research.google.com/github/yufree/OSSpRe/blob/main/OSSpRe_workflow.ipynb

# Docker Usage

OSSpRe workflow is also available through the [xcmsrocker](https://github.com/yufree/xcmsrocker) project via the specialized `deepspace` container image.

1.  Install Docker and run it on your system.
2.  Pull the image from DockerHub: `docker pull yufree/deepspace:latest`
3.  Run the image: `docker run -e PASSWORD=xcmsrocker -p 8787:8787 yufree/deepspace:latest`
4.  Open your browser and go to `http://localhost:8787` or `http://[your-ip-address]:8787`.
5.  Log in to RStudio Server with username `rstudio` and password `xcmsrocker`.
6.  To use the OSSpRe workflow, go to `File -> New File -> R Markdown...`, select `From Template`, choose `Mass Spectrometry Imaging Workflow`, and click `OK`.