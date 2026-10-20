library(data.table)
library(Rcpp)
source('helper.r')

# set the parameters
refpath <- 'ref2d.csv'
lib_path <- 'libtimsdata.so'
path <- 'PATH_TO_RAW.d'
method <- 'tic'
zero_proportion_cutoff <- 0.95
coordpath <- 'coord.csv'
normpath <- 'ticmzccs.csv'
xrange <- c(2900, 3000)
yrange <- c(0, 10000)
# set the uncertainty of mz and ccs (must match the values used in refpeakpicking.R)
mz_ppm <- 20
ccs_tolerance <- 0.03
batch_size <- 25000

# run the function
getqlist(refpath = refpath,
         lib_path = lib_path,
         path = path,
         method = method,
         zero_proportion_cutoff = zero_proportion_cutoff,
         coordpath = coordpath,
         normpath = normpath,
         xrange = xrange,
         yrange = yrange,
         mz_ppm = mz_ppm,
         ccs_tolerance = ccs_tolerance,
         batch_size = batch_size)
