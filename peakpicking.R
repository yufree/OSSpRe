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

# run the function
getqlist(refpath = refpath, 
         lib_path = lib_path, 
         path = path, 
         method = method, 
         zero_proportion_cutoff = zero_proportion_cutoff, 
         coordpath = coordpath, 
         normpath = normpath,
         xrange = xrange,
         yrange = yrange)
