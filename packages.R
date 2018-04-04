# commands to install needed R packages
# sudo R -f packages.R
install.packages('getopt', repos='http://cran.us.r-project.org')
install.packages('ggplot2', repos='http://cran.us.r-project.org')

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")

if(!require('glmnet')) {
  install.packages('glmnet')
  library(glmnet)
}
if(!require('sva')) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("sva")
  library(sva)
}

biocLite()

if(!require('gplots')) {
  install.packages('gplots')
  library(gplots)
}

if(!require('samr')) {
 install.packages('samr')
  library(samr)
}

if(!require('plyr')) {
  install.packages('plyr')
  library(plyr)
}
