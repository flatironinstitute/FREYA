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
  biocLite("sva")
  library(sva)
}
if(!require('gplots')) {
  install.packages('gplots')
  library(gplots)
}
if(!require('samr')) {
 biocLite("impute")
 install.packages('samr')
 library(samr)
}
if(!require('plyr')) {
  install.packages('plyr')
  library(plyr)
}
if(!require('topGO')) {
  biocLite("topGO")
  library(topGO)
}
if(!require('goseq')) {
  biocLite("goseq")
  library(goseq)
}
if(!require('dplyr')) {
  install.packages('dplyr')
  library(dplyr)
}
if(!require('knitr')) {
  install.packages('knitr')
  library(knitr)
}
if(!require('broom')) {
  install.packages('broom')
  library(broom)
}
if(!require('biobroom')) {
  biocLite("biobroom")
  library(biobroom)
}
if(!require('tidyr')) {
  install.packages('tidyr')
  library(tidyr)
}
if(!require('qvalue')) {
  biocLite("qvalue")
  library(qvalue)
}
if(!require('edgeR')) {
  biocLite("edgeR")
  library(edgeR)
}
if(!require('biomaRt')) {
  biocLite("biomaRt")
  library(biomaRt,pos = "package:base")
}
if(!require('reshape2')) {
  install.packages('reshape2')
  library(reshape2)
}
if(!require('devtools')) {
  install.packages('devtools')
  library(devtools)
}
if(!require('org.Hs.eg.db')) {
  biocLite("org.Hs.eg.db")
}
