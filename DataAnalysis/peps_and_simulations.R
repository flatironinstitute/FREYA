#!/usr/bin/env Rscript

## Kiley Graim
## March 2018

## Generates Supplementary Table 1 and Supplementary Fig S5 
## Creates the PEP lists and runs PEP list generation simulations.
## IMPORTANT: REQUIRES run_genPEPs.R

# The Multiple-Histology-Per-Patient Design is Required to Effectively Capture PEPs
# We explicitly designed this study to take advantage the fact that CMTs frequently present with multiple primary tumors; we selected dogs presenting with at least one tumor of each histologic group.
# In order to test the importance of this type of design, we performed simulations using two types of subsampling approaches.
# The first approach mirrors our own study: a total of 30 samples from 10 patients (each with at least one normal, adenoma, and carcinoma sample) are selected at random.
# The second approach simulates cases where only two histologic categories can be gathered from a patient (as is typical in normal versus disease studies).
# Using these two subsample groups, we repeat the steps used to generate the PEPs and compare the maximum q-values characterizing the PEPs to those of the complete run (See Methods).
# Figure S5 presents Spearman correlations of 300 simulation runs for each subsample group.
# For all three PEPs, the three-histology approach significantly outperformed the paired approach (Wilcoxon rank-sum test p < 7e-10) when using a comparable number of samples.
# This result reiterates the importance of using patients presenting with all three histologic groups to capture these expression patterns, and underscores the power of the canine model.


############
### Set up the command line arguments
########

if(!require('getopt')) {
  install.packages('getopt')
  library(getopt)
}

## usage, options and doc 
argspec <- paste(get_Rscript_filename(), c(paste('Generates PEP lists and runs simulation tests.

  Usage: 
    ',get_Rscript_filename(),' -s <phenotype data> -m <mapping file>
  Options:
    -o <output directory>     Directory to write/store results
    -i <iterations>           Number of simulations to run
    -d <data directory>       Directory containing any input files
    -c <counts directory>     Directory containing counts files
        ')))

args <- commandArgs(TRUE)

## Print help if requested
if( length(args)==0 ) { args = '--help' } # If run without arguments, assume the user wants help
if ( '--help' %in% args | '-h' %in% args ) {
  write(argspec, stderr())
  quit()
}

## Set up input specs (long flag, short flag, required/optional, type)
spec <- matrix( c(
       'samplesCanine',            's', 1, 'character',
       'iterations',               'i', 2, 'integer',
       'datadir',                  'd', 2, 'character',
       'mapping',                  'm', 1, 'character',
       'outdir',                   'o', 2, 'character',
       'countdir',                 'c', 2, 'character'
      ),
    ncol=4,
    byrow=TRUE
         )

opt <- getopt( spec=spec )

## Set up reasonable defaults
if( is.null(opt$iterations) ) { opt$iterations = 300 }
if( is.null(opt$datadir) ) { opt$datadir = './' }
if( is.null(opt$countdir) ) { opt$countdir = '../' }
if( is.null(opt$outdir) ) { opt$outdir = './results' }

if( is.null(opt$samplesCanine) ) {
  print("ERROR: samplesCanine is a required argument"); flush.console()
  quit(save='no',status=1)
}


## If the output directory doesn't exist, create it
if(!dir.exists(opt$outdir)) {
  print(paste('Creating output directory',opt$outdir))
  system(paste('mkdir -p',opt$outdir))
}


############
### Functions for the script
########

## Given a patient number, return 1 sample ID of each histology type
pick.3.samples <- function(pat.num) {
  dat.pat <- dat.hist[dat.hist$PatientNumber==pat.num,]
  return( list(
    N=sample(dat.pat[dat.pat$Hist=='N','Qlabel'],1),
    B=sample(dat.pat[dat.pat$Hist=='B','Qlabel'],1),
    M=sample(dat.pat[dat.pat$Hist=='M','Qlabel'],1)
  ) )
}

## Given a list of patient numbers, return 2 samples for each. 
##    For half of the patients return 1 N and 1 B histology sample, for the other half return 1 N and 1 M histology sample
pick.2.samples <- function(pat.num) {

  # Split into 2 groups: ones with N&M histology and ones with N&B histology
  # Randomly sort the patient list since we use all patients every time
  pat.num <- sample(pat.num, length(pat.num), replace=F)
  pat.num.nb <- pat.num[1:(length(pat.num)/2)]
  pat.num.nm <- pat.num[(floor(length(pat.num)/2)):length(pat.num)]

  dat.pat.NB <- dat.hist[dat.hist$PatientNumber %in% pat.num.nb,]
  dat.pat.NM <- dat.hist[dat.hist$PatientNumber %in% pat.num.nm,]

  s.nb <- sapply(pat.num.nb, function(x) {
    dat.pat.NB <- dat.hist[dat.hist$PatientNumber %in% x,]
    list(
      sample(dat.pat.NB[dat.pat.NB$Hist=='N','Qlabel'],1),
      sample(dat.pat.NB[dat.pat.NB$Hist=='B','Qlabel'],1))
    } )
  s.nm <- sapply(pat.num.nm, function(x) {
    dat.pat.NM <- dat.hist[dat.hist$PatientNumber %in% x,]
    list(
      sample(dat.pat.NM[dat.pat.NM$Hist=='N','Qlabel'],1),
      sample(dat.pat.NM[dat.pat.NM$Hist=='M','Qlabel'],1))
    } )
  return( unlist(list(unlist(s.nb), unlist(s.nm))) )
}

############
### Begin analysis
########

## Load the function to generate PEP lists.
##   NOTE: This takes some time to load
#load(paste(opt$datadir,'humanmapping.rda', sep='/')) #TODO
load(opt$mapping) # TODO
source('run_genPEPs.R') 

## Load phenotype data for the dogs
## Make sure the required phenotype columns are in the data, quit if any are missing
print('Loading phenotype data...');flush.console()
print(paste(opt$datadir,opt$samplesCanine, sep='/'))
dat.hist <- read.csv(paste(opt$datadir,opt$samplesCanine, sep='/'))
if( all(c('Qlabel','Hist') %in% colnames(dat.hist)) ) {
  print(paste('Success,',nrow(dat.hist),'dog samples loaded.'));flush.console()
} else {
  print('ERROR: Qlabel and Hist columns required in the phenotype data.')
  quit(save='no',status=1) 
}

## Create the 'true' PEP lists using the full data
print('Calculating PEPs.'); flush.console()
peps.real <- gen_PEPs(dat.hist, opt$countdir)
qlim <- 0.05 # Minimum value to be included in a given PEP list
peps <- list( Adenoma=peps.real[peps.real$Adenoma_Expression_Pattern < qlim,'HumanSymbol'], Carcinoma=peps.real[peps.real$Carcinoma_Expression_Pattern < qlim,'HumanSymbol'], Tumor=peps.real[peps.real$Tumor_Expression_Pattern < qlim,'HumanSymbol'])
print('Full PEP lengths:'); flush.console()
print(sapply(peps, length)); flush.console() # Print num genes in each PEP 
write.table(peps.real, file=paste(opt$outdir,'CMT_PEPs.csv',sep='/'), sep=',', col.names=TRUE, row.names=FALSE, quote=FALSE)


## Run simulations using subsets of the data -- 2 versions of simulations will be run

## first approach: a total of 30 samples from 10 patients (each with at least one normal, adenoma, and carcinoma sample) are selected at random.
print('Running PEP simulations version 1.'); flush.console()
n.iters <- opt$iterations
print(paste('Running',n.iters,'iterations for each version.')); flush.console()
peps.cor.3hist <- matrix(NA, nrow=n.iters, ncol=3)
colnames(peps.cor.3hist) <- c('Tumor_Expression_Pattern','Carcinoma_Expression_Pattern','Adenoma_Expression_Pattern')

for(i in 1:n.iters) {
  ## Pick 10 random patients and 10 random samples(1 of each N,B,M from each patient)
  ids <- sort(sample(unique(dat.hist$PatientNumber), 10, replace=FALSE))
  s.ids <- unlist(sapply(ids, pick.3.samples))

  ## Generate the PEPs again
  peps.new <- gen_PEPs(dat.hist[dat.hist$Qlabel %in% s.ids,], opt$countdir)

  ## Calculate Spearman correlation
  ids.genes <- intersect( peps.new$EnsGene, peps.real$EnsGene )
  peps.cor.3hist[i,] <- sapply(colnames(peps.cor.3hist), function(x) { cor(peps.new[peps.new$EnsGene %in% ids.genes,x], peps.real[peps.real$EnsGene %in% ids.genes,x], method='spearman', use='complete') } )
  cat('.')
}
print('done!'); flush.console()

# second approach: simulates cases where only two histologic categories can be gathered from a patient (as is typical in normal versus disease studies).
# we used one normal and one adenoma sample from each of 8 randomly selected patients and 
#   used one normal and one carcinoma from each of the remaining patients, resulting in 32 samples per simulation.
print('Running PEP simulations version 2.'); flush.console()
peps.cor.2hist <- matrix(NA, nrow=n.iters, ncol=3)
colnames(peps.cor.2hist) <- c('Tumor_Expression_Pattern','Carcinoma_Expression_Pattern','Adenoma_Expression_Pattern')
for(i in 1:n.iters) {
  ## Pick 10 random patients and 10 random samples(1 of each N,B,M from each patient)
  ids <- unique(dat.hist$PatientNumber)
  s.ids <- pick.2.samples(ids)

  ## Generate the PEPs again
  peps.new <- gen_PEPs(dat.hist[dat.hist$Qlabel %in% s.ids,], opt$countdir)

  ## Calculate Spearman correlation
  ids.genes <- intersect( peps.new$EnsGene, peps.real$EnsGene ) 
  peps.cor.2hist[i,] <- sapply(colnames(peps.cor.2hist), function(x) { cor(peps.new[peps.new$EnsGene %in% ids.genes,x], peps.real[peps.real$EnsGene %in% ids.genes,x], method='spearman', use='complete') } )
  cat('.')
}
print('done!'); flush.console()

save(peps.cor.2hist, peps.cor.3hist, file=paste(opt$outdir,'Hist_Sims.RData',sep='/'))

peps.cor.2hist <- as.data.frame(peps.cor.2hist)
peps.cor.3hist <- as.data.frame(peps.cor.3hist)
peps.cor.2hist$Hist <- rep('Hist2', nrow(peps.cor.2hist))
peps.cor.3hist$Hist <- rep('Hist3', nrow(peps.cor.3hist))

# Wilcoxon rank-sum test of comparable number of samples
print("Calculating Wilcoxon rank-sum tests on each pep.");flush.console()
#print(paste(colnames(peps.cor.2hist)[1:3], sapply(colnames(peps.cor.2hist)[1:3], function(x) {wilcox.test(peps.cor.2hist[,x], peps.cor.3hist[,x], na.ignore=T)$p.value} )));flush.console()
print(paste('Wilcox rank-sum test across all PEPs:', wilcox.test(unlist(peps.cor.2hist[,1:3]), unlist(peps.cor.3hist[,1:3]), na.ignore=T)$p.value) );flush.console()

### Make the simulations boxplot

## Color Palette for manuscript
print('Setting color palettes.'); flush.console()
cols.palette <- c('#9DC7D8','#7FA1BE','#EBDA8C','#01B3CA','#4F6E35','#965354','#7DD1B9','#808040','#C6CBCE','#1D4E88','#C78C6C','#F69256','#D2B29E','#8B868C','#E38691','#B490B2')
cols.hist    <- c('#7DD1B9','#EBDA8C','#965354') # order = healthy, benign, malignant
cols.peps    <- c('#7FA1BE','#F69256','#E38691') # order = tumor, adenoma, carcinoma
cols         <- cols.palette[! (cols.palette %in% c(cols.hist, cols.peps)) ] # Don't resue the pep/histology colors

library(ggplot2)
print('Generating simulations correlation plot.'); flush.console()
tmp <- melt( rbind(peps.cor.2hist, peps.cor.3hist) )
ggplot(tmp) + geom_boxplot((aes(x=variable,y=value,fill=Hist))) + coord_flip() + scale_fill_manual(values=sample(cols,2)) + ylab('Spearman Correlation to true PEPs') + xlab('') + theme_bw(base_size = 18) 
ggsave(paste(opt$outdir,'PEP_hist_simulations.pdf',sep='/'), width=8,height=4) 

print('Finished! Success!'); flush.console()
