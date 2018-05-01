#!/usr/bin/env Rscript

## NOTE: I pulled out the PAM50 BRCA label wrangling, since it only needed to be done once. Code is in: compile_subtypes.R and output is in: BRCA_PAM50_labels.csv

## NOTE: This works on rusty but not on Kiley's local mac. 

# Example: ./PAM50_refactored.R --clinical=./data/UserDogData_Phenotype.csv --expression=./data/UserDogData_RNASeq.tab

#########################################
## Set up command line arguments
#####################################
if(!require('getopt')) {
  install.packages('getopt')
  library(getopt)
}

require(methods) # glmnet needs this but it isn't loaded by default in Rscript

## usage, options and doc 
argspec <- paste(get_Rscript_filename(), c(paste('predicts PAM50 subtypes based on TCGA BRCA RNA-Seq data and user-provided CMT RNA-Seq data.

  Usage: 
    ',get_Rscript_filename(),' -c <dog clinical/phenotype data> -e <dog expression data> -p <PEPs list>
  Options:
    -d <integer>      Delimiter used in CMT data files
    -o <directory>    Directory for results files
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
       'clinical',       'c', 1, 'character',
       'expression',     'e', 1, 'character',
       'PEP',            'p', 1, 'character',
       'outdir',         'o', 2, 'character',
       'delim',          'd', 2, 'character'
      ),
    ncol=4,
    byrow=TRUE
         )

opt <- getopt( spec=spec )

# Set defaults for optional parameters
#   Use tab delimiter if not provided
if( is.null(opt$delim) ) { opt$delim = '\t' }
if( is.null(opt$outdir) ) { opt$outdir = './results' }


## Set the arguments to easy-to-read names
data.dir <- './data/'      # Working directory - we should be providing this with the requisite files # TODO: This will need to match the layout we give the whole pipeline. Have it be wherever we store the data
fn.hist  <- opt$clinical   # './data/UserDogData_Phenotype.csv' # User-provided
fn.dd    <- opt$expression # './data/UserDogData_RNASeq.tab'    # Dog RNASeq expression data filename
fn.peps  <- opt$PEP        # './data/CMT_peps.csv'  # PEP lists created in the expression pipeline. This is also Supp Table 1 in the manuscript 
fn.sep   <- opt$delim      # delimiter in their files, default should be tab
fn.out   <- opt$outdir     # Results folder

## Filenames that don't change between users
fn.pam50.genes  <- paste0(data.dir,'PAM50_genes.csv')       # This list will never change, no need to have as input
fn.brca         <- paste0(data.dir,'BRCA_rnaseq_paired_noMets.t.txt') # TCGA paired BRCA RNA-Seq data
fn.brca.pam50   <- paste0(data.dir,'BRCA_PAM50_labels.csv') # TCGA BRCA PAM50 labels, from: http://www.cell.com/cell/fulltext/S0092-8674(15)01195-2

## Make sure all of the required files exist - quit if any are missing
for( fn in c(fn.hist, fn.dd, fn.peps, fn.pam50.genes, fn.brca, fn.brca.pam50) ) {
  if(!file.exists(fn)) { print(paste('ERROR: Unable to locate',fn)); quit(save='no',status=1) }
}


#########################################
## Analysis!
#####################################

## Load libraries, install them if they aren't already
print('Loading packages...');flush.console()
if(!require('glmnet')) {
  install.packages('glmnet')
  library(glmnet)
}
if(!require('sva')) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("sva")
  library(sva)
}
if(!require('gplots')) {
  install.packages('gplots')
  library(gplots)
}

## PAM50 genes list - this never changes
genes           <- rownames(read.table(fn.pam50.genes, sep=',', row.names=1))

## Load PEP lists
print('Loading PEPs. PEP list lengths:'); flush.console()
peps <- read.table(fn.peps, sep=',', header=TRUE, stringsAsFactors=FALSE)
peps <- list( Adenoma=peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], Carcinoma=peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'], Tumor=peps[peps$Tumor_Expression_Pattern < 0.05,'HumanSymbol'])
print(sapply(peps, length)); flush.console() # Print num genes in each PEP 

## Load the BRCA RNASeq data, extract tumor/normal labels & rescale
##    From: http://www.cell.com/cell/fulltext/S0092-8674(15)01195-2
print('Loading human data...');flush.console()
dat           <- read.table(fn.brca, sep='\t', header=T, row.names=1, check.names=F)
dat           <- dat[,genes] # Subset to PAM50 genes
dat           <- log(dat+0.01, base=2) # Log scale
dat.status    <- substr(rownames(dat),14,15) # Tumor/Normal labels
rownames(dat) <- substr(rownames(dat),1,15) # Remove the 'A' from labels
print(paste('Success,',nrow(dat),'human samples loaded.'));flush.console()

## Create the PAM50 labels
pam50.labs <- read.table(fn.brca.pam50, sep=',', header=TRUE, row.names=1, check.names=FALSE)

## Subset expression & labels data to overlapping sets
ids        <- intersect( rownames(dat), rownames(pam50.labs) )
dat        <- dat[ids,]
pam50.labs <- pam50.labs[ids,,drop=FALSE]

## Remove NA subtypes
if( any(is.na(pam50.labs)) ) {
  ids.na     <- rownames(pam50.labs)[which(!is.na(pam50.labs))] # Names for non-NA labels
  dat        <- dat[ids.na,]
  pam50.labs <- pam50.labs[ids.na,,drop=FALSE]
}

## Load the canine data 
print('Loading dog data...');flush.console()
dat.dog <- read.table(fn.dd, header=TRUE, row.names=1, check.names=FALSE, sep=fn.sep)
print(paste('Success,',nrow(dat.dog),'dog samples loaded.'));flush.console()

## Drop both data matrices to the overlapping genes
genes     <- genes[genes %in% colnames(dat.dog)]
dat.dog   <- dat.dog[,genes]
dat.pam50 <- dat[,genes]

## Phenotype data for the dogs
dog.clin <- read.table(fn.hist, sep=fn.sep, header=T, row.names=1) # This should be user input

## Subset dog expression & phenotype data to overlapping samples
ids <- intersect( rownames(dog.clin), rownames(dat.dog) )
if(  !all(rownames(dog.clin) %in% ids & rownames(dat.dog) %in% ids) ) { 
  print("WARNING: Dog expression and phenotype data have unmatched samples!");flush.console()
  dat.dog  <- dat.dog[ids,]
  dog.clin <- dog.clin[ids,]
}


#########################################
## Color Palette for manuscript
#####################################
print('Setting color palettes.'); flush.console()
 
cols.palette <- c('#9DC7D8','#7FA1BE','#EBDA8C','#01B3CA','#4F6E35','#965354','#7DD1B9','#808040','#C6CBCE','#1D4E88','#C78C6C','#F69256','#D2B29E','#8B868C','#E38691','#B490B2')
cols.hist    <- c('#7DD1B9','#EBDA8C','#965354') # order = healthy, benign, malignant
cols.peps    <- c('#7FA1BE','#F69256','#E38691') # order = tumor, adenoma, carcinoma
cols         <- cols.palette[! (cols.palette %in% c(cols.hist, cols.peps)) ] # Don't resue the pep/histology colors


#########################################
## Remove batch effects
#####################################
print('Training models, making predictions, etc...');flush.console()

## Remove batch effects
dat.both     <- rbind( dat.dog, dat.pam50 )
pheno        <- as.factor(c( rep('Dog', nrow(dat.dog)), rep('Human', nrow(dat.pam50)) ) )
names(pheno) <- c(rownames(dat.dog), rownames(dat.pam50))
combat.edata <- t(ComBat(dat=t(dat.both), batch=pheno))

## Mean center the data
combat.edata <- sweep(combat.edata, 1, rowMeans(combat.edata))


#########################################
## Train PAM50 model & predict subtypes for dog
#####################################

## Irain subtype predictor, then predict dog samples
print('Training the predictor');flush.console()
res             <- cv.glmnet(data.matrix(combat.edata[rownames(pam50.labs),]), pam50.labs[,1], family='multinomial')
print('Predicting new labels');flush.console()
preds           <- predict(res, combat.edata, type='class')
rownames(preds) <- rownames(combat.edata)

## Write predictions to file (dog only)
write.table( preds[rownames(dat.dog),,drop=FALSE], sep=',', col.names=FALSE, row.names=TRUE, quote=FALSE, file=paste(fn.out,'PAM50_dog.csv',sep='/'))

## Barplots of the predictions - not used in manuscript
print('User-check barplots');flush.console()
pdf(paste(fn.out,'PAM50_barplots.pdf',sep='/'))
barplot(table(preds[rownames(pam50.labs),1], pam50.labs[,1]), legend=TRUE, col=cols) # Human sample predictions
barplot( table(preds[rownames(combat.edata) %in% rownames(dog.clin),1], dog.clin$Hist ), legend=TRUE, col=cols ) # Dog sample predictions
dev.off()


#########################################
## Figure 6 in the manuscript - combined heatmaps
#####################################

## Identify PAM50 Genes in the PEP lists
genes.cols <- rep('grey60', ncol(combat.edata))
names(genes.cols) <- colnames(combat.edata)
genes.cols[intersect(names(genes.cols), peps$Tumor)]     <- cols.peps[1] # Tumor PEP genes
genes.cols[intersect(names(genes.cols), peps$Adenoma)]   <- cols.peps[2] # Adenoma PEP genes
genes.cols[intersect(names(genes.cols), peps$Carcinoma)] <- cols.peps[3] # Carcinoma PEP genes

## Colors for the PAM50 subtypes
preds.cols <- matrix(nrow=nrow(preds),ncol=ncol(preds)+1, dimnames=list(rownames(preds)))
preds.cols[preds[,1]=='Normal',1] <- cols.palette[4]
preds.cols[preds[,1]=='Her2',1]   <- cols.palette[13]
preds.cols[preds[,1]=='LumA',1]   <- cols.palette[1]
preds.cols[preds[,1]=='LumB',1]   <- cols.palette[10]
preds.cols[preds[,1]=='Basal',1]  <- cols.palette[9]

## Color bar for the species
preds.cols[names(pheno),2]        <- pheno
preds.cols[preds.cols[,2]=='1',2] <- cols.palette[10] # Purple
preds.cols[preds.cols[,2]=='2',2] <- cols.palette[3] # Dark green

print('Generating heatmaps');flush.console()

## Add line breaks for each PAM50 subtype gorup in the heatmaps
tbl.preds <- table(preds)
rowsep.nums <- c( sum(tbl.preds[1]), sum(tbl.preds[1:2]), sum(tbl.preds[1:3]), sum(tbl.preds[1:4]) )

# Heatmap with PAM50 Subtype color bar
pdf(paste(fn.out,'PAM50_combined_heatmap_subtypes_bar.pdf',sep='/'),width=8,height=20) 
heatmap.2(combat.edata[names(sort(preds[,1])),], trace='none', col=colorRampPalette(c('blue','blue','grey20','yellow','yellow')), RowSideColors=preds.cols[names(sort(preds[,1])),1], dendrogram='column', ColSideColors=genes.cols, Rowv='none', rowsep=rowsep.nums, sepcolor='grey60')
dev.off()
 
# Heatmap with species color bar
pdf(paste(fn.out,'PAM50_combined_species_bar.pdf',sep='/'),width=8,height=20)
heatmap.2(combat.edata[names(sort(preds[,1])),], trace='none', col=colorRampPalette(c('blue','blue','grey20','yellow','yellow')), RowSideColors=preds.cols[names(sort(preds[,1])),2], dendrogram='column', ColSideColors=genes.cols, Rowv='none', rowsep=rowsep.nums, sepcolor='grey60')
dev.off()

print('Finished with PAM50 analysis!');flush.console()

