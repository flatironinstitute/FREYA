#!/usr/bin/env Rscript

## Load the BRCA data, subset to the genes in canine (all genes, not just Carcinoma&Tumor PEPs
## Calculate differential expression, make table of PEP memberships
## Plot Figure 5a
## Then print to file
##
## Kiley Graim
## Created: Jan 2017
## Updated: Dec 2017, converted to standalone script

#########################################
## Set up command line arguments
#####################################
if(!require('getopt')) {
  install.packages('getopt')
  library(getopt)
}

## usage, options and doc 
argspec <- paste(get_Rscript_filename(), c(paste('compares differential expression in human tumor/normal samples against the PEP genes from dog.

  Usage: 
    ',get_Rscript_filename(),'-p <PEPs filename>

  Options:
    -o <output directory>     Directory to write/store enrichment results

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
       'PEP',            'p', 1, 'character',
       'outdir',         'o', 2, 'character'
      ),
    ncol=4,
    byrow=TRUE
         )

opt <- getopt( spec=spec )

# Set defaults for optional parameters
if( is.null(opt$outdir) ) { opt$outdir = './results' }

## Set the arguments to easy-to-read names
fn.peps <- opt$PEP      # './data/CMT_peps.csv'  # PEP lists created in the expression pipeline. This is also Supp Table 1 in the manuscript 
out.dir  <- opt$outdir


#########################################
## Load all data needed for this R script
#####################################

## Filenames that don't change between users
data.dir <- './data/' # Working directory - we should be providing this with the requisite files # TODO: This will need to match the layout we give the whole pipeline. Have it be wherever we store the data
fn.brca.seq   <- paste0(data.dir,'BRCA_rnaseq_paired_noMets.t.txt')
fn.brca.pam50 <- paste0(data.dir,'BRCA_tumorVSnormal_paired.txt')
fn.geneconvs  <- paste0(data.dir,'Canine_Human_Gene_Conversion.txt') # Human-dog conversion - TODO should read this from the dog genom version (see Dima's code)
fn.cosmic       <- paste0(data.dir,'roles_in_cancer.csv') # COSMIC genes list, should download most recent version instead of using included one?

## Make sure all of the required files exist - quit if any are missing
for( fn in c(fn.peps, fn.brca.seq, fn.brca.pam50, fn.geneconvs, fn.cosmic) ) {
  if(!file.exists(fn)) { print(paste('ERROR: Unable to locate',fn)); quit(save='no',status=1) }
}

## Load libraries
print('Loading packages...');flush.console()
if(!require('samr')) {
 install.packages('samr')
  library(samr)
}
if(!require('ggplot2')) {
  install.packages('ggplot2')
  library(ggplot2)
}
if(!require('plyr')) {
  install.packages('plyr')
  library(plyr)
}

## Load the PEPs
peps <- read.table(fn.peps, sep=',', header=TRUE, stringsAsFactors=FALSE)
pep.genes <- c( peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], peps[peps$Tumor_Expression_Pattern < 0.05,'HumanSymbol'], peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'] )
pep.lists <- data.frame(Gene=pep.genes, PEP=rep(NA, length(pep.genes)))
pep.lists[ pep.lists$Gene %in% peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], 'PEP'] <- 'Adenoma'
pep.lists[ pep.lists$Gene %in% peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'], 'PEP'] <- 'Carcinoma'
pep.lists[ pep.lists$Gene %in% peps[peps$Tumor_Expression_Pattern < 0.05,'HumanSymbol'], 'PEP'] <- 'Tumor'
peps <- pep.lists

## Commented out for now- not sure this is the best place to be investigating the PEPs themselves.
## If you want to check overlaps in the lists 
#if(FALSE) {
#  length(intersect( peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], peps[peps$Tumor_Expression_Pattern < 0.05,'HumanSymbol'] ))
#  length(intersect( peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'] ))
#  length(intersect( peps[peps$Tumor_Expression_Pattern   < 0.05,'HumanSymbol'], peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'] ))
#}

## Load the data
dat <- read.table(fn.brca.seq, sep='\t', header=TRUE, row.names=1, check.names=FALSE)
dat <- as.matrix(dat) + 1
dat <- log(dat, base=2)

## Load the labels and subset to the labeled data
labs <- read.table(fn.brca.pam50, header=FALSE, row.names=1, check.names=FALSE)
ids <- intersect( rownames(labs), rownames(dat) )
dat <- dat[ids,]
labs <- labs[ids,]

## Subset to homologous genes
genes.all <- read.table(fn.geneconvs, sep='\t', header=TRUE, stringsAsFactors=FALSE)

## Subset PEP list to those also in human data, for our stats
dat.ovlp <- t(dat[,colnames(dat) %in% genes.all$Hum_Symb])
peps     <- peps[peps$Gene %in% rownames(dat.ovlp),]

## Run SAM
sam.seed <- sample(seq(1,1e5),1) 
res <- SAM(dat.ovlp, labs, resp.type='Two class unpaired', genenames=rownames(dat.ovlp), logged2=TRUE, nperm=100, testStatistic='wilcoxon', fdr.output=5e-6, random.seed=sam.seed)
print(paste('Using seed:', sam.seed)); flush.console()

## Get counts for sig up/down genes in each of the PEP lists
sig.genes <- as.data.frame(rbind(res$siggenes.table$genes.up, res$siggenes.table$genes.lo))
colnames(sig.genes)[6] <- 'FoldChange'
sig.genes$FoldChange <- as.numeric(as.character(sig.genes$FoldChange))
rownames(sig.genes) <- sig.genes[,'Gene ID']

sig.genes$PEP <- rownames(sig.genes) %in% peps$Gene
sig.genes$Log10FoldChange <- log10(sig.genes$FoldChange)
sig.genes$AbsLog10FoldChange <- abs(sig.genes$Log10FoldChange)

peps$AbsLog10FoldChange <- sig.genes[ peps$Gene , 'AbsLog10FoldChange']

## Add non-PEP rows from sig.genes to peps, that way we can have genes in 2+ peps
ids.nonPEP <- rownames(sig.genes)[!rownames(sig.genes) %in% peps$Gene]
tmp  <- data.frame( AbsLog10FoldChange=sig.genes[ids.nonPEP,'AbsLog10FoldChange'], PEP='nonPEP', row.names=ids.nonPEP )
tmp  <- data.frame( Gene=ids.nonPEP, PEP='nonPEP', AbsLog10FoldChange=sig.genes[ids.nonPEP,'AbsLog10FoldChange'] )
peps <- rbind(peps, tmp)
rm(tmp)


#########################################
## Color Palette for manuscript
#####################################
print('Setting color palettes.'); flush.console()

cols <- c('#9DC7D8','#7FA1BE','#EBDA8C','#01B3CA','#4F6E35','#965354','#7DD1B9','#808040','#C6CBCE','#1D4E88','#C78C6C','#F69256','#D2B29E','#8B868C','#E38691','#B490B2') # All colors in palette
cols.hist <- c('#7DD1B9','#EBDA8C','#965354') # order = healthy, benign, malignant
cols.peps <- c('#7FA1BE','#F69256','#E38691') # order = tumor, adenoma, carcinoma



#ggplot(na.omit(peps)) +  # Use if you want to show the nonPEP bar (colors will be wrong!)
ggplot(na.omit(peps[peps$PEP!='nonPEP',])) + 
  geom_boxplot(aes(PEP, AbsLog10FoldChange, fill=PEP) ,outlier.colour = 'white', outlier.alpha=0) + 
  scale_fill_manual(values=c(cols.peps[c(1,3,2)],'grey')) +
  coord_flip() + 
  ylim(0,1) +
  theme_minimal() +
  theme(axis.title = element_text(size=20), axis.text  = element_text(vjust=0.5, size=16)) +
  geom_hline( aes(yintercept=median(peps[ peps$PEP=='nonPEP', 'AbsLog10FoldChange'])),color='black',linetype='dotdash',lwd=1.5 ) + # Non-PEP genes 
  guides(fill = "none")
ggsave(paste(out.dir,'PEP_BRCA_AbsLogFoldChange.pdf',sep='/'),width=8,height=6)

mean(sig.genes$AbsLog10FoldChange)

print('Wilcoxon text absolute log10 fold change of each PEP'); flush.console()
print(wilcox.test(AbsLog10FoldChange ~ PEP=='Adenoma', data = peps)); flush.console()
print(wilcox.test(AbsLog10FoldChange ~ PEP=='Carcinoma', data = peps)); flush.console()
print(wilcox.test(AbsLog10FoldChange ~ PEP=='Tumor', data = peps)); flush.console()


### Print out the list of oncogenes/TSGs, PEP membership, SAM score, and fold change
## This is stored in Supp Table 5 for now.
print('Looking for oncogenes and tumor suppressor genes in the PEPs.'); flush.console()
genes.drivers <- read.table(fn.cosmic, sep=',', header=TRUE, stringsAsFactors=FALSE)
genes.drivers[genes.drivers$Role.in.Cancer=='','Role.in.Cancer'] <- NA
genes.drivers$Role.in.Cancer <- as.factor(genes.drivers$Role.in.Cancer)
colnames(genes.drivers)[1]   <- 'Gene'

peps <- peps[peps$Gene %in% rownames(sig.genes),] # Only keep differentially expressed genes
diff.eq <- join( genes.drivers, peps, by='Gene', type='inner')
diff.eq <- diff.eq[complete.cases(diff.eq),]

## NOTE: This is the results from when we ran it. I don't have the seed anymore :/
#diff.eq <- read.table('/Users/kgraim/Documents/CanineCancer/BRCA/TSG_Oncogene_SAM_Scores_by_PEP.csv', sep=',', header=T) 

## Make the table
tbl <- table( diff.eq$Role.in.Cancer, diff.eq$PEP )
tbl <- tbl[,which(colnames(tbl)%in%c('Adenoma','Carcinoma','Tumor'))] # Don't show non-PEP groups
print(tbl); flush.console()

print('Storing to file'); flush.console()
write.table(diff.eq[diff.eq$PEP!='nonPEP',], file=paste(out.dir,'PEP_Roles_in_Cancer_summary.csv',sep='/'), col.names=TRUE, row.names=FALSE, quote=FALSE, sep=',')

## Plot & save
pdf(paste(out.dir,"PEP_Roles_in_Cancer_Barplot.pdf",sep='/'))
barplot(tbl, legend=TRUE, ylim=c(0,max(apply(tbl,2,sum))+8), col=cols[c(4,5,7)]) # Add vertical space for legend
dev.off()

print('Success!')
