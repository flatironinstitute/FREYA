#!/usr/bin/env Rscript

## Kiley Graim
## Created: August 2017
## Last Updated: Dec 2017

## Example: ./mut_analysis_refactored.R --clinical=./data/UserDogData_Phenotype.csv --mutation=./data/mutations_genesOnly.csv --PEP=./data/CMT_peps.csv --num=5



#########################################
## Set up command line arguments
#####################################
if(!require('getopt')) {
  install.packages('getopt')
  library(getopt)
}


## usage, options and doc 
argspec <- paste(get_Rscript_filename(), c(paste('runs mutations analysis on the CMT data.

  Usage: 
    ',get_Rscript_filename(),'-c <clinical/phenotype filename> -m <mutations filename>  -p <PEPs filename>
  Options:
    -n <integer>      (Minimum number of mutations per gene to be included in analysis (default 5)
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
       'mutation',       'm', 1, 'character',
       'PEP',            'p', 1, 'character',
       'num',            'n', 2, 'integer'
      ),
    ncol=4,
    byrow=TRUE
         )

opt <- getopt( spec=spec )

#set some reasonable defaults for optional parameters
if( is.null(opt$num) ) { opt$num = 5 }

## Set the arguments to easy-to-read names
fn.hist <- opt$clinical # './data/UserDogData_Phenotype.csv' # User-provided
fn.muts <- opt$mutation # './data/mutations_genesOnly.csv' # Created by the mutations pipeline (see shell scripts)
fn.peps <- opt$PEP      # './data/CMT_peps.csv'  # PEP lists created in the expression pipeline. This is also Supp Table 1 in the manuscript  

## TODO: Make the PEP analysis optional- if the user doesn't provide a PEP file, then don't run it!

#########################################
## Load all data needed for this R script
#####################################

## Filenames that don't change between users
data.dir <- './data/' # Working directory - we should be providing this with the requisite files # TODO: This will need to match the layout we give the whole pipeline. Have it be wherever we store the data
fn.pam50        <- paste0(data.dir,'PAM50_dog.csv') # Created by PAM50_refactored.R, used doesn't need to provide (but should be optional???)
fn.cosmic       <- paste0(data.dir,'genes_COSMIC.csv') # COSMIC genes list, should download most recent version instead of using included one?
fn.pam50.genes  <- paste0(data.dir,'PAM50_genes.csv') # This list will never change, no need to have as input

## Make sure all of the required files exist - quit if any are missing
for( fn in c(fn.hist, fn.muts, fn.peps, fn.pam50, fn.cosmic, fn.pam50.genes) ) {
  if(!file.exists(fn)) { print(paste('ERROR: Unable to locate',fn)); quit(save='no',status=1) }
}

## Load the clinical data, extract dog IDs
print('Loading script data'); flush.console()
dat.hist <- read.table(fn.hist, sep=',', header=TRUE, row.names=1)
if('9A' %in% rownames(dat.hist)) {rownames(dat.hist)[which(rownames(dat.hist)=='9A')] <- '9A1' } # TODO: Specific to our dataset only! Remove this line once we're done testing (or if we switch 9A1.bam to 9A.bam)

## Alphabetic IDs for each dog instead of numeric
## Generally don't need this- only CMTGA changes up the names halfway through
dat.hist$Patient <- as.character(dat.hist$Patient) # Ensure patient names are character strings for plotting consistency

## Load the PAM50 subtypes 
##   This file is created by PAM50_refactored.R, should just port straight over (don't need user to specify)
pam50 <- read.table(fn.pam50, sep=',', row.names=1)

print('PAM50 subtype counts per patient:'); flush.console()
print(table(pam50[,1], dat.hist[rownames(pam50),'Patient']))
dat.hist$PAM50 <- pam50[rownames(dat.hist),1]

## Load the list of COSMIC genes
genes.cosmic <- rownames(read.table(fn.cosmic, sep=',', header=TRUE, row.names=1))

## Load the mutations data, make 0/1 calls instead of # calls per gene
dat <- read.table(fn.muts, sep=',', header=TRUE, row.names=1, check.names=FALSE)

dat.bin <- dat
dat.bin[dat.bin>0] <- 1

## Create 2 matrices: Benign and Malignant samples
print('Creating the mutations matrices.'); flush.console()
dat.m <- dat.bin[,rownames(dat.hist)[dat.hist$Hist=='M']]
dat.b <- dat.bin[,rownames(dat.hist)[dat.hist$Hist=='B']]

dat.m <- t(aggregate(t(dat.m), by=list(dat.hist[colnames(dat.m),'Patient']), FUN=sum))
colnames(dat.m) <- dat.m[1,]
dat.m <- dat.m[-1,]

dat.b <- t(aggregate(t(dat.b), by=list(dat.hist[colnames(dat.b),'Patient']), FUN=sum))
colnames(dat.b) <- dat.b[1,]
dat.b <- dat.b[-1,]

## Convert from character to numeric
class(dat.b) <- 'numeric'
class(dat.m) <- 'numeric'

# For now we don't care about # samples mutated in each gene per patient, just that at least 1 sample is mutated
# So set >1 values to 1
dat.m[dat.m>0] <- 1
dat.b[dat.b>0] <- 1

#########################################
## Color Palette for manuscript
#####################################
print('Setting color palettes.'); flush.console()

cols <- c('#9DC7D8','#7FA1BE','#EBDA8C','#01B3CA','#4F6E35','#965354','#7DD1B9','#808040','#C6CBCE','#1D4E88','#C78C6C','#F69256','#D2B29E','#8B868C','#E38691','#B490B2') # All colors in palette
cols.hist <- c('#7DD1B9','#EBDA8C','#965354') # order = healthy, benign, malignant
cols.peps <- c('#7FA1BE','#F69256','#E38691') # order = tumor, adenoma, carcinoma


#########################################
## Figure 2b - red&blue histogram
#####################################
print('Generating mutation histogram'); flush.console()

require(ggplot2)

## Create data frame for patient summaries, converting to alphabet patient IDs instead of numeric
mut.rates <- data.frame(Muts=apply(dat.bin, 2, sum), Hist=dat.hist[colnames(dat.bin),'Hist'], Dog=dat.hist[colnames(dat.bin),'Patient'], Sample=colnames(dat.bin))

## Plot the subfigure
ggplot(mut.rates, aes(Dog, Muts)) + geom_bar(aes(fill = Hist), position = "dodge", stat="identity") + scale_fill_manual(values=cols.hist[2:3]) + theme_minimal() + coord_flip() + theme(axis.text.x=element_text(angle = -325, hjust = 1), text = element_text(size=30))
ggsave('Sample_Mut_Rates.pdf',width=4,height=10)


#########################################
## Figure 2a - red&blue density plot
#####################################
print('Generating density plot'); flush.console()

## Count mutations in each benign & malignant sample, create and save density plot
samples.freq <- data.frame(Mutations=apply(dat, 2, sum), Hist=dat.hist[colnames(dat),'Hist'])
ggplot(samples.freq) + geom_density(aes(Mutations,group=Hist,col=Hist),lwd=3) + scale_color_manual(values=cols.hist[2:3]) + theme_bw() + theme(text = element_text(size=20))
ggsave('Sample_Mutation_Counts_Density.pdf',width=12, height=4)

## Print the median number of mutated genes per histology
print( paste('Median mutations in benign samples:', median( samples.freq[samples.freq$Hist=='B','Mutations']) )); flush.console()
print( paste('Median mutations in malignant samples:', median( samples.freq[samples.freq$Hist=='M','Mutations']) )); flush.console()


#########################################
## Figure 2b - navy&white dot plot
#####################################
print('Generating pooled mutations plot'); flush.console()

require(reshape2)

## Calculate most frequently mutated (by % of samples of each type) to get balanced frequently mutated genes
##    Otherwise will give mostly benign mutations, since we have 2x benign samples
ids.benign <- rownames(dat.hist)[ dat.hist$Hist=='B' ]
ids.tumor <- rownames(dat.hist)[ dat.hist$Hist=='M' ]

benign.ratios <- apply(dat.bin[,ids.benign], 1, function(x){sum(x==1)/length(x)})
tumor.ratios <- apply(dat.bin[,ids.tumor], 1, function(x){sum(x==1)/length(x)})
max.ratios <- apply(cbind(benign.ratios,tumor.ratios), 1, max)

## INPUT: This should be an optional parameter (file with list of genes OR statistic to use for picking genes) with default 30 genes w/max ratios
genes <- names(sort(apply(dat[rownames(dat) %in% genes.cosmic,], 1, function(x){sum(x>0)}),decreasing=TRUE))[1:30] # Use this for the COSMIC plot
#genes            <- names(max.ratios)[max.ratios>0.15]  # Another option - pick some cutoff of mutated ratios for benign/malignant
#genes <- sample(rownames(dat)[!rownames(dat) %in% genes.cosmic], 30) #Use this for the random sampling plot (randomly samples from non-cosmic genes)

## Melt the malignant sample matrix
dat.m.melted <- melt( as.matrix(dat.m[rownames(dat.m) %in% genes,]) )
dat.m.melted$value <- as.numeric(as.character(dat.m.melted$value))
dat.m.melted$value[dat.m.melted$value>0] <- 1
dat.m.melted$value <- as.factor(dat.m.melted$value) # For color scales

## Melt the benign sample matrix
dat.b.melted <- melt( as.matrix(dat.b[rownames(dat.b) %in% genes,]) )
dat.b.melted$value <- as.numeric(as.character(dat.b.melted$value))
dat.b.melted$value[dat.b.melted$value>0] <- 1
dat.b.melted$value <- as.factor(dat.b.melted$value) # For color scales


## Combine the 2 melted matrices
dat.melted<- cbind(dat.m.melted, dat.b.melted$value)
colnames(dat.melted) <- c('Gene','Dog','Tumor','Benign')
dat.melted$Dog <-  as.character(dat.melted$Dog) # So the plot sorts them alphabetically

## Plot the result
ggplot(dat.melted) + geom_point(aes(Gene, Dog, col=Tumor), size=8, pch=15) +
  geom_point(aes(Gene,Dog,col=Benign), size=4, pch=16) +
  theme(axis.text.x=element_text(angle = -325, hjust = 1)) +
  scale_color_manual(values=c('white',cols[10]))
ggsave('COSMIC_Genes_Mutations.pdf', width=10, height=5.5)


#########################################
## Supplemental Figure 4 - Frequently mutated genes
#####################################

print('Generating per-sample mutations plot'); flush.console()

genes            <- names(max.ratios)[max.ratios>0.15]
s.counts         <- table(dat.hist[colnames(dat.bin),'Patient'])

# For our dataset only, reorder the names (because we used numeric patient names :/
# TODO: Nick once we've finished testing for our data, we should remove this next line
s.counts         <- s.counts[sort(paste0(names(s.counts),'A'),index=TRUE)$ix] 

dat.bin.melted           <- melt(as.matrix(dat.bin[genes,]))
dat.bin.melted           <- dat.bin.melted[,c(2,1,3)]
colnames(dat.bin.melted) <- c('Sample','Gene','Alteration')
dat.bin.melted$Sample    <- as.character(dat.bin.melted$Sample)
dat.bin.melted$Hist      <- dat.hist[ dat.bin.melted$Sample, 'Hist' ]

ggplot(dat.bin.melted) +
  geom_point(aes(Sample, Gene, color=interaction(factor(Alteration),Hist)),pch=15,size=3) +
  scale_color_manual(values=c('white',cols.hist[2],'white',cols.hist[3])) +
  theme_classic() +
  theme(legend.position='none',axis.text.x=element_text(angle = -325, hjust = 1)) +
  geom_vline(xintercept=cumsum(s.counts[-length(s.counts)])+0.5,col=cols[14],size=2)
ggsave('MutationConsistency.pdf',width=13,height=7)


#########################################
## Added after sharing with SF
#####################################

## Do the subtypes have different numbers of mutations (total, not just in PAM50 genes)
##   For samples of each subtype, print median num mutations in the samples
num.muts <- apply(dat.bin, 2, sum)
colnames(pam50)[1] <- 'PAM50'
pam50$Muts <- NA
pam50[colnames(dat.bin),'Muts'] <- num.muts
print('PAM50 sample counts:'); flush.console()
print(sapply( levels(pam50$PAM50), function(x) {median( pam50[pam50$PAM50==x,'Muts'], na.rm=TRUE )} )); flush.console() 

### Are COSMIC genes more frequently mutated than non-COSMIC?
print(res.ttest <- t.test( apply(dat.bin, 1, sum) ~ factor(rownames(dat.bin) %in% genes.cosmic) ))
if( res.ttest$p.value < 0.05) {
  print( paste('COSMIC genes are significantly more frequently mutated than non-COSMIC genes, p-value =', signif(res.ttest$p.value,digits=3)) )
} else {
  print( paste('COSMIC genes are NOT significantly more frequently mutated than non-COSMIC genes, p-value =', signif(res.ttest$p.value,digits=3)) )
}
flush.console()
rm(res.ttest)

### Are PAM50 genes more frequently mutated than non-PAM50?
genes.pam50 <- rownames(read.table(fn.pam50.genes, sep=',', row.names=1))
print( res.ttest <- t.test( apply(dat.bin, 1, sum) ~ factor(rownames(dat.bin) %in% genes.pam50) ) ) 
if( res.ttest$p.value < 0.05) {
  print( paste('PAM50 genes are significantly more frequently mutated than non-COSMIC genes, p-value =', signif(res.ttest$p.value,digits=3)) )
} else {
  print( paste('PAM50 genes are NOT significantly more frequently mutated than non-COSMIC genes, p-value =', signif(res.ttest$p.value,digits=3)) )
}
rm(res.ttest)
flush.console()

## Correlate mutations w/clinical factors of interest- this will return a matrix of dat.hist columns by genes, filled with corrected pvals 
print(paste('Calculating correlations between mutations and phenotype data (Patient, Location, Histology, etc) in genes with >',opt$num,'mutations in the cohort.'))
get.pvals <- function(id) {
  phen.cols <- c('Patient','Location','Goldschmidt','Hist','SimHist','DetHist','PAM50') # Which clinical factors we care about
  phen.cols <- phen.cols[ phen.cols %in% colnames(dat.hist) ] # Make sure these are in the provided phenotype/clinical data
  p.adjust(apply(dat.hist[colnames(dat.bin),phen.cols], 2, function(x) {try(chisq.test(table( factor(x), unlist(dat.bin[id,])))$p.value)}))
}
genes       <- names(which(apply(dat.bin, 1, sum)>opt$num)) # Only care about frequently mutated genes
genes.pvals <- sapply(genes, get.pvals)
write.table(signif(t(genes.pvals),digits=5), file='FreqMutatedGenes_ClinicalCorrelations.csv', sep=',', col.names=TRUE, row.names=TRUE, quote=FALSE)
print('Phenotype/Clinical correlations stored to file.')

## Are PEP list genes more frequently mutated?
## Load the PEPs & print PEP genes that are frequently mutated
print('Loading PEPs. PEP list lengths:'); flush.console()
peps <- read.table(fn.peps, sep=',', header=TRUE, stringsAsFactors=FALSE)

peps <- list( Adenoma=peps[peps$Adenoma_Expression_Pattern < 0.05,'HumanSymbol'], Carcinoma=peps[peps$Carcinoma_Expression_Pattern < 0.05,'HumanSymbol'], Tumor=peps[peps$Tumor_Expression_Pattern < 0.05,'HumanSymbol'])
print(sapply(peps, length)); flush.console() # Print num genes in each PEP 
print(paste('Checking for frequently mutated PEP genes (>',opt$num,'mutations):'))
pep.mut.counts <- sapply(peps, function(x) { apply(dat.bin[rownames(dat.bin) %in% x,], 1, sum)} )
print(sapply(pep.mut.counts, function(x){ names(x)[which(x>opt$num)] }))
print(sapply(pep.mut.counts, summary))
flush.console()

print('Done with mutation analysis.') 
