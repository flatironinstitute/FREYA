#!/usr/bin/env Rscript

## Kiley Graim
## March 2018
## 
## maps between human and dog, normalizes the data, etc

###########################
###   Script setup
####################

if(!require('getopt')) {
  install.packages('getopt')
  library(getopt)
}

## usage, options and doc
argspec <- paste(get_Rscript_filename(), c(paste('generates rda and delimited files for the canine data. This preprocessing script creates the input files for the rest of the analysis scripts (except CMT_PEPs.csv)

  Usage: 
    ',get_Rscript_filename(),' 
  Options:
    -o <output directory>     Directory to write/store results
    -m <mapping>              Set this flag to re-create the humanmapping.rda, human mapping file. Generally unnecessary.
    -d <data directory>       Directory containing any input files
        ')))

args <- commandArgs(TRUE)

## Print help if requested
if ( '--help' %in% args | '-h' %in% args ) {
  write(argspec, stderr())
  quit()
}

## Set up input specs (long flag, short flag, required/optional, type)
spec <- matrix( c(
       'outdir',     'o', 2, 'character',
       'mapping',    'm', 2, 'logical',
       'datadir',    'd', 2, 'character'
      ),
    ncol=4,
    byrow=TRUE
         )

opt <- getopt( spec=spec )

### set some reasonable defaults
if( is.null(opt$outdir) ) { opt$outdir = 'output/' }
if( is.null(opt$map.human) ) { opt$map.human = FALSE } # Set to TRUE if you want to regenerate the human mapping file.
if( is.null(opt$datadir) ) { opt$datadir = './' }

## If the output directory doesn't exist, create it
if(!dir.exists(opt$outdir)) {
  print(paste('Creating output directory',opt$outdir))
  system(paste('mkdir -p',opt$outdir))
}

library(knitr)
library(plyr)
library(dplyr)
library(broom)
library(biobroom)
library(tidyr)
library(qvalue)
library(edgeR)
library(biomaRt,pos = "package:base")
library(reshape2)
#library(bindrcpp)

library(devtools)
knitr::opts_chunk$set(warning = FALSE,
                      message = FALSE,
                      prompt = FALSE)

############
### Prepare Data & Identify Differentially Expressed Genes
########

## Load the phenotype data 
print('Loading phenotype data'); flush.console()
samplesCanine <- read.csv(paste(opt$datadir,'samples_canine_updated.csv',sep='/'))

## Read in Count Data
print('Loading count data'); flush.console()
counts <- readDGE(samplesCanine$File, paste(opt$datadir,"dexseq_count/",sep='/'), header = FALSE)
minCount <- round((nrow(samplesCanine)/length(levels(samplesCanine$Hist))), 0)
cpms <- cpm(counts)
noint <- rownames(counts) %in% c("_ambiguous_readpair_position","_ambiguous","_empty","_lowaqual","_notaligned") 
counts_filtered <- counts[rowSums(cpms > 1) >= minCount & !noint,]
samplesCanine$Hist <- relevel(samplesCanine$Hist,"N")
samplesCanine <- samplesCanine %>% mutate(Patient = factor(PatientNumber))

## Estimate Dispersion
print('Estimating dispersion'); flush.console()
design <- model.matrix( ~Patient+Hist, samplesCanine)
d <- DGEList(counts = counts_filtered, group = samplesCanine$Hist)
d <- calcNormFactors(d)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

## Perform DE Analysis
print('Performing differential expresion analysis'); flush.console()
d_fit <- glmFit(d,design)
lrt_list <- list(
    b_n = glmLRT(d_fit,coef = "HistB"),
    m_n = glmLRT(d_fit,coef = "HistM"),
    m_b = glmLRT(d_fit,contrast = makeContrasts(HistM-HistB, levels = d_fit$design)),
    n_b_m = glmLRT(d_fit,coef = c("HistM", "HistB"))
    )
tidy.DGELRT <- function(x, ...) {
  ret <- fix_data_frame(x$table, newcol = "gene")
}
LRTtidied <- ldply(lrt_list, tidy,.id = "contrast") %>% 
  group_by(contrast) %>% 
  mutate(qval = qvalue(PValue)$qvalues)

## Load mappings to human (if option is set)
##    NOTE: Sometimes the bioMart datasets won't load. Just wait a bit and try again, they usually come back up quickly
if(!opt$map.human) {
  print('Using previously generated human mapping file'); flush.console()
  load('humanmapping.rda')
} else {
  print('Generating new human mapping file'); flush.console()
  allgenes <- LRTtidied %>% ungroup %>% select(gene) %>% distinct %>% .$gene

  Map_CanEns2Info <-  getBM(attributes = c('ensembl_gene_id','external_gene_name','description'),
                          filters = list(ensembl_gene_id = allgenes),
                          mart = useMart("ensembl",dataset = "cfamiliaris_gene_ensembl"))

  Map_CanEns2HumEns <- getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),
                           filters = list(ensembl_gene_id = allgenes, with_hsapiens_homolog = TRUE),
                           mart = useMart("ensembl",dataset = "cfamiliaris_gene_ensembl"))

  Map_HumEns2Entrez <- getBM(attributes = c('ensembl_gene_id','entrezgene'),
                           filters = 'ensembl_gene_id',
                           values = Map_CanEns2HumEns$hsapiens_homolog_ensembl_gene,
                           mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl"))

  Map_HumEns2Symb <- getBM(attributes = c('ensembl_gene_id','external_gene_name'),
                         filters = 'ensembl_gene_id',
                         values = Map_CanEns2HumEns$hsapiens_homolog_ensembl_gene,
                         mart = useMart("ensembl",dataset = "hsapiens_gene_ensembl"))

  Map_CanEns2HumEnt <- inner_join(Map_CanEns2HumEns,
                                Map_HumEns2Entrez,
                                by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id"))

  Map_CanEns2HumEnt_unique <- Map_CanEns2HumEnt %>% 
    select(Can_Ens = ensembl_gene_id, Hum_Ent = entrezgene) %>% distinct %>% 
    group_by(Can_Ens) %>% 
    filter(length(Can_Ens) == 1) %>% 
    group_by(Hum_Ent) %>% 
    filter(length(Hum_Ent) == 1) %>%
    ungroup %>% mutate(Hum_Ent = as.character(Hum_Ent))

  Map_CanEns2HumSymb <- inner_join(Map_CanEns2HumEns,
                                 Map_HumEns2Symb,
                                 by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id"))

  Map_CanEns2HumSymb_unique <- Map_CanEns2HumSymb %>% 
    select(Can_Ens = ensembl_gene_id,Hum_Symb = external_gene_name) %>% distinct %>%
    group_by(Can_Ens) %>% 
    filter(length(Can_Ens) == 1) %>% 
    group_by(Hum_Symb) %>% 
    filter(length(Hum_Symb) == 1)

  save(Map_CanEns2HumEns,
       Map_CanEns2Info,
       Map_CanEns2HumEnt,
       Map_CanEns2HumEnt_unique,
       Map_CanEns2HumSymb,
       Map_CanEns2HumSymb_unique,file="humanmapping.rda")

} # End human mapping

## Generate the profile metrics (later used to generate PEPs)
print('Generating profile metrics'); flush.console()
profileMetrics <- LRTtidied %>%
    filter(contrast!="n_b_m") %>%
 select(contrast, gene, logFC, PValue, qval) %>%
   gather(type, value, -contrast, -gene) %>%
   unite(contrast_type, contrast, type, sep=".") %>%
   spread(contrast_type, value) %>%
inner_join(Map_CanEns2HumSymb_unique,by=c("gene"="Can_Ens"))

getCPMWithoutPatient <- function(DGEFit) {
  cpm = cpm(DGEFit$counts,normalized.lib.sizes = TRUE, log = FALSE)
  cpm = cpm + .25 #Add prior
  logcpm = log(cpm)
  patient.cols = grep("Patient", colnames(DGEFit$coefficients))
  due.to.patient = DGEFit$coefficients[, patient.cols] %*% t(DGEFit$design[, patient.cols])
  without.patient = logcpm - due.to.patient
  without.patient = log(exp(without.patient), base = 2) #To make comparable to cpm()
  return(without.patient)
}

## Generate the canine expression data matrices
print('Generating canine data matrices'); flush.console()
CPMmatrices <- list(
  raw = cpm(d_fit$counts, normalized.lib.sizes = FALSE, log = FALSE),
  norm = cpm(d_fit$counts, normalized.lib.sizes = TRUE, log = FALSE),
  znorm = t(scale(t(cpm(d_fit$counts, normalized.lib.sizes = TRUE, log = FALSE)),center = TRUE, scale = TRUE)),
  patregressed = getCPMWithoutPatient(d_fit),
  zpatregressed = t(scale(t(getCPMWithoutPatient(d_fit)),center = TRUE,scale = TRUE))
  )

## Save updated files
print('Saving rda files'); flush.console()
save(CPMmatrices, file=paste(opt$outdir,'CPMmatrices.rda',sep='/'))
save(LRTtidied, file=paste(opt$outdir,'LRTtidied.rda',sep='/'))
save(profileMetrics, file=paste(opt$outdir,'profileMetrics.rda',sep='/'))


##########################
##  Create csv files for expression data (Healthy & Carcinoma only and all samples)
#################

## Save the dog matrices
print('Saving delimited files'); flush.console()
dat.dog <- CPMmatrices$zpatregressed
names.dog <- intersect( Map_CanEns2HumSymb_unique$Can_Ens, rownames(dat.dog) ) # Only care about genes that map to human
Map_CanEns2HumSymb_unique <- Map_CanEns2HumSymb_unique[ Map_CanEns2HumSymb_unique$Can_Ens %in% names.dog, ]
rownames( Map_CanEns2HumSymb_unique ) <- Map_CanEns2HumSymb_unique$Can_Ens

## Map to human, make dog ID's more dog-obvious
dat.dog <- dat.dog[names.dog,]
rownames(dat.dog) <- Map_CanEns2HumSymb_unique$Hum_Symb

## Save to file
write.table(t(dat.dog), file=paste(opt$outdir,'Canine_RNASeq.tab',sep='/'), sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)
write.table(t(dat.dog), file=paste(opt$outdir,'Canine_RNASeq.csv',sep='/'), sep=',', col.names=TRUE, row.names=TRUE, quote=FALSE)

## Load the canine samples data
#samplesCanine <- read.csv('samples_canine_updated.csv')
dog.labs <- data.frame( samplesCanine$Qlabel, samplesCanine$Hist )
colnames(dog.labs) <- c('ID','Hist')
write.table(dog.labs, file=paste(opt$outdir,'Canine_Labels.tab',sep='/'), sep='\t', col.names=TRUE, row.names=FALSE, quote=FALSE)

## Now grab malignant vs normal - remove benign samples (for tumor-normal analyses)
x <- t(dat.dog[ ,as.character(dog.labs[ dog.labs$Hist !='B', 'ID' ]),drop=FALSE])
write.table(x, file=paste(opt$outdir,'Canine_RNASeq_TumorVSNormal.t.tab',sep='/'), sep='\t', col.names=TRUE, row.names=TRUE, quote=FALSE)

print('Finished. Success!'); flush.console()
