## Kiley Graim
## March 2018
## 
## Function to create PEP lists from a dataset and sample list

## Generates the PEPs


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


## Load the previously prepared human mapping file
load('humanmapping.rda')

############
### Prepare Data & Identify Differentially Expressed Genes
########
library(devtools)
knitr::opts_chunk$set(warning = FALSE, message = FALSE, prompt = FALSE)

## Read in Count Data
#samplesCanine <- read.csv('samples_canine_updated.csv')
samplesCanine <- read.csv(opt$samplesCanine)
counts <- readDGE(samplesCanine$File, paste0(opt$countdir,'/'), header = FALSE)
minCount <- round((nrow(samplesCanine)/length(levels(samplesCanine$Hist))), 0)
cpms <- cpm(counts)
noint <- rownames(counts) %in% c("_ambiguous_readpair_position","_ambiguous","_empty","_lowaqual","_notaligned") 
counts_filtered <- counts[rowSums(cpms > 1) >= minCount & !noint,]
samplesCanine$Hist <- relevel(samplesCanine$Hist,"N")
samplesCanine <- samplesCanine %>% mutate(Patient = factor(PatientNumber))
## TODO: Add error handling - make sure we have N,B,M Histology levels

## Estimate Dispersion
design <- model.matrix( ~Patient+Hist, samplesCanine)
d <- DGEList(counts = counts_filtered, group = samplesCanine$Hist)
d <- calcNormFactors(d)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

## TODO: For some reason next step this crashes if you don't first run it in the terminal... why???
## Perform DE Analysis
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

## Create the stats
profileMetrics <- LRTtidied %>%
    filter(contrast!="n_b_m") %>%
  select(contrast, gene, logFC, PValue, qval) %>%
    gather(type, value, -contrast, -gene) %>%
    unite(contrast_type, contrast, type, sep=".") %>%
    spread(contrast_type, value) %>% 
  inner_join(Map_CanEns2HumSymb_unique,by=c("gene"="Can_Ens")) 

qval.lim <- 0.05 # Set this to be the cutoff for q-value

## Now generate profiles for each PEP based on their definitions
## Version 0 - defined based on what the paper methods section says we do
# Tumor-Specific : A gene differentially expressed in both the normal-adenoma and the normal-carcinoma comparisons; the sign of the change is the same for both comparisons.
tumor.sp  <- profileMetrics[ profileMetrics$b_n.qval < qval.lim & profileMetrics$m_n.qval < qval.lim & (sign(profileMetrics$b_n.logFC) == sign(profileMetrics$m_n.logFC)),]
# Carcinoma-Specific : A gene differentially expressed in both the normal-carcinoma and the adenoma-carcinoma comparisons; the sign of the change is the same for both comparisons.
malign.sp <- profileMetrics[ profileMetrics$m_n.qval < qval.lim & profileMetrics$m_b.qval < qval.lim & (sign(profileMetrics$m_n.logFC) == sign(profileMetrics$m_b.logFC)),]
# Adenoma-Specific : A gene differentially expressed in both the normal-adenoma and the adenoma-carcinoma comparisons; the sign of the change is different between the two comparisons.
benign.sp <- profileMetrics[ profileMetrics$b_n.qval < qval.lim & profileMetrics$m_b.qval < qval.lim & (sign(profileMetrics$b_n.logFC) != sign(profileMetrics$m_n.logFC)),]

## Using the PEP profiles, create the PEPs
peps <- data.frame(EnsGene=profileMetrics$gene, HumanSymbol=profileMetrics$Hum_Symb, Tumor_Expression_Pattern=1, Carcinoma_Expression_Pattern=1, Adenoma_Expression_Pattern=1)
peps$Tumor_Expression_Pattern[peps$EnsGene %in% tumor.sp$gene] <- min(tumor.sp$b_n.qval, tumor.sp$m_n.qval)
peps$Adenoma_Expression_Pattern[peps$EnsGene %in% benign.sp$gene] <- min(benign.sp$b_n.qval, benign.sp$m_b.qval)
peps$Carcinoma_Expression_Pattern[peps$EnsGene %in% malign.sp$gene] <- min(malign.sp$m_n.qval, malign.sp$m_b.qval)

gen_PEPs <- function(samplesCanine, fn.cDir='dexseq_count') {
  ############
  ### Prepare Data & Identify Differentially Expressed Genes
  ########

  ## Read in Count Data
  counts <- readDGE(samplesCanine$File, fn.cDir, header = FALSE) ## TODO: Make folder an optional argument 
  minCount <- round((nrow(samplesCanine)/length(levels(samplesCanine$Hist))), 0)
  cpms <- cpm(counts)
  noint <- rownames(counts) %in% c("_ambiguous_readpair_position","_ambiguous","_empty","_lowaqual","_notaligned") 
  counts_filtered <- counts[rowSums(cpms > 1) >= minCount & !noint,]
  samplesCanine$Hist <- relevel(samplesCanine$Hist,"N")
  samplesCanine <- samplesCanine %>% mutate(Patient = factor(PatientNumber))

  ## Make sure we have N,B,M Histology levels, quit otherwise TODO: make this return with an error, not quit entirely
  if(!all(c('B','N','M') %in% levels(samplesCanine$Hist))) {
    print('ERROR: Must have all histologies: N,B,M'); quit(save='no',status=1)
  }

  ## Estimate Dispersion
  design <- model.matrix( ~Patient+Hist, samplesCanine)
  d <- DGEList(counts = counts_filtered, group = samplesCanine$Hist)
  d <- calcNormFactors(d)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)

## TODO: For some reason next step this crashes if you don't first run it in the terminal... why???
  ## Perform DE Analysis
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

  ## Create the stats
  profileMetrics <- LRTtidied %>%
      filter(contrast!="n_b_m") %>%
    select(contrast, gene, logFC, PValue, qval) %>%
      gather(type, value, -contrast, -gene) %>%
      unite(contrast_type, contrast, type, sep=".") %>%
      spread(contrast_type, value) %>% 
    inner_join(Map_CanEns2HumSymb_unique,by=c("gene"="Can_Ens")) 

  ## Set the q-value cutoff for PEP profiles
  qval.lim <- 0.05 

  ## Now generate profiles for each PEP based on their definitions
  ## Version 0 - defined based on what the paper methods section says we do
    # Tumor-Specific : A gene differentially expressed in both the normal-adenoma and the normal-carcinoma comparisons; the sign of the change is the same for both comparisons.
  tumor.sp  <- profileMetrics[ profileMetrics$b_n.qval < qval.lim & profileMetrics$m_n.qval < qval.lim & (sign(profileMetrics$b_n.logFC) == sign(profileMetrics$m_n.logFC)),]
    # Carcinoma-Specific : A gene differentially expressed in both the normal-carcinoma and the adenoma-carcinoma comparisons; the sign of the change is the same for both comparisons.
  malign.sp <- profileMetrics[ profileMetrics$m_n.qval < qval.lim & profileMetrics$m_b.qval < qval.lim & (sign(profileMetrics$m_n.logFC) == sign(profileMetrics$m_b.logFC)),]
    # Adenoma-Specific : A gene differentially expressed in both the normal-adenoma and the adenoma-carcinoma comparisons; the sign of the change is different between the two comparisons.
  benign.sp <- profileMetrics[ profileMetrics$b_n.qval < qval.lim & profileMetrics$m_b.qval < qval.lim & (sign(profileMetrics$b_n.logFC) != sign(profileMetrics$m_n.logFC)),]

  ## Using the PEP profiles, create the PEPs
  peps <- data.frame(EnsGene=profileMetrics$gene, HumanSymbol=profileMetrics$Hum_Symb, Tumor_Expression_Pattern=1, Carcinoma_Expression_Pattern=1, Adenoma_Expression_Pattern=1)
  peps$Tumor_Expression_Pattern[peps$EnsGene %in% tumor.sp$gene] <- min(tumor.sp$b_n.qval, tumor.sp$m_n.qval)
  peps$Adenoma_Expression_Pattern[peps$EnsGene %in% benign.sp$gene] <- min(benign.sp$b_n.qval, benign.sp$m_b.qval)
  peps$Carcinoma_Expression_Pattern[peps$EnsGene %in% malign.sp$gene] <- min(malign.sp$m_n.qval, malign.sp$m_b.qval)

  ## Return PEPs 
  return(peps)
} # end gen_PEPs
