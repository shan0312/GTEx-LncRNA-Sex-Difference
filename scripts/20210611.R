library(readr)
library(tidyverse)
if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("biomaRt")
library(biomaRt)
library(dplyr)
 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 
 BiocManager::install("rtracklayer")
library(rtracklayer)

#-------------------------------------------------------------oringinal GTEx data-------------------------------------------------------------
#GTEx tissues to process
tissulist <- read.delim('data/tissuelist.txt', stringsAsFactors = FALSE, header = FALSE)[,1]

#Download the GTEx read count matrix for all genes and all sampels (copied from Dr. Stranger's Science Paper)
counts <- as.data.frame(read_delim('data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',col_names=T, 
              comment="!", skip=2, delim='\t', col_types = cols(.default = col_double(), Name = col_character(),
              Description = col_character())))

#Figure out the function of each command above
#as.data.frame is not required for this input data file
#counts1 <- read_delim('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',
                      #col_names=T, #not required, In read_delim, default col_names=T
                      #comment="!", #not required for this file
                      #skip=2, #required, skip the top two rows
                      #delim='\t' #required, specify the delimiter type
                      #col_types = cols(.default = col_double(), Name = col_character(), Description = col_character())) #not required, default 
#----------------------------------------------------------------------------------------------------------------------------------------------

#Set ENSEMBL ID as rownames, and delete gene names
rownames(counts) <- counts[,1]
counts <- counts[,-c(1,2)]
#-------------------------------------------------------------oringinal GTEx data-------------------------------------------------------------

#----------------------------------------------extract genes by protein coding and non-protein coding------------------------------------------
#METHOD 1: Download GTEx_Liver_counts from Huang lab drive uploaded by Yingbo (package: biomaRt) [a total of 13704 lncRNAs were defined]
GTEx_Liver_counts <- as.data.frame(read_delim('data/GTEx_Liver_counts'))
colnames(GTEx_Liver_counts)[1] <- "Ensembl_ID"
GTEx_Liver_counts$Ensembl_ID <- gsub("\\..*","",GTEx_Liver_counts$Ensembl_ID)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genemap <- getBM(attributes= c("ensembl_gene_id", "gene_biotype", "chromosome_name","external_gene_name",
                               "ensembl_gene_id_version"), mart=ensembl) 
dim(genemap)
  #extract lncRNA in genemap
  genemap_lnc <- genemap[genemap$gene_biotype == "lncRNA",]
  #extract microRNA in genemap
  genemap_micro <- genemap[genemap$gene_biotype == "miRNA",]
  
# For all the genes screened in the GTEx, add gene_biotype from the corresponding genemap, extract lncRNA
GTEx_Liver_counts <- merge(GTEx_Liver_counts, genemap, by.x = "Ensembl_ID", by.y="ensembl_gene_id", all = FALSE)
GTEx_Liver_counts_lnc <- GTEx_Liver_counts[GTEx_Liver_counts$gene_biotype == "lncRNA",]
GTEx_lnc <- GTEx_Liver_counts_lnc %>% 
  dplyr::select(Ensembl_ID, gene_biotype, external_gene_name)

#METHOD 2: using the gencode annotation used in Olivia M. de Goede cell paper (2021) (package: rtracklayer) [a total of 7520 lncRNAs were defined]
gtf <- rtracklayer::import('data/gencode.v26.annotation.gtf')
gtf_df=as.data.frame(gtf)
linc <- gtf_df[gtf_df$gene_type == "lincRNA",]
#--------------------------------------------------------Method 1 was used for further analysis----------------------------------------