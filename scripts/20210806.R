library(readxl)
library(dplyr)
library(tidyverse) 

#Import top 500 sex diff expressed genes in each tissue as diff data frames. Data source: Meritxell Oliva, science, 2020. 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- excel_sheets(filename)
  x<- lapply(sheets,function(X) read_excel(filename, sheet = X))
  names(x) <- sheets
  if(!tibble) x <- lapply(x, as.data.frame)
  x
}
mysheets <- read_excel_allsheets("data/top 500 diff ex genes in each tissue by science.xlsx")

#Delete README and Replication in mysheets
mysheets <- mysheets[-c(1,2)]

#Extract microRNA from top 500 genes in each tissue
mysheets <- lapply(mysheets, function(X){
  X<- 
    X %>% dplyr::mutate(ENSEMBL_gene_id = gsub("\\..*","",ENSEMBL_gene_id))
})

#load all the microRNAs 
load("data/genemap_micro.RData")
mysheets2 <- lapply(mysheets, function(X) {
  X <- merge(X, genemap_micro, by.x = "ENSEMBL_gene_id", by.y = "ensembl_gene_id", all =FALSE)
}
  )
#------------------------------------------------No microRNA is in top 500 sex bias genes------------------------------------------