library(reshape2)
library(plyr)
library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggthemes)

#Import top 500 sex diff expressed genes in each tissue as diff data frames. Data source: Meritxell Oliva, science, 2020. 
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename) #sheets contain all the names of sheet
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X)) #list apply:laaply(d, FUN)
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

mysheets <- read_excel_allsheets() #A total of 44 tissues

#Delete the README and Replication in mysheets
mysheets <- mysheets[-c(1,2)]

#Extract lncRNA from top 500 genes in each tissue
mysheets <- lapply(mysheets, function(X) {
  X <- 
    X %>% 
    mutate(ENSEMBL_gene_id = gsub("\\..*","", ENSEMBL_gene_id))
  })

  ##load all the lncs 
  load("data/genemap_lnc.RData")

mysheets1 <- lapply(mysheets, function(X) {
  X <- merge(X, genemap_lnc, by.x = "ENSEMBL_gene_id", by.y = "ensembl_gene_id", all = FALSE)
}
  )

#Plot the number of lncs in each tissue into one bar plot
number <- sapply(mysheets1, function (X) {
   NROW(X)
}
  )
number <- as.data.frame(number)
number$tissue <- rownames(number)
number <- number %>%  arrange(desc(number))
  ggplot(data = number, aes(x = reorder(tissue, number), y = number)) +
  ylab("Number of Sex Bias lncRNA") + xlab("Tissue")+ ylim(0,85)+
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label=number), hjust = -0.3, size = 3)+
  coord_flip()

#Total number of sex bias lncRNA that is significant in at least one tissue
f<- as.data.frame(unlist(sapply(mysheets1, function(X) {
  X$ENSEMBL_gene_id})))
  total_num <- unique(f)
  
#Extract lncRNAs that are significant in at least 20 tissues)
df <- ldply(mysheets1, data.frame)

output <- data.frame(NULL)
for (i in unique(df$ENSEMBL_gene_id)){
  temp <- subset(df, ENSEMBL_gene_id == i)
  if (nrow(temp) >= 20){
    output <- rbind(output, temp)
  }
}

#Check how many lncRNAs were selected in the previous step
output %>% distinct(ENSEMBL_gene_id) %>% nrow()

#Plot heatmap using significance/effect size of each lncRNA 
beta<- 
  ggplot(data = output[output$ENSEMBL_gene_id != "ENSG00000229807",], aes(x=.id, y=ENSEMBL_gene_id, fill=MASH.beta)) +
#ggplot(data = output, aes(x=.id, y=ENSEMBL_gene_id, fill=MASH.beta)) +   #exclude XIST
  geom_tile() +
  theme_bw() +
  theme_tufte()+
  theme(axis.text.x = element_text(angle = 90), aspect.ratio = 1.05/1, axis.title = element_blank(), legend.position = "left",
        legend.title = element_text(size = 9)) +
  scale_fill_gradient2(midpoint=0, low="#B2182B", high="#2166AC") + 
  guides(fill = guide_colorbar(title = "MASH \n beta", ticks = FALSE, label = FALSE))

output <- 
  output %>%
  mutate(LFSR_new = if_else(MASH.LFSR <= 10^(-40), -log10(10^(-40)), -log10(MASH.LFSR)))

LFSR<- 
  #ggplot(data=output, aes(x=.id, y=ENSEMBL_gene_id, fill=LFSR_new)) +   
  ggplot(data = output[output$ENSEMBL_gene_id != "ENSG00000229807",], aes(x=.id, y=ENSEMBL_gene_id, fill=LFSR_new)) + #exclude XIST
  geom_tile() +
  theme_bw() +
  theme_tufte()+
  scale_fill_gradient(high = "darkseagreen4", low = "white")+
  theme(axis.text.x = element_text(angle = 90), aspect.ratio = 1.05/1, axis.title = element_blank(), axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), legend.title = element_text(size = 9)) +
  #scale_fill_gradient2() +
  guides(fill = guide_colorbar(title = "-log10\n(MASH.LFSR)", ticks = FALSE))

beta +LFSR
