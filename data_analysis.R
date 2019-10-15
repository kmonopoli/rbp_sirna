install.packages("readxl")
install.packages("cowplot")
install.packages("tidyverse")
library("readxl")
library(cowplot)
library("ggplot2")
library(tidyverse)
library(gridExtra)

setwd("~/Dropbox (UMass Medical School)/research/sirna_rbp/data/")
bdna_data <- read_excel("bdna.xlsx")

# make new column with gene name info
names1 <- (sapply(strsplit(bdna_data$`Duplex Name`,"_"), `[`, 1))
names1 <- (sapply(strsplit(names1," "), `[`, 1)) # remove HTT intron/exon (since same gene)

names1[(grepl("mSTAT3", names1))] <- "mSTAT3" #change mSTAT3 duplicates
names1
bdna_data$gene_names <-names1

l <- lengths(split(bdna_data$gene_names, bdna_data$gene_names))
col <- rep(seq_along(l), l)

# barplot(bdna_data$`mRNA expression` ,col = col,border=NA,legend.text =unique(names) )

q1 <- ggplot(bdna_data)+ 
  geom_col(mapping = aes(x = bdna_data$`Duplex Name`,
                         # y = bdna_data$`mRNA expression`,
                         y = bdna_data$`Normalized mRNA Expression`,
                         
                         fill=bdna_data$gene_names))+
  
  theme_cowplot(12)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title= element_blank())+
  labs(y = "target expression (% of control)",title = "BDNA")
  

# Create summary table
SummaryTable <- data.frame(
  total_siRNAs = nrow(bdna_data),
  number_of_genes = length(unique(names1)),
  # max_target_expression= max(bdna_data$`mRNA expression`),
  # min_target_expression= min(bdna_data$`mRNA expression`),
  max_target_expression= max(bdna_data$`Normalized mRNA Expression`),
  min_target_expression= min(bdna_data$`Normalized mRNA Expression`),
  avg_number_of_sequences_per_gene = round(length(names1)/length(unique(names1)),digits=0)
)
mytheme <- ttheme_default(base_size = 8, colhead=list(fg_params = list(parse=TRUE)), padding = unit(c(5, 5), "mm") )
tbl_nms = c("total siRNAs","num of genes","max target expression","min target expression","number of sequences per gene ~")
tbl1 <- tableGrob(SummaryTable, cols=tbl_nms,rows=NULL, theme = mytheme)

######## qPCR Data ########

qpcr_fm_data <- read_excel("old/qPCR_fm.xlsx")
# remove NAs
r <- seq(from=524,to=nrow(qpcr_fm_data)+1,by=1)
qpcr_fm_data<-qpcr_fm_data[-r, ]

# make new column with gene name info
names2 <- (sapply(strsplit(qpcr_fm_data$`Duplex Name`,"_"), `[`, 1))
names2 <- (sapply(strsplit(names2," "), `[`, 1)) # remove HTT intron/exon (since same gene)



qpcr_fm_data$gene_names<-names2
l <- lengths(split(bdna_data$gene_names, bdna_data$gene_names))
col <- rep(seq_along(l), l)

# barplot(bdna_data$`mRNA expression` ,col = col,border=NA,legend.text =unique(names) )
dat_cutoff <-(qpcr_fm_data[qpcr_fm_data$Expression...15 < 1.5,])

q2 <- ggplot(dat_cutoff)+ 
  geom_col(mapping = aes(x = dat_cutoff$`Duplex Name`,y = dat_cutoff$Expression...15,fill=dat_cutoff$gene_names))+
  theme_cowplot(12)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title= element_blank())+
  labs(x = "siRNAs",y = "target expression (% of control)",title = "qPCR fully-modified")+
  theme(legend.position = "none")+
  ylim(0,1.5)


q2norm <- ggplot(qpcr_fm_data)+ 
  geom_col(mapping = aes(x = qpcr_fm_data$`Duplex Name`,y = qpcr_fm_data$`Normalized Per Exp`,fill=qpcr_fm_data$gene_names))+
  theme_cowplot(12)+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title= element_blank())+
  labs(y = "target expression (% of control)",title = "qPCR fully-modified")


# Create summary table
SummaryTable <- data.frame(
  total_siRNAs = nrow(qpcr_fm_data),
  number_of_genes = length(unique(names2)),
  max_target_expression= max(qpcr_fm_data$Expression...15),
  min_target_expression= min(qpcr_fm_data$Expression...15),
  avg_number_of_sequences_per_gene = round(length(names2)/length(unique(names2)),digits=0)
)
mytheme <- ttheme_default(base_size = 8, colhead=list(fg_params = list(parse=TRUE)), padding = unit(c(5, 5), "mm") )
tbl_nms = c("total siRNAs","num of genes","max target expression","min target expression","number of sequences per gene ~")
tbl2 <- tableGrob(SummaryTable, cols=tbl_nms,rows=NULL, theme = mytheme)




lay <- rbind(c(1,1,1),#3,3,3),
             c(1,1,1),#3,3,3),
             c(2,2,2))#4,4,4))

grid.arrange(q2, 
             tbl2,
             layout_matrix=lay)

grid.arrange(q1, 
             tbl1,
             layout_matrix=lay)


View(bdna_data)




