install.packages("readxl")
install.packages("cowplot")
install.packages("tidyverse")
install.packages("RCurl")
library("readxl")
library(cowplot)
library("ggplot2")
library(tidyverse)
library(gridExtra)
library(RCurl)


# get list of siRNA human gene targets
# from:
#TODO:
#   luciferase
setwd("~/Google_Drive/research/sirna_rbp/data/")
# BDNA  ########
bdna_data <- read_excel("sirna_data/bdna.xlsx")
# make new column with gene name info
names1 <- (sapply(strsplit(bdna_data$`Duplex Name`,"_"), `[`, 1))
names1 <- (sapply(strsplit(names1," "), `[`, 1)) # remove HTT intron/exon (since same gene)
names1[(grepl("mSTAT3", names1))] <- "mSTAT3" #change mSTAT3 duplicates
names1
bdna_data$gene_names <-names1
#get human only
h_bdna_data <- bdna_data[substring(bdna_data$gene_names,1,1)=="h",]
nrow(h_bdna_data)
h_bdna_genes <- unique(substring(h_bdna_data$gene_names,2))

# qPCR ########
qpcr_fm_data <- read_excel("sirna_data/qPCR_fm.xlsx")
# remove NAs
r <- seq(from=524,to=nrow(qpcr_fm_data)+1,by=1)
qpcr_fm_data<-qpcr_fm_data[-r, ]

# make new column with gene name info
names2 <- (sapply(strsplit(qpcr_fm_data$`Duplex Name`,"_"), `[`, 1))
names2 <- (sapply(strsplit(names2," "), `[`, 1)) # remove HTT intron/exon (since same gene)
qpcr_fm_data$gene_names<-names2
l <- lengths(split(bdna_data$gene_names, bdna_data$gene_names))
col <- rep(seq_along(l), l)

#get human only

h_qpcr_fm_data <- qpcr_fm_data[toupper(qpcr_fm_data$gene_names) == qpcr_fm_data$gene_names,]
nrow(h_qpcr_fm_data)
h_qpcr_fm_genes <- unique(h_qpcr_fm_data$gene_names)

# combine genes and output encori fetch script #########
all_genes = unique(c(h_bdna_genes,h_qpcr_fm_genes))
# number of genes:
length(all_genes)

# get RBP binding data for each gene
# make file with command:
web <- 'http://starbase.sysu.edu.cn/api/RBPTarget/?assembly=hg19&geneType=mRNA&RBP=all&clipExpNum=5&pancancerNum=0&target='
cell <- '&cellType=all'
out_file_name <-" > ENCORI_hg19_RBPTarget_"

encori_cmds <- c()
for(g in all_genes){
  f <-  paste("curl ","'",web,g,cell,"'",out_file_name,g,".csv", sep = "")
  encori_cmds <- c(encori_cmds,f)
}

write.table(encori_cmds, file="fetch_encori.sh", append = FALSE, sep = " ", dec = ".",
            row.names = FALSE, col.names = FALSE,quote=FALSE)
# Read in encori RBP data ############
# TODO: MANUALLY run the bash script (fetch_encori.sh) from the command line

setwd("~/Google_Drive/research/sirna_rbp/data/encori_data/")
temp = list.files(pattern="*.csv")
encori_data = lapply(temp, read.csv, header=TRUE, sep="\t", skip = 2, comment.char = "#")
setwd("~/Google_Drive/research/sirna_rbp/data/")
# get sequences that correspond to chromosome coordinates ##########
web2<-'http://genome.ucsc.edu/cgi-bin/das/hg18/dna?segment='

# for(d in encori_data){
  # df <- data.frame(d)
  for(r in 1:nrow(df)){
    chr<- df$chromosome[r]
    start<- df$broadStart[r]
    stop<-df$broadEnd[r]
    out<-getURL(paste(web2,chr,":",start,",",stop, sep = ""))
    # get sequence
    f<- sapply(str_locate_all(out, '\">\n'), tail, 1)[2]
    l<-sapply(str_locate_all(out, '\n</DNA>\n</SEQUENCE>\n</DASDNA>\n'),head,1)[1]
    seq <- gsub("[\n]", "",substr(out, f,l)) 
    # add sequence to dataframe
    df$rbp_target_sequence[r] = seq
  }
# }

# get mRNA sequences (ALL isoforms) for each gene
setwd("~/Desktop/sirna_pred/LocalSequencePredictions/data/")
# test with ACAD8
# get accession numbers rna.gbk_gene2acc.txt


# get mRNA sequences formatted_rna.fa.concatenated_new.txt
# TODO: only look at regions in 3' UTR (need rna.gbk.cds.txt)

# compile into single file