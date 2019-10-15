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
# places all info into encori_data, which is a list of dataframes for each gene
web2<-'http://genome.ucsc.edu/cgi-bin/das/hg18/dna?segment='
for(d in encori_data){
  df <- data.frame(d)
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
}

# get mRNA sequences (ALL isoforms) for each gene #######
setwd("~/Desktop/sirna_pred/LocalSequencePredictions/data/")

# get accession numbers rna.gbk_gene2acc.txt

accessions <- c()
for(g in all_genes){
  l <- grep(paste('^',g,'\t',sep=""),readLines("rna.gbk_gene2acc.txt"))
  accessions<- c(accessions, lapply(readLines("rna.gbk_gene2acc.txt")[l],gsub,pattern=paste(g,"\t",sep=""),replacement=""))
}
# get mRNA sequences formatted_rna.fa.concatenated_new.txt
## TODO: change so actually gets all accessions (not just first 3)
acc2 <- as.character(accessions[1:3])
mrna_sequences <- c()
gene_names <- c()
for(a in acc2){
  l <- grep(paste('^',a,'\t',sep=""),readLines("formatted_rna.fa.concatenated_new.txt"))
  seq<-readLines("formatted_rna.fa.concatenated_new.txt")[l]
  pos_seq_start<-sapply(str_locate_all(seq, '\t'), tail, 1)[1]
  pos_gene_start<-sapply(str_locate_all(seq, '\t'), head, 1)[1]
  sq <- gsub("[\t]", "",substr(seq, pos_seq_start,nchar(seq)))
  gene <-gsub("[\t]", "",substr(seq, pos_gene_start, pos_seq_start))
  mrna_sequences <- c(mrna_sequences, seq)
  gene_names <- c(gene_names, gene)
}
# acc2 <- as.list(acc2)
# mrna_sequences <- as.list(mrna_sequences)
mrna_dat<-do.call(rbind, Map(data.frame, gene = gene_names, accessions=acc2, seq=mrna_sequences,  stringsAsFactors = FALSE))
View(mrna_dat)
# TODO: only look at regions in 3' UTR (need rna.gbk.cds.txt)

# for each gene, compare rbp targeting sequence with the mRNA sequenence ###########
comp <- mrna_dat[1,]


# get appropriate encori data for the gene name
#    encori data has all genes listed alphabetically
all_genes<-sort(all_genes) # sort all genes
# check that length of all_genes matches encori data length (so not missing any data)
length(all_genes) == length(encori_data)
# get index of gene in all_genes (will match that in encori)
i <- match(comp$gene,all_genes)
c_encori_df <- data.frame(encori_data[6])
View(c_encori_df)
cur_df$geneName[1] == comp$gene



# compile into single file