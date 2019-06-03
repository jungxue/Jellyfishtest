#-------------------------------------Premeables--------------------------------------------------#

set.seed(12345678)

setwd("~/Desktop/Tue")                                                                  # Home
# setwd("H:/Desktop/7.1 No time to waste/Stats 798A Research Master/R-codes/BoxJellyfish") # Uni


# Install Bioconductor (for genetics packages)----------#
# https://www.bioconductor.org/

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

## Installing ------------------------------------------#

### Genetics packages

# goTools         # GO functions
# GOpro           # Find the most characteristic gene ontology terms for groups of human genes
# GOexpress       # Visualise microarray and RNAseq data using gene ontology annotations
# biomaRt         # BioMart databases 
# clusterProfiler # statistical analysis and visualization of functional profiles for genes and gene clusters
# Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609
# org.Hs.eg.db    # gene data base
# GO.db           # gene data base
# phenoTest       # Tools to test association between gene expression and phenotype

BioPackages <- c("GOpro","GOexpress", "goTools","org.Hs.eg.db","biomaRt","clusterProfiler","GO.db",
                 "GOFunction","phenoTest")
BiocManager::install(BioPackages)



### Other packages

Packages <- c("glue","readxl", "tidyverse","xtable","rjags","ggplot2")
lapply(Packages, install.packages, character.only = TRUE)

## Loading --------------------------------------------#

Packages2 <- c(BioPackages, Packages)

lapply(Packages2, library, character.only = TRUE)


#----------------------------------Data Cleanning-----------------------------------------------------#

# Import data
### https://www.nature.com/articles/s41467-019-09681-1#Sec32

genedata <- read_excel("41467_2019_9681_MOESM2_ESM.xlsx")
#genedata <- genedata[complete.cases(genedata), ] # get rid of NA in Fold change
genedata[1:20,]
genedata <- genedata[order(genedata$Gene), ] #Sort by genename
genedata[1:20,]

dim(genedata)
summary(genedata)

nRNA  = nrow(genedata)                       # number of sgRNAs
nnRNA = unname(unlist(table(genedata$Gene))) # number of sgRNAs per gene
genenames = names(table(genedata$Gene))      # name of genes
nGene = length(genenames)                    # number of genes
location  = cumsum(nnRNA)-5                  # location of first obs of sgRNA of certain gene

### Find avg of log Fold change ---------------------#

table(nRNA)

genedata_sum <- genedata %>% 
  group_by(Gene) %>% 
  summarise(avg_fold_changes = mean(log(`Fold changes`))) #Use mean log value of each gene   

genedata_sum <- genedata_sum[order(genedata_sum$avg_fold_changes,decreasing=T),] #Order by change decs

write.csv(genedata_sum, 'genedata_sum.csv') 
genedata_sum = read.csv("genedata_sum.csv")
genedata_sum = genedata_sum [,-1]                                           # delete index
genedata_sum <- genedata_sum[-which(is.na(genedata_sum$avg_fold_changes)),] # delete NA

dim(genedata_sum )
head(genedata_sum)

### plot density
pdf("GeneDensity.pdf") 

plot = ggplot(genedata_sum, aes(x=avg_fold_changes)) + 
  geom_histogram(aes(y=..density..),binwidth=.2,colour="dark grey", fill="light blue") +
  geom_density(alpha=.2, fill="light blue")  
print(plot + ggtitle("Density of mean log fold change for all Genes")+
        labs(y="Density", x = "mean of log fold change"))

plot(density(genedata_sum$avg_fold_changes), main = "Density of mean log fold change for all Genes")
dev.off()
####-----------------#####

table(nRNA)

genedata_sum <- genedata %>% 
  group_by(Gene) %>% 
  summarise(avg_fold_changes = mean(`Fold changes`)) #Use mean log value of each gene   

genedata_sum <- genedata_sum[order(genedata_sum$avg_fold_changes,decreasing=T),] #Order by change decs

write.csv(genedata_sum, 'genedata_sum.csv') 
genedata_sum = read.csv("genedata_sum.csv")
genedata_sum = genedata_sum [,-1]                                           # delete index
genedata_sum <- genedata_sum[-which(is.na(genedata_sum$avg_fold_changes)),] # delete NA

dim(genedata_sum )
head(genedata_sum)




# source("bayes-t-test.R")

############################################################
#    avg of Fold change
############################################################


head(genedata_sum)
summary(genedata_sum)

### Use 100% sample

num = 1:nrow(genedata_sum)
smp1 = sample(num,round(nrow(genedata_sum)*1,0))
sample1 = genedata_sum[smp1,]
summary(sample1)

### Create levels

lv1 = groupGO(as.character(sample1$Gene), keyType = "SYMBOL", ont = "BP", level = 1,OrgDb = org.Hs.eg.db, readable = FALSE) # assign levels
lv2 = groupGO(as.character(sample1$Gene), keyType = "SYMBOL", ont = "BP", level = 2,OrgDb = org.Hs.eg.db, readable = FALSE)
lv3 = groupGO(as.character(sample1$Gene), keyType = "SYMBOL", ont = "BP", level = 3,OrgDb = org.Hs.eg.db, readable = FALSE)
lv4 = groupGO(as.character(sample1$Gene), keyType = "SYMBOL", ont = "BP", level = 4,OrgDb = org.Hs.eg.db, readable = FALSE)
lv5 = groupGO(as.character(sample1$Gene), keyType = "SYMBOL", ont = "BP", level = 5,OrgDb = org.Hs.eg.db, readable = FALSE)

lv1_members = as.data.frame(summary(lv1)) # output summary as dataframe
lv2_members = as.data.frame(summary(lv2))
lv3_members = as.data.frame(summary(lv3))
lv4_members = as.data.frame(summary(lv4))
lv5_members = as.data.frame(summary(lv5))

lv1_members = lv1_members[order(lv1_members$Count,decreasing = T),] # order highest count 
lv2_members = lv2_members[order(lv2_members$Count,decreasing = T),]
lv3_members = lv3_members[order(lv3_members$Count,decreasing = T),]
lv4_members = lv4_members[order(lv4_members$Count,decreasing = T),]
lv5_members = lv5_members[order(lv5_members$Count,decreasing = T),]

lv1_members = lv1_members[!lv1_members$Count == 0,] # delete all 0 counts
lv2_members = lv2_members[!lv2_members$Count == 0,]
lv3_members = lv3_members[!lv3_members$Count == 0,]
lv4_members = lv4_members[!lv4_members$Count == 0,]
lv5_members = lv5_members[!lv5_members$Count == 0,]

lv1_members$Level = 1 # assign levels
lv2_members$Level = 2
lv3_members$Level = 3
lv4_members$Level = 4
lv5_members$Level = 5


lv_nrow = c(nrow(lv1_members),
            nrow(lv2_members),
            nrow(lv3_members),
            nrow(lv4_members),
            nrow(lv5_members))
lv_nrow 

### for full data ###########################

# Level 1  2   3    4     5 
#       1 32 558 3890 11017   # all

# Level 1  2   3    4     5 
#       1 29 383 2304 6100   # no 0 counts

#############################################

# Assign sum of avg folds for all Genes

### lists containing all Genes of group for each level

lv1_gene.list = list()
lv2_gene.list = list()
lv3_gene.list = list()
lv4_gene.list = list()
lv5_gene.list = list()

for (i in 1:lv_nrow[1]){lv1_gene.list [[i]]<-unlist(strsplit(as.character(lv1_members$geneID[i]), "/"))}
for (i in 1:lv_nrow[2]){lv2_gene.list [[i]]<-unlist(strsplit(as.character(lv2_members$geneID[i]), "/"))}
for (i in 1:lv_nrow[3]){lv3_gene.list [[i]]<-unlist(strsplit(as.character(lv3_members$geneID[i]), "/"))}
for (i in 1:lv_nrow[4]){lv4_gene.list [[i]]<-unlist(strsplit(as.character(lv4_members$geneID[i]), "/"))}
for (i in 1:lv_nrow[5]){lv5_gene.list [[i]]<-unlist(strsplit(as.character(lv5_members$geneID[i]), "/"))}

lv1_gene_loc.list = list()
lv2_gene_loc.list = list()
lv3_gene_loc.list = list()
lv4_gene_loc.list = list()
lv5_gene_loc.list = list()

for (i in 1:lv_nrow[1]){lv1_gene_loc.list [[i]]<-which(rowSums(outer(genedata_sum$Gene,lv1_gene.list[[i]],"=="))==1)}
for (i in 1:lv_nrow[2]){lv2_gene_loc.list [[i]]<-which(rowSums(outer(genedata_sum$Gene,lv2_gene.list[[i]],"=="))==1)}
for (i in 1:lv_nrow[3]){lv3_gene_loc.list [[i]]<-which(rowSums(outer(genedata_sum$Gene,lv3_gene.list[[i]],"=="))==1)}
for (i in 1:lv_nrow[4]){lv4_gene_loc.list [[i]]<-which(rowSums(outer(genedata_sum$Gene,lv4_gene.list[[i]],"=="))==1)}
for (i in 1:lv_nrow[5]){lv5_gene_loc.list [[i]]<-which(rowSums(outer(genedata_sum$Gene,lv5_gene.list[[i]],"=="))==1)}

lv1_Folds = list()
lv2_Folds = list()
lv3_Folds = list()
lv4_Folds = list()
lv5_Folds = list()

for (i in 1:lv_nrow[1]){lv1_Folds[[i]]<-sum(genedata_sum$avg_fold_changes[unlist(lv1_gene_loc.list[[i]])])}
for (i in 1:lv_nrow[2]){lv2_Folds[[i]]<-sum(genedata_sum$avg_fold_changes[unlist(lv2_gene_loc.list[[i]])])}
for (i in 1:lv_nrow[3]){lv3_Folds[[i]]<-sum(genedata_sum$avg_fold_changes[unlist(lv3_gene_loc.list[[i]])])}
for (i in 1:lv_nrow[4]){lv4_Folds[[i]]<-sum(genedata_sum$avg_fold_changes[unlist(lv4_gene_loc.list[[i]])])}
for (i in 1:lv_nrow[5]){lv5_Folds[[i]]<-sum(genedata_sum$avg_fold_changes[unlist(lv5_gene_loc.list[[i]])])}

lv1_members$TotalFold = unlist(lv1_Folds)
lv2_members$TotalFold = unlist(lv2_Folds)
lv3_members$TotalFold = unlist(lv3_Folds)
lv4_members$TotalFold = unlist(lv4_Folds)
lv5_members$TotalFold = unlist(lv5_Folds)


lv1_members2 = lv1_members[order(lv1_members$TotalFold,decreasing = T),]
lv2_members2 = lv2_members[order(lv2_members$TotalFold,decreasing = T),]
lv3_members2 = lv3_members[order(lv3_members$TotalFold,decreasing = T),]
lv4_members2 = lv4_members[order(lv4_members$TotalFold,decreasing = T),]
lv5_members2 = lv5_members[order(lv5_members$TotalFold,decreasing = T),]

members =            rbind(lv1_members, lv2_members, lv3_members, lv4_members, lv5_members)
membersbyTotalFold = rbind(lv1_members2,lv2_members2,lv3_members2,lv4_members2,lv5_members2)

### select top 5 & top 10 

levellocation = cumsum(c(0,nrow(lv1_members),
                           nrow(lv2_members),
                           nrow(lv3_members),
                           nrow(lv4_members)))+1

top3list <- list()
top5list <- list()
top10list <- list()

for (i in 1:5){
  top3list[[i]] <- levellocation[i]:(levellocation[i]+3-1)
  top3list[[1]] <- 1
}

for (i in 1:5){
  top5list[[i]] <- levellocation[i]:(levellocation[i]+5-1)
  top5list[[1]] <- 1
}

for (i in 1:5){
  top10list[[i]] <- levellocation[i]:(levellocation[i]+10-1)
  top10list[[1]] <- 1
}

top3list  <- unlist(top3list)
top5list  <- unlist(top5list)
top10list <- unlist(top10list)


### create table

top3P  = members[top3list,-c(4,5)]
top3TF = membersbyTotalFold[top3list,-c(4,5)]
top5P  = members[top5list,-c(4,5)]
top5TF = membersbyTotalFold[top5list,-c(4,5)]

rownames(top5P) <- c()
rownames(top5TF ) <- c()

top5Ptable <- xtable(top5P, type = "latex", file = "top5Ptable.tex")
top5TFtable <- xtable(top5TF, type = "latex", file = "top5TFtable.tex")


print(xtable(top5P, type = "latex"), file = "top5Ptable.tex")
print(xtable(top5TF, type = "latex"), file = "top5TFtable.tex")


top5Ptable 
top5TFtable 

write.table(top5Ptable, file = "top5Ptable.txt")
write.table(top5TFtable, file = "top5TFtable.txt")

#####################

# DAG Plots


pdf("GOGroup1.pdf") 
layout.matrix <- matrix(c(1, 0,
                          1, 3,
                          2, 3,
                          2, 0), ncol = 4)

layout(mat = layout.matrix,
       heights = c(5, 6),    # Heights of the two rows
       widths  = c(1, 2, 2, 1)) # Widths of the two columns


# https://rdrr.io/bioc/GOFunction/src/R/createGODAG.R

plot(createGODAG(c("GO:0044260")),main = "Biological Pathway for GO:0044260 (Level 5)")
plot(createGODAG(c("GO:0050794")),main = "Biological Pathway for GO:0050794 (level 4)")
plot(createGODAG(c("GO:0050794","GO:0044260")),main = "Combined Biological Pathway for GO:0050794 and GO:0044260")
dev.off()


pdf("GOGroup2.pdf") 
par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
##lv2
plot(createGODAG(c("GO:0009987","GO:0065007","GO:0050789")),
     main = "level 2 GO Groups")
##lv3
plot(createGODAG(c("GO:0050789","GO:0071704","GO:0050794")),
     main = "level 3 GO Groups")
##lv4
plot(createGODAG(c("GO:0050794","GO:0043170","GO:0044260")),
     main = "level 4 GO Groups")
##lv5
plot(createGODAG(c("GO:0044260","GO:0031323","GO:0060255")),
     main = "level 5 GO Groups")

# mtext("Combined Biological Pathway for top 3 most significant GO Groups", outer = TRUE)
par(mfrow=c(1,1))
dev.off()

# https://www.biostars.org/p/1272/
# create DAG and do cluster analysis

### cluster analysis example https://guangchuangyu.github.io/
### https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12628










