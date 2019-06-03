#################################
# Drafts
#########################################

### add numeric ratio
top = as.numeric(unlist(strsplit(as.character(members$GeneRatio),"/"))[c(T,F)])
base = as.numeric(unlist(strsplit(as.character(members$GeneRatio),"/"))[c(F,T)])

members$GeneRatioN  = top/base
head(members)[-5]

### add avg_folds and per group





#END HERE



###SELECT 10% MOST SIG GROUPS AND COUNT NUMBER OF REPEATED GENE?  MAYBE WE WILL GET SIMILAR RESULT WITH GENE ENRICHMENT ANALYSIS








#########################################################################

gene.df <- bitr(genedata_sum$Gene, fromType = "SYMBOL", toType = c("GOALL"),OrgDb = org.Hs.eg.db)
gene.df <- subset(gene.df, ONTOLOGYALL == 'BP')
head(gene.df)
nrow(gene.df)




plotGOgraph #need enrichment output

data(geneList, package = "DOSE")
de <- names(geneList)[1:100]
yy <- enrichGO(de, 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
head(yy)







####################################################################

# https://www.bioconductor.org/packages/devel/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
# http://jdblischak.github.io/nw/analysis/mouse/go.html
# https://support.bioconductor.org/p/96577/
# mart = useMart('ensembl')
# listDatasets(mart)
listAttributes
# Get ensembl gene ids and GO terms

human = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
test <- getBM(attributes = c("description","ensembl_gene_id","go_id"), values = genedata$Gene, mart = human)
go_id <- getBM(attributes = c("go_id"), values = genedata$Gene, mart = human)

genedata$ensembl_gene_id  = ensembl_gene_id 
genedata$go_id = go_id
genedata [1:50,]


############################################################################


# classify into go term

# draw a go tree

# gene set enrichment???
#https://www.bioconductor.org/help/course-materials/2010/SeattleIntro/Bioconductor-Introduction.pdf

# 2 sample t test, using hierarchical bayes

# get posterior

# comapre theta =0 and theta > 0 
