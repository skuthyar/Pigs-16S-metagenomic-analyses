
# SNP analyses 
setwd("~/Documents/PIG/PIG_Analysis_Chapter_1/SNP")
library(tidyverse) 
library(reshape2)
library(SNPRelate)
library(gdsfmt)
library(SeqArray)
library(readxl)

vcf_file <- snpgdsVCF2GDS(vcf.fn = "plink.vcf", out.fn = "PIG2.gds", method = "copy.num.of.ref")
(genofile <- snpgdsOpen("PIG2.gds"))

meta <- read_excel("SNP_metadata.xlsx")
dom_cat <- meta$dom_cat
genetics <- meta$genetics
location <- meta$location
sampleid <- meta$sampleid

# phylogenetic tree
set.seed(100)
ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram,main="Dendrogram based on SNPs")

# PCA
pca <- snpgdsPCA(genofile, snp.id = NULL, autosome.only = FALSE)

dissMatrix  <-  snpgdsDiss(genofile , sample.id=NULL, snp.id=NULL,
                           autosome.only=TRUE,remove.monosnp=TRUE, maf=NaN, missing.rate=NaN,
                           num.thread=10, verbose=TRUE)
hc <- snpgdsHCluster(dissMatrix)
names(hc)
plot(hc$dendrogram)
SNP_dissMatrix <- as.data.frame(dissMatrix$diss)
write_csv(SNP_dissMatrix, "SNP_dissMatrix.csv")

pc.percent <- pca$varprop*100
pc.percent.one <- round(pc.percent, 2)[1]
pc.percent.two <- round(pc.percent, 2)[2]

ggplot_Data <- pca$eigenvect %>% as.data.frame()

# by domestication context
SNP_domcat_location <- ggplot(ggplot_Data, aes(x=V1, y=V2, color=dom_cat, shape=location)) + geom_point(size=4.5) + scale_color_manual(values = c("industrial"="#963a3a","research"="#ac6d2c","free-ranging domestic"="#52360b","wild"="#36651e"), labels=c("industrial","research","free-ranging domestic","free-ranging wild")) + 
  scale_shape_manual(values=c(0,15,2,16,17,23,12,1,8)) + theme_bw() + xlab("PC1: 50.95%") + ylab("PC2: 6.64%") + theme(legend.title=element_blank()) + theme(legend.text=element_text(size=15))
  
SNP_domcat_location

# by genetics
pdf("PCA.SNPRelate.genetics.pdf")
pca <- snpgdsPCA(genofile, snp.id = NULL, autosome.only = FALSE)
pc.percent <- pca$varprop*100

pc.percent.one <- round(pc.percent, 2)[1]
pc.percent.two <- round(pc.percent, 2)[2]

plot(pca$eigenvect[,1], pca$eigenvect[,2], pch=16, main = "SNP PCA", xlab=paste("PC1", pc.percent.one), ylab=paste("PC2", pc.percent.two), col=as.numeric(as.factor(genetics)))
legend(x="bottomright",legend=unique(as.factor(genetics)), col=unique(factor(genetics)), pch=16)
dev.off()

# by population
pdf("PCA.SNPRelate.location.pdf")
pca <- snpgdsPCA(genofile, snp.id = NULL, autosome.only = FALSE)
pc.percent <- pca$varprop*100

pc.percent.one <- round(pc.percent, 2)[1]
pc.percent.two <- round(pc.percent, 2)[2]

plot(pca$eigenvect[,1], pca$eigenvect[,2], pch=16, main = "SNP PCA", xlab=paste("PC1", pc.percent.one), ylab=paste("PC2", pc.percent.two), col=as.numeric(as.factor(location)))
legend(x="bottomright",legend=unique(as.factor(location)), col=unique(factor(location)), pch=16)
dev.off()

# by individual
pdf("PCA.SNPRelate.sample.pdf")
pca <- snpgdsPCA(genofile, snp.id = NULL, autosome.only = FALSE)
pc.percent <- pca$varprop*100

pc.percent.one <- round(pc.percent, 2)[1]
pc.percent.two <- round(pc.percent, 2)[2]

plot(pca$eigenvect[,1], pca$eigenvect[,2], pch=16, main = "SNP PCA", xlab=paste("PC1", pc.percent.one), ylab=paste("PC2", pc.percent.two), col=as.numeric(as.factor(sampleid)))
legend(x="bottomright",legend=unique(as.factor(sampleid)), col=unique(factor(sampleid)), pch=16)
dev.off()

# structure
structure2 <- read.delim("Pigs.2.structure.2.meanQ", header=TRUE)
structure2_reshape <- melt(structure2)
ggplot(structure2_reshape, aes(fill=variable, y=value, x=genetics)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("grey", "dark blue")) + theme(axis.title.x=element_blank()) + theme(legend.title= element_blank())

structure3 <- read.delim("Pigs.3.structure.3.meanQ", header=TRUE)
structure3_reshape <- melt(structure3)
ggplot(structure3_reshape, aes(fill=variable, y=value, x=genetics)) + 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values = c("grey", "dark blue","purple")) + theme(axis.title.x=element_blank()) + theme(legend.title= element_blank())
