setwd("/proj/b2014034/nobackup/POPULATION_RESEQ/PCA_ADMIXTURE/all/")
library("SNPRelate")
library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)
vcf.fn<-"all_pca.vcf"
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- openfn.gds("ccm.gds")
snpgdsSummary("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- read.table(".pop.group", header = F)


ccm_pca<-snpgdsPCA(genofile)
#names(ccm_pca)
head(cbind(sample.id, pop_code))

snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
names(snpset)
head(snpset$chr1)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=1)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

sample.id = pca$sample.id
pop = factor(pop_code$V1)

tab <- data.frame(sample.id = pca$sample.id, pop = factor(pop_code$V1)[match(pca$sample.id, sample.id)],colour=factor(pop_code$V2)[match(pca$sample.id, sample.id)] , EV1 = pca$eigenvect[,1], EV2 = pca$eigenvect[,2], stringsAsFactors = FALSE)
head(tab)


plot(tab$EV2~tab$EV1, col=as.character(tab$colour), xlab="eigenvector 2", ylab="eigenvector 1",pch=19, cex.lab=1.5)
#with(tab, text(tab$EV2~tab$EV1, labels = as.character(tab$sample.id)), pos = 4)
legend("bottomright", legend=levels(tab$pop), pch=19, col=c("#006400","#C8C800","#0000FF","#FF8C00","#FF0000","#C04000"), cex=1)


lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:5], col=as.character(tab$colour), labels=lbls, pch=19)


