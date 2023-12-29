# load library
library(affy)
library(limma)
library(simpleaffy)
library(scatterplot3d)
library(ggfortify)
library(hgu133a.db)
library(annotate)
library(EnhancedVolcano)
library(clusterProfiler)

# save data
save.image('result.RData')

# set working directory
setwd('~/onine_test_1')

# annotation file
adf <- read.csv('annotation.csv')

# load data
affybatch <- ReadAffy()
colnames(affybatch)
head(exprs(affybatch[,c(1,5)]))
adf$accession
head(exprs(affybatch[,c(5,17)]))

# annotate the sample names
exprs(affybatch) <- exprs(affybatch)[,sapply(adf$accession,grep,colnames(affybatch))]
sampleNames(affybatch) <- adf$status
sampleNames(affybatch)
head(exprs(affybatch[,c(1,5)]))

# quality control of raw data
# QC stats
png('QC stats.png',700,700)
qc_stats <- qc(affybatch)
plot(qc_stats)
dev.off()
# histogram
png('hist_raw_QC.png')
par(mar=c(5,5,2,2))
hist(affybatch, cex.lab = 1.3, cex.axis = 1.3, main='QC: Histogram of log2 Expression (raw data)')
dev.off()
# boxplot
colours <- c(rep('green',9),rep('red',12))
png('box_raw_QC.png',800,800)
par(mar=c(5,5,2,2))
boxplot(affybatch, col=colours, cex.lab = 1.3, cex.axis = 1.3, boxwex=0.5,
        ylab='Log2 Expression', main='QC: Boxplot of Log2 Expression (raw data)')
dev.off()
# PMA flag
calls <- mas5calls(affybatch)
calls <- exprs(calls)
absent <- colSums(calls == 'A')
present <- colSums(calls == 'P')
marginal <- colSums(calls == 'M')
PMA_df <- data.frame(absent,present,marginal)
PMA_df
write.csv(PMA_df,'PMA.csv')
# mva plot
# png('mva_before_norm.png',1000,1000)
# mva.pairs(affy::pm(affybatch), main='Before Normalisation')
# dev.off()
png('mva_NC1_RA1_BeforeNorm.png',700,700)
mva.pairs(affy::pm(affybatch[,c(1,10)]), main='MvA Plot: NC1 vs RA1 (raw data)')
dev.off()

# normalization
eset <- rma(affybatch)  # log2 is performed with RMA
values <- exprs(eset)

# QC after normalization
# boxplot
png('box_RMA_QC.png',800,800)
par(mar=c(5,5,2,2))
boxplot(values, col=colours, cex.lab = 1.3, cex.axis = 1.3, boxwex=0.5, 
        ylab='Log2 Expression',main='QC: Boxplot of Log2 expression (normalised data)')
dev.off()
# mva plot
# png('mva_after_norm.png',1000,1000)
# mva.pairs(values, main='After Normalisation')
# dev.off()
png('mva_NC1_RA1_AfterNorm.png',700,700)
mva.pairs(values[,c(1,10)], main='MvA Plot: NC1 vs RA1 (normalised data)')
dev.off()

# sample similarities
# hierarchical clustering
hc <- hclust(as.dist(1-cor(values, method="pearson")), method="average")
png('hclust.png',700,700)
plot(hc)
dev.off()
# PCA
pca <- prcomp(t(values), scale=T)
pca_sum <- summary(pca)
# 2D plot
cluster <- c(rep('NC',9),rep('RA',12))
autoplot(pca,data=data.frame(t(values),cluster),colour='cluster', 
         main='2D PCA') + 
  geom_text(label=colnames(values),nudge_x=0.015, nudge_y=0.015, size=2.5)
ggsave('PCA_2D.jpg',device='jpg', width = 6, height = 5)
# 3D plot
# calculate variation proportion of each PC
pca_1 <- round(pca_sum$importance[2,'PC1'],4)*100
pca_1 <- paste('PC1 (',pca_1,'%)',sep='')
pca_2 <- round(pca_sum$importance[2,'PC2'],4)*100
pca_2 <- paste('PC2 (',pca_2,'%)',sep='')
pca_3 <- round(pca_sum$importance[2,'PC3'],4)*100
pca_3 <- paste('PC3 (',pca_3,'%)',sep='')
png('PCA_3D.png',700,700)
s3d <- scatterplot3d(pca$x[,1:3], pch=19, color=colours,
                     xlab = pca_1, ylab = pca_2, zlab = pca_3, main='3D PCA')
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(values),
     pos = 3,offset = 0.5)
dev.off()

# gene names annotation
tmp <- AnnotationDbi::select(hgu133a.db,keys=featureNames(eset),
              columns=c('SYMBOL','GENENAME'),keytypes='PROBEID')
#sum(duplicated(tmp$PROBEID))
tmp <- tmp[match(featureNames(eset), tmp$PROBEID),]  # only use the first matched ID and drop duplicated ID
nrow(tmp[rowSums(is.na(tmp))>0,])
tmp[tmp=="NA"] <- NA
tmp[is.na(tmp)] <- 'FailAnno' 
fData(eset) <- tmp

# top 10 fold change genes
vals10 <- 2^values
NC_mean <- apply(vals10[,paste0('NC',1:9)],1,mean)
NC_sd <- apply(vals10[,paste0('NC',1:9)],1,sd)
RA_mean <- apply(vals10[,paste0('RA',1:12)],1,mean)
RA_sd <- apply(vals10[,paste0('RA',1:12)],1,sd)
cal_p <- function(row){
  return(t.test(row[paste0('NC',1:9)], row[paste0('RA',1:12)],
                "two.sided",paired = FALSE)$p.value)
}
p_val <- apply(vals10,1,cal_p)
fdrs <- p.adjust(p_val, method="BH")
RA_NC_ratio <- log2(RA_mean/NC_mean)
FC_df <- as.data.frame(cbind(RA_mean,RA_sd,NC_mean,NC_sd,RA_NC_ratio,p_val,fdrs))
FC_df <- FC_df[order(abs(FC_df$RA_NC_ratio),decreasing=TRUE),]
FC_df <- cbind(tmp[match(rownames(FC_df), tmp$PROBEID),],FC_df)
head(FC_df,10)
write.csv(head(FC_df,10),'FC_ratio.csv')
FC_df_fdr <- FC_df[order(FC_df$fdrs,decreasing=FALSE),]
head(FC_df_fdr,10)
write.csv(head(FC_df_fdr,10),'FC_fdr.csv')


# Differentially Expressed Genes
# design matrix
design <- model.matrix(~0+factor(c(rep(1,9),rep(2,12))))
colnames(design) <- c("NC","RA")
design
# contrast matrix
contrastmatrix <- makeContrasts(RA-NC,levels=design)
contrastmatrix
# DE genes by linear model
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrastmatrix)
# Find DE by function 'treat' (more strict compared to function 'ebayes')
fit_treat <- treat(fit2,lfc = log2(1.1), trend = TRUE, robust = TRUE)  # fold change threshold: log2(1.1)
DE_treat <- topTreat(fit_treat, coef = 1, number = Inf, p.value = 0.05, adjust.method = "BH")
head(DE_treat)
nrow(DE_treat)
all_treat <- topTreat(fit_treat, coef = 1, number = Inf, adjust.method = "BH")
nrow(all_treat)
write.csv(DE_treat,"DE_treat.csv")
# volcano plot
png('volcano_treat.png',1000,1000)
EnhancedVolcano(all_treat, lab = all_treat$SYMBOL, x = 'logFC', y = 'adj.P.Val',
                pCutoff = 0.05, FCcutoff = log2(1.1),drawConnectors = TRUE,
                widthConnectors = 0.75,title = "Volcano plot (thres: log2(1.1), p.adjust=0.05)",
                subtitle = "by function 'treat'",
                max.overlaps = 10)
dev.off()

# enrichment analysis
# add 'EntrezID' to feature annotation map
tmp <- select(hgu133a.db,keys=featureNames(eset),
              columns=c('ENTREZID','SYMBOL','GENENAME'),keytypes='PROBEID')
tmp <- tmp[match(featureNames(eset), tmp$PROBEID),]
tmp[tmp=="NA"] <- NA
fData(eset) <- tmp
# remove probes that cannot find EntrezID
sum(is.na(fData(eset)$ENTREZID))
eset_t <- eset[!is.na(fData(eset)$ENTREZID),]
head(fData(eset_t))
length(fData(eset_t)$ENTREZID)
DEG_id <- fData(eset_t)$ENTREZID[fData(eset_t)$PROBEID %in% DE_treat$PROBEID]
length(DEG_id)
head(DEG_id)
DE_treat$PROBEID[!DE_treat$PROBEID %in% fData(eset_t)$PROBEID]  # Probes without EntrezID
# kegg pathway
enrich_kegg <- enrichKEGG(gene = DEG_id,
                         keyType = "kegg",
                         organism = 'hsa',  # homo sapiens
                         pAdjustMethod = "fdr",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.1  # default is 0.2
)
enrich_kegg <- data.frame(enrich_kegg)
enrich_kegg$ratio <- sapply(enrich_kegg$GeneRatio, FUN = function(x){eval(parse(text=x))})
write.csv(enrich_kegg,'enrich_kegg.csv',row.names = TRUE)
# plot KEGG enrichment
ggplot(enrich_kegg %>%
         dplyr::slice_min(n = 30, order_by = p.adjust),
       aes(y=reorder(Description,-p.adjust),x=ratio))+
  geom_point(aes(size=Count,color=p.adjust))+
  scale_color_gradient(low = "red",high ="blue")+
  labs(color=expression(p.adjust,size="Count"),
       x="Gene Ratio",y="KEGG term",title="Top 30 KEGG Enrichment")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12))
ggsave('enrich_kegg.jpg',device='jpg', width = 10, height = 8)


# get Human gene list signatures from Molecular Signatures Database (MSigDB)
system("wget http://bioinf.wehi.edu.au/software/MSigDB/human_H_v5p2.rdata") #hallmark
load("human_H_v5p2.rdata")
#convert to indexes
H.indices <- ids2indices(Hs.H,fData(eset_t)$ENTREZID)

# self-contained enrichment
result_mroast <- mroast(eset_t, index=H.indices, design=design,
                        contrast=contrastmatrix[,1], adjust.method = "BH")
head(result_mroast,20)
write.csv(head(result_mroast,20),'mroast.csv')

# competitive enrichment
result_camera <- camera(eset_t, index=H.indices, design=design, 
                        contrast=contrastmatrix[,1])
head(result_camera,20)
write.csv(head(result_camera,20),'camera.csv')

# GSEA like competitive enrichment
result_romer <- romer(eset_t, index=H.indices, design=design, 
                      contrast=contrastmatrix[,1])
romer_df <- as.data.frame(result_romer)
head(romer_df[order(romer_df$Mixed,decreasing=FALSE),],20)
write.csv(head(romer_df[order(romer_df$Mixed,decreasing=FALSE),],20),'romer.csv')
