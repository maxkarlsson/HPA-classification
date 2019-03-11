library('pathview')
library('DESeq2')
library('RUVSeq')
library('gage')
library('gageData')
library('VennDiagram')
library("AnnotationDbi")
library("org.Hs.eg.db")
#library('tidyverse')
library('EnhancedVolcano')
library('ggthemr')
library('ggbeeswarm')
library('gridExtra')
library('vsn')
library('RColorBrewer')
library('pcaExplorer')
library('ggrepel')
library('stringr')
library('clusterProfiler')


setwd('/Users/wen.zhong/Work/BloodAtlas/')
result_folder <- paste0("./result/", Sys.Date())
dir.create(result_folder, showWarnings = FALSE)

gene.info <- read.delim("../database/ensembl92/geneinfo_92.txt",header = TRUE, sep='\t',check.names = T)
#protein.class <- read.table(paste(workpath,'RNAseq','annotation','gene.classes.txt',sep='/'),header = TRUE, sep='\t')
#rownames(protein.class) <- protein.class$rna.genes

## kegg
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]

## GO
data(go.sets.hs)
data(go.subs.hs)

## blood counts
counts <- read.table('./data/summarizedTPM_count_gene_level_rename.txt',sep='\t',header=TRUE,row.names=1)
pdata <- read.table('./data/sample_info.txt',sep='\t',header=TRUE,row.names = 1)
counts.pseudo <- round(counts)
blood.Set <- newSeqExpressionSet(as.matrix(counts.pseudo),
                                 phenoData = data.frame(pdata))

blood.atlas.tpm <- read.table('./data/summarizedTPM_proteinCoding_cellTypes_ensmbl92.txt',sep='\t',header=TRUE,row.names=1)
blood.atlas.tpm <- blood.atlas.tpm[,20:132]
pdata <- pdata[match(colnames(blood.atlas.tpm),rownames(pdata)),]
blood.tpm.Set <- newSeqExpressionSet(as.matrix(blood.atlas.tpm),phenoData = data.frame(pdata))

## sample remove
qcfailed <- c('bd6_myeloid_dc','bd2_non_classical_monocytes','bd5_memory_cd8','bd1_basophils')
blood.Set <- blood.Set[,!(sampleNames(blood.Set) %in% qcfailed)]
blood.tpm.Set <- blood.tpm.Set[,!(sampleNames(blood.tpm.Set) %in% qcfailed)]

## sample select
# combn(levels(pdata$Cell), 2, function(x){
#   sample1 <- x[1]
#   sample2 <- x[2]
# })
sample1 <- 'tregs'
sample2 <- 'memory_cd8'
#sample1 <- 'B_cells'
#sample2 <- 'T_cells'

sample1_sample2_folder <- paste0("./result/", Sys.Date(), "/", sample1,'_',sample2)
dir.create(sample1_sample2_folder, showWarnings = FALSE)

blood.Set.select <- blood.Set[,blood.Set$Cell == sample1|blood.Set$Cell == sample2]
#blood.Set.select <- blood.Set[,blood.Set$Type == sample1|blood.Set$Type == sample2]
filter <- apply(counts(blood.Set.select),1,function(x) length(x[x>5])>=2)
genes.filter <- counts(blood.Set.select)[filter,]
pheno.select <- pData(blood.Set.select)
pheno.select$Cell <- factor(pheno.select$Cell)
#pheno.select$Type <- factor(pheno.select$Type)
blood.Set.select <- newSeqExpressionSet(genes.filter,
                                        phenoData = pheno.select)
idxbycell <- order(blood.Set.select$Cell)
#idxbytype<- order(blood.Set.select$Type)
blood.Set.select <- blood.Set.select[,idxbycell]
#blood.Set.select <- blood.Set.select[,idxbytype]

blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Cell == sample1|blood.tpm.Set$Cell == sample2]
#blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Type == sample1|blood.tpm.Set$Type == sample2]
pheno.tpm.select <- pData(blood.tpm.Set.select)
pheno.tpm.select$Cell <- factor(pheno.tpm.select$Cell)
#pheno.tpm.select$Type <- factor(pheno.tpm.select$Type)
blood.tpm.Set.select <- newSeqExpressionSet(counts(blood.tpm.Set.select),
                                            phenoData = pheno.tpm.select)

## PCA
colors <- c("#f1a340","#998ec3")
x <- blood.Set.select$Cell
#x <- blood.Set.select$Type
#plotRLE(blood.Set.select,outline=FALSE, ylim=c(-4, 4),col=colors[x])
#plotPCA(blood.Set.select, col=colors[x], cex=1.2)

## RUV
design <- model.matrix(~Cell, data=pData(blood.Set.select))
#design <- model.matrix(~Type, data=pData(blood.Set.select))
y <- DGEList(counts=counts(blood.Set.select), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
top <- topTags(lrt, n=nrow(blood.Set.select))$table
#empirical <- rownames(top[top$FDR>0.5,])
empirical <- rownames(blood.Set.select)[which(!(rownames(blood.Set.select) %in% rownames(top)[1:5000]))]

blood.Set.select.RUV <- RUVg(blood.Set.select, empirical, k=1)
#plotRLE(blood.Set.select.RUV, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#plotPCA(blood.Set.select.RUV, col=colors[x], cex=1.2)

## set cutoff
cutoff_fc <- 1
cutoff_adjP <- 0.01

## DESeq2 original
dds <- DESeqDataSetFromMatrix(countData = counts(blood.Set.select),
                              colData = pData(blood.Set.select),
#                              design = ~ Type)
                              design = ~ Cell)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
colnames(vsd) <- factor(colnames(vsd), levels = colnames(vsd))

#meanSdPlot(assay(vsd))
res <- results(dds, contrast=c("Cell",sample2,sample1))
#res <- results(dds, contrast=c("Type",sample2,sample1))
res <- na.omit(res)
res.symbol <- merge(data.frame(res), gene.info, by.x="row.names", by.y="ensg_id",all.x=TRUE)
res.Ordered <- res.symbol[order(res.symbol$pvalue),]
res.Ordered$entrez <- mapIds(org.Hs.eg.db,
                             keys=res.Ordered$Row.names, 
                             column="ENTREZID",
                             keytype="ENSEMBL",
                             multiVals="first")
res.sig <- res.Ordered[res.Ordered$padj<cutoff_adjP & abs(res.Ordered$log2FoldChange)>cutoff_fc,]
sig_DEGs <- res.sig[order(abs(res.sig$log2FoldChange),decreasing = TRUE),]$Row.names
sig.up <- res.sig[res.sig$log2FoldChange>0,]
sig.down <- res.sig[res.sig$log2FoldChange<0,]

DESeq_result_folder <- paste0(sample1_sample2_folder,"/DESeq")
dir.create(DESeq_result_folder, showWarnings = FALSE)
write.table(data.frame("ensgid"=rownames(res.Ordered),res.Ordered),
            file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_pvalue.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

write.table(data.frame("ensgid"=rownames(res.sig),res.sig),
            file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_sig_pvalue.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_basemean <- sapply(levels(dds$Cell), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$Cell == lvl]))
#DESeq_basemean <- sapply(levels(dds$Type), function(lvl) rowMeans(counts(dds,normalized=TRUE)[,dds$Type == lvl]))
colnames(DESeq_basemean) <- c('sample2','sample1') 
DESeq_basemean <- merge(DESeq_basemean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
DESeq_basemean$sig <- ifelse(DESeq_basemean$Row.names %in% res.sig$Row.names,'DEG','no')
pdf(file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_basemean.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_basemean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
        geom_text(data=DESeq_basemean[DESeq_basemean$Row.names %in% sig_DEGs,],
                        aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                        size=2.8,check_overlap = TRUE)+
        labs(title = 'BaseMean: DESeq normalized counts',
             subtitle = paste0('up_DEGs:',nrow(sig.up),'  ','down_DEGs:',nrow(sig.down)))+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        scale_color_manual(values=c("#e41a1c","#999999"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

DESeq_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cell), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cell == lvl]))
#DESeq_tpm_mean <- sapply(levels(blood.tpm.Set.select$Type), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Type == lvl]))
colnames(DESeq_tpm_mean) <- c('sample2','sample1') 
DESeq_tpm_mean <- merge(DESeq_tpm_mean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
DESeq_tpm_mean$sig <- ifelse(DESeq_tpm_mean$Row.names %in% res.sig$Row.names,'DEG','no')
pdf(file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_tpm_mean.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_tpm_mean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
        geom_text(data=DESeq_tpm_mean[DESeq_tpm_mean$Row.names %in% sig_DEGs,],
                        aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                        size=2.8,segment.size=0.2,check_overlap = TRUE)+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        labs(title = 'TPM: DESeq',
             subtitle = paste0('up_DEGs:',nrow(sig.up),'  ','down_DEGs:',nrow(sig.down)))+
        scale_color_manual(values=c("#e41a1c","#999999"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()


pdf(file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_vsd_pca.pdf'),width=7,height=7, useDingbats = F)
print(pcaplot(vsd,intgroup=c('Cell'),ntop=1000,
#print(pcaplot(vsd,intgroup=c('Type'),ntop=1000,
              pcX=1, pcY=2,title='PCA: DESeq',text_labels = FALSE, ellipse = TRUE)+
        geom_text_repel(aes(label=paste(vsd$Cell,vsd$Subject,sep='.')))+
        coord_fixed(ratio=0.8)+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

vsd.counts <- tidyr::gather(as.data.frame(assay(vsd)))
names(vsd.counts) <- c("Sample", "vstExpression")
vsd.counts$Sample <- factor(vsd.counts$Sample, levels = colnames(vsd))
des <- colData(vsd)
vsd.counts <- as.data.frame(merge(vsd.counts, des, by.x="Sample", by.y="row.names"))

pdf(file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_vsd_boxplot.pdf'),width=5,height=5, useDingbats = F)
print(ggplot(vsd.counts, aes_string(x = "Sample", y = "vstExpression")) + 
        geom_boxplot(aes_string(color= "Cell", fill = "Cell"), 
#        geom_boxplot(aes_string(color= "Type", fill = "Type"), 
                     alpha = 0.5)+
        ggtitle('Boxplot: DESeq')+
        scale_fill_manual(values=c("#fc8d59", "#99d594"))+
        scale_color_manual(values=c("#fc8d59", "#99d594"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text.x = element_text(angle=45,hjust=1),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.7, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

#ggthemr_reset()
p.volcano <- EnhancedVolcano(res.symbol,
                             lab = res.symbol$gene_name,
                             x = "log2FoldChange",
                             y = "padj",
                             xlab = bquote(~Log[2]~ "fold change"),
                             ylab = bquote(~-Log[10]~adjusted~italic(P)),
                             pCutoff = cutoff_adjP,
                             FCcutoff = cutoff_fc,
                             #xlim=c(-6,6),
                             transcriptLabSize = 3.0,
                             colAlpha = 1,
                             #title = paste0("Method:DESeq2",'\n','cutoff: ','log2FC, ',cutoff_fc,' adjP, ',cutoff_adjP),
                             legend=c("NS","Log2 FC","Adjusted p-value",
                                      "Adjusted p-value & Log2 FC"),
                             legendPosition = "bottom",
                             legendLabSize = 12,
                             legendIconSize = 3.0)+
  geom_text(aes(label=sample2, x= max(res.symbol$log2FoldChange)/2,y=2*max(-log10(res.symbol$padj))/3),size=5,color='#d95f02',alpha=0.5)+
  geom_text(aes(label=sample1, x= min(res.symbol$log2FoldChange)/2,y=2*max(-log10(res.symbol$padj))/3),size=5,color='#d95f02',alpha=0.5)
pdf(file=paste0(DESeq_result_folder,"/",sample1,'_',sample2,'_DESeq_volcano.pdf'),width=8,height=8, useDingbats = F)
print(grid.arrange(p.volcano, top=textGrob(paste0("Method: DESeq2, up_DEGs: ",nrow(sig.up),", down_DEGs:",nrow(sig.down),'\n','cutoff:  ','|log2FC|>',cutoff_fc,', adjP<',cutoff_adjP), gp=gpar(fontsize=15,fontface="bold"))))
dev.off()
#2*max(-log10(res.symbol$padj))/3
#2*max(-log10(res.symbol$padj))/3
# ## enrichment 
DESeq_result_enrich_folder <- paste0(sample1_sample2_folder,"/DESeq/enrichment/")
dir.create(DESeq_result_enrich_folder, showWarnings = FALSE)
DEGs_sig_all <- res.sig$entrez
DEGs_sig_up <- res.sig[res.sig$log2FoldChange>0,]$entrez
DEGs_sig_down <- res.sig[res.sig$log2FoldChange<0,]$entrez
DEGs_df <- data.frame(Entrez=res.sig$entrez,FC=res.sig$log2FoldChange)
DEGs_df$group <- ifelse(DEGs_df$FC>0,"upregulated","downregulated")

# enrich KEGG
DESeq_result_enrich_KEGG_folder <- paste0(sample1_sample2_folder,"/DESeq/enrichment/KEGG")
dir.create(DESeq_result_enrich_KEGG_folder , showWarnings = FALSE)
DESeq_enrich_KEGG <- compareCluster(Entrez~group, data=DEGs_df, fun="enrichKEGG",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_enrich_KEGG)),as.data.frame(DESeq_enrich_KEGG)),
            file=paste0(DESeq_result_enrich_KEGG_folder,"/",sample1,'_',sample2,'_enrich_KEGG.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

pdf(file=paste0(DESeq_result_enrich_KEGG_folder,"/",sample1,'_',sample2,'_enrich_KEGG_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_enrich_KEGG,showCategory = 30)+
  labs(title='KEGG enrichment anlaysis',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_enrich_KEGG))))
dev.off()

# enrich GO
DESeq_result_enrich_GO_folder <- paste0(sample1_sample2_folder,"/DESeq/enrichment/GO")
dir.create(DESeq_result_enrich_GO_folder , showWarnings = FALSE)

DESeq_enrich_GO_BP <- compareCluster(Entrez~group, data=DEGs_df, 'org.Hs.eg.db',ont = "BP",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_enrich_GO_BP)),as.data.frame(DESeq_enrich_GO_BP)),
            file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_BP.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_BP_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_enrich_GO_BP,showCategory = 30)+
  labs(title='GO enrichment anlaysis, BP',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_enrich_GO_BP))))
dev.off()

DESeq_enrich_GO_MF <- compareCluster(Entrez~group, data=DEGs_df, 'org.Hs.eg.db',ont = "MF",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_enrich_GO_MF)),as.data.frame(DESeq_enrich_GO_MF)),
            file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_MF.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_MF_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_enrich_GO_MF,showCategory = 30)+
  labs(title='GO enrichment anlaysis, MF',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_enrich_GO_MF))))
dev.off()

DESeq_enrich_GO_CC <- compareCluster(Entrez~group, data=DEGs_df, 'org.Hs.eg.db',ont = "CC",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_enrich_GO_CC)),as.data.frame(DESeq_enrich_GO_CC)),
            file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_CC.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_CC_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_enrich_GO_CC,showCategory = 30)+
  labs(title='GO enrichment anlaysis, CC',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_enrich_GO_CC))))
dev.off()

## GAGE
DESeq_result_gage_folder <- paste0(sample1_sample2_folder,"/DESeq/gage/")
dir.create(DESeq_result_gage_folder, showWarnings = FALSE)
DESeq_foldchanges = res.Ordered$log2FoldChange
DESeq_foldchanges_sorted <- sort(DESeq_foldchanges,decreasing = TRUE)
names(DESeq_foldchanges) = res.Ordered$entrez

# GAGE KEGG
DESeq_result_gage_KEGG_folder <- paste0(sample1_sample2_folder,"/DESeq/gage/KEGG/")
dir.create(DESeq_result_gage_KEGG_folder, showWarnings = FALSE)

DESeq_keggres = gage(DESeq_foldchanges, gsets=kegg.sets.hs, same.dir=FALSE)
write.table(data.frame("keggid"=rownames(DESeq_keggres$greater),DESeq_keggres$greater),
            file=paste0(DESeq_result_gage_KEGG_folder,"/",sample1,'_',sample2,'_gage_KEGG.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_keggres_select <- DESeq_keggres$greater[, "q.val"] < 0.05 & !is.na(DESeq_keggres$greater[, "q.val"])
DESeq_keggres_sigids <- rownames(DESeq_keggres$greater)[DESeq_keggres_select]
DESeq_keggres_sigids <- substr(DESeq_keggres_sigids, 1, 8)

setwd(DESeq_result_gage_KEGG_folder)
dir.create(paste('/Users/wen.zhong/Work/BloodAtlas/result',Sys.Date(),paste0(sample1,'_',sample2),'DESeq/gage/KEGG/','tmp',sep='/'), showWarnings = FALSE)
for(pwid in DESeq_keggres_sigids){
  pathview(gene.data=DESeq_foldchanges, 
           pathway.id=pwid, 
           kegg.native = TRUE,
           out.suffix=paste0('DESeq_',sample1,'_',sample2),
           kegg.dir = paste('/Users/wen.zhong/Work/BloodAtlas/result',Sys.Date(),paste0(sample1,'_',sample2),'DESeq/gage/KEGG/','tmp',sep='/'),
           species="hsa")}
setwd('/Users/wen.zhong/Work/BloodAtlas/')

# GAGE GO
DESeq_result_gage_GO_folder <- paste0(sample1_sample2_folder,"/DESeq/gage/GO/")
dir.create(DESeq_result_gage_GO_folder, showWarnings = FALSE)

DESeq_GO_MF <- gage(DESeq_foldchanges, gsets = go.sets.hs[go.subs.hs$MF],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_GO_MF$greater),DESeq_GO_MF$greater),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_MF_greater.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_GO_MF$less),DESeq_GO_MF$less),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_MF_less.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
 
DESeq_GO_BP <- gage(DESeq_foldchanges, gsets = go.sets.hs[go.subs.hs$BP],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_GO_BP$greater),DESeq_GO_BP$greater),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_BP_greater.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_GO_BP$less),DESeq_GO_BP$less),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_BP_less.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
 
DESeq_GO_CC <- gage(DESeq_foldchanges, gsets = go.sets.hs[go.subs.hs$CC],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_GO_CC$greater),DESeq_GO_CC$greater),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_CC_greater.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_GO_CC$less),DESeq_GO_CC$less),
             file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_CC_less.txt'),
             sep='\t',
             row.names=FALSE,
             quote=FALSE)
 
cutoff_GO_p <- 0.01
# 
DESeq_GO_select_BP_great <- as.data.frame(subset(DESeq_GO_BP$greater,DESeq_GO_BP$greater[,'p.val'] < cutoff_GO_p))
DESeq_GO_select_BP_great$type <- rep('up',nrow(DESeq_GO_select_BP_great))
DESeq_GO_select_BP_great$category <- rep('BP',nrow(DESeq_GO_select_BP_great))
DESeq_GO_select_BP_less <- as.data.frame(subset(DESeq_GO_BP$less,DESeq_GO_BP$less[,'p.val'] < cutoff_GO_p))
DESeq_GO_select_BP_less$type <- rep('down',nrow(DESeq_GO_select_BP_less))
DESeq_GO_select_BP_less$category <- rep('BP',nrow(DESeq_GO_select_BP_less))
 
DESeq_GO_select_MF_great <- as.data.frame(subset(DESeq_GO_MF$greater,DESeq_GO_MF$greater[,'p.val'] < cutoff_GO_p))
DESeq_GO_select_MF_great$type <- rep('up',nrow(DESeq_GO_select_MF_great))
DESeq_GO_select_MF_great$category <- rep('MF',nrow(DESeq_GO_select_MF_great))
DESeq_GO_select_MF_less <- as.data.frame(subset(DESeq_GO_MF$less,DESeq_GO_MF$less[,'p.val'] < cutoff_GO_p))
DESeq_GO_select_MF_less$type <- rep('down',nrow(DESeq_GO_select_MF_less))
DESeq_GO_select_MF_less$category <- rep('MF',nrow(DESeq_GO_select_MF_less))
# 
DESeq_GO_select_CC_great <- as.data.frame(subset(DESeq_GO_CC$greater,DESeq_GO_CC$greater[,'p.val'] < cutoff_GO_p))
DESeq_GO_select_CC_great$type <- rep('up',nrow(DESeq_GO_select_CC_great))
DESeq_GO_select_CC_great$category <- rep('CC',nrow(DESeq_GO_select_CC_great))
DESeq_GO_select_CC_less <- as.data.frame(subset(DESeq_GO_CC$less,DESeq_GO_CC$less[,'q.val'] < cutoff_GO_p))
DESeq_GO_select_CC_less$type <- rep('down',nrow(DESeq_GO_select_CC_less))
DESeq_GO_select_CC_less$category <- rep('CC',nrow(DESeq_GO_select_CC_less))
DESeq_GO_select <- rbind(DESeq_GO_select_BP_great,
                         DESeq_GO_select_BP_less,
                         DESeq_GO_select_MF_great,
                         DESeq_GO_select_MF_less,
                         DESeq_GO_select_CC_great,
                         DESeq_GO_select_CC_less)
DESeq_GO_select$GO <- str_split_fixed(rownames(DESeq_GO_select),' ',2)[,2]
 
pdf(file=paste0(DESeq_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO.pdf'),width=10,height=10, useDingbats = F)
# #  ggthemr('light')
print(ggplot(DESeq_GO_select)+
        geom_point(aes(x=stat.mean, y=-log10(p.val),size=-log10(p.val), fill=category, color=category))+
        geom_text_repel(mapping = aes(x=stat.mean, y=-log10(p.val),label=ifelse(q.val<0.05,GO,''),color=category),size=3)+
        ggtitle(paste0('GO: DESeq','\n','cutoff: p<',cutoff_adjP,'(dots), adjp<0.05(text)'))+
        geom_text(aes(label=sample2, x= max(stat.mean)/2,y=5*max(-log10(p.val))/6),size=5,color='#984ea3')+
        geom_text(aes(label=sample1, x= min(stat.mean)/2,y=5*max(-log10(p.val))/6),size=5,color='#984ea3')+
        geom_vline(xintercept = 0, color='black',linetype="dotted",size=1)+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text.x = element_text(angle=45,hjust=1),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
               legend.key.size= unit(0.7, "cm"),
               legend.title = element_text(face="italic",size=12),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")))
 # ggthemr_reset()
 dev.off()



## DESeq2 + RUV
dds.RUV <- DESeqDataSetFromMatrix(countData = counts(blood.Set.select.RUV),
                                  colData = pData(blood.Set.select.RUV),
#                                  design = ~ W_1+Type)
                                  design = ~ W_1+Cell)
dds.RUV <- DESeq(dds.RUV)
vsd.RUV <- vst(dds.RUV, blind=FALSE)
res.RUV <- results(dds.RUV, contrast=c("Cell",sample2,sample1))
#res.RUV <- results(dds.RUV, contrast=c("Type",sample2,sample1))
res.RUV <- na.omit(res.RUV)
res.RUV.symbol <- merge(data.frame(res.RUV), gene.info, by.x="row.names", by.y="ensg_id",all.x=TRUE)
res.RUV.Ordered <- res.RUV.symbol[order(res.RUV.symbol$pvalue),]
res.RUV.Ordered$entrez <- mapIds(org.Hs.eg.db,
                                 keys=res.RUV.Ordered$Row.names, 
                                 column="ENTREZID",
                                 keytype="ENSEMBL",
                                 multiVals="first")
res.RUV.sig <- res.RUV.Ordered[res.RUV.Ordered$padj<cutoff_adjP & abs(res.RUV.Ordered$log2FoldChange)>cutoff_fc,]
sig_DEGs_RUV <- res.RUV.sig[order(abs(res.RUV.sig$log2FoldChange),decreasing = TRUE),]$Row.names
sig_RUV.up <- res.sig[res.RUV.sig$log2FoldChange>0,]
sig_RUV.down <- res.sig[res.RUV.sig$log2FoldChange<0,]


DESeq_RUV_result_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/")
dir.create(DESeq_RUV_result_folder, showWarnings = FALSE)
write.table(data.frame("ensgid"=rownames(res.RUV.Ordered),res.RUV.Ordered),
            file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_pvalue.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
write.table(data.frame("ensgid"=rownames(res.RUV.sig),res.RUV.sig),
            file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_sig_pvalue.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_RUV_basemean <- sapply(levels(dds.RUV$Cell), function(lvl) rowMeans(counts(dds.RUV,normalized=TRUE)[,dds.RUV$Cell == lvl]))
#DESeq_RUV_basemean <- sapply(levels(dds.RUV$Type), function(lvl) rowMeans(counts(dds.RUV,normalized=TRUE)[,dds.RUV$Type == lvl]))
colnames(DESeq_RUV_basemean) <- c('sample2','sample1') 
DESeq_RUV_basemean <- merge(DESeq_RUV_basemean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
DESeq_RUV_basemean$sig <- ifelse(DESeq_RUV_basemean$Row.names %in% res.RUV.sig$Row.names,'DEG','no')
pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_basemean.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_RUV_basemean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
        geom_text(data=DESeq_RUV_basemean[DESeq_RUV_basemean$Row.names %in% sig_DEGs_RUV,],
                        aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                        size=2.8, check_overlap = TRUE)+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        labs(title = 'BaseMean: DESeq+RUV normalized counts',
             subtitle = paste0('up_DEGs:',nrow(sig_RUV.up),'  ','down_DEGs:',nrow(sig_RUV.down)))+
        scale_color_manual(values=c("#e41a1c","#999999"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

DESeq_DEGs_nonRUV <- sig_DEGs[!(sig_DEGs %in% sig_DEGs_RUV)]
DESeq_RUV_DEGs_nonDEseq <- sig_DEGs_RUV[!(sig_DEGs_RUV %in% sig_DEGs)]
DEGs_nonoverlap <- append(DESeq_DEGs_nonRUV,DESeq_RUV_DEGs_nonDEseq)
DEGs_overlap <- intersect(sig_DEGs, sig_DEGs_RUV)
DESeq_RUV_basemean$compare <- ifelse(DESeq_RUV_basemean$Row.names %in% DESeq_DEGs_nonRUV, "DESeq",
                                     ifelse(DESeq_RUV_basemean$Row.names %in% DESeq_RUV_DEGs_nonDEseq, "DESeq_RUV",
                                            ifelse(DESeq_RUV_basemean$Row.names %in% DEGs_overlap,"overlap_DEGs","none")))
DESeq_RUV_basemean$compare <- factor(DESeq_RUV_basemean$compare,levels=c('DESeq','DESeq_RUV','overlap_DEGs','none'))


pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_non_overlapDEG_basemean_adjp0.01.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_RUV_basemean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=compare),size=1.5,alpha=0.8)+
        geom_text(data=DESeq_RUV_basemean[DESeq_RUV_basemean$Row.names %in% DEGs_nonoverlap,],
                  aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                  size=3, check_overlap = TRUE)+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        labs(title = 'BaseMean: DESeq vs DESeq+RUV',
             subtitle = paste0('DEG_specific:',length(DESeq_DEGs_nonRUV),'  ','DEG_RUV_specific:',length(DESeq_RUV_DEGs_nonDEseq)))+
        scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","lightgrey"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                           size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cell), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cell == lvl]))
#DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Type), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Type == lvl]))
colnames(DESeq_RUV_tpm_mean) <- c('sample2','sample1') 
DESeq_RUV_tpm_mean <- merge(DESeq_RUV_tpm_mean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
DESeq_RUV_tpm_mean$sig <- ifelse(DESeq_RUV_tpm_mean$Row.names %in% res.RUV.sig$Row.names,'DEG','no')
pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_tpm_mean.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_RUV_tpm_mean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
        geom_text(data=DESeq_RUV_tpm_mean[DESeq_RUV_tpm_mean$Row.names %in% sig_DEGs_RUV,],
                        aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                        size=2.8, check_overlap = TRUE)+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        labs(title = 'TPM: DESeq+RUV',
             subtitle = paste0('up_DEGs:',nrow(sig_RUV.up),'  ','down_DEGs:',nrow(sig_RUV.down)))+
        scale_color_manual(values=c("#e41a1c","#999999"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

DESeq_RUV_tpm_mean$compare <- ifelse(DESeq_RUV_tpm_mean$Row.names %in% DESeq_DEGs_nonRUV, "DESeq",
                                     ifelse(DESeq_RUV_tpm_mean$Row.names %in% DESeq_RUV_DEGs_nonDEseq, "DESeq_RUV",
                                            ifelse(DESeq_RUV_tpm_mean$Row.names %in% DEGs_overlap,"overlap_DEGs","none")))
DESeq_RUV_tpm_mean$compare <- factor(DESeq_RUV_tpm_mean$compare,levels=c('DESeq','DESeq_RUV','overlap_DEGs','none'))

pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_non_overlapDEG_tpm_mean_adjp0.01.pdf'),width=7,height=7, useDingbats = F)
print(ggplot(DESeq_RUV_tpm_mean)+
        geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=compare),size=1.5,alpha=0.8)+
        geom_text(data=DESeq_RUV_tpm_mean[DESeq_RUV_tpm_mean$Row.names %in% DEGs_nonoverlap,],
                   aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
                   size=3, check_overlap = TRUE)+
        xlab(paste0('log2(',sample1,'+1)'))+
        ylab(paste0('log2(',sample2,'+1)'))+
        labs(title = 'TPM: DESeq vs DESeq+RUV',
             subtitle = paste0('DEG_specific:',length(DESeq_DEGs_nonRUV),'  ','DEG_RUV_specific:',length(DESeq_RUV_DEGs_nonDEseq)))+
        scale_color_manual(values=c("#e41a1c","#4daf4a","#377eb8","lightgrey"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              plot.subtitle = element_text(face = "bold",
                                           size = rel(1.2), hjust = 0.5),
              text = element_text(),
              panel.background = element_blank(),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_pca.pdf'),width=7,height=7, useDingbats = F)
print(pcaplot(vsd.RUV,intgroup=c('Cell'),ntop=1000,
#print(pcaplot(vsd.RUV,intgroup=c('Type'),ntop=1000,
              text_labels = FALSE,
              pcX=1, pcY=2,title='PCA: DESeq + RUV', ellipse = TRUE)+
        geom_text_repel(aes(label=paste(vsd.RUV$Cell,vsd.RUV$Subject,sep='.')))+
        coord_fixed(ratio=0.8)+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.8, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

vsd.RUV.counts <- tidyr::gather(as.data.frame(assay(vsd.RUV)))
names(vsd.RUV.counts) <- c("Sample", "vstExpression")
vsd.RUV.counts$Sample <- factor(vsd.RUV.counts$Sample, levels = colnames(vsd.RUV))
des.RUV <- colData(vsd.RUV)
vsd.RUV.counts <- as.data.frame(merge(vsd.RUV.counts, des, by.x="Sample", by.y="row.names"))

pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_boxplot.pdf'),width=5,height=5, useDingbats = F)
print(ggplot(vsd.RUV.counts, aes_string(x = "Sample", y = "vstExpression")) + 
        #geom_boxplot(aes_string(color= "Type", fill = "Type"), 
        geom_boxplot(aes_string(color= "Cell", fill = "Cell"), 
                     alpha = 0.5)+
        ggtitle('Boxplot: DESeq + RUV')+
        scale_fill_manual(values=c("#fc8d59", "#99d594"))+
        scale_color_manual(values=c("#fc8d59", "#99d594"))+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text.x = element_text(angle=45,hjust=1),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.7, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
dev.off()

DESeq_RUV.sig.genes <- res.symbol[abs(res.RUV.symbol$log2FoldChange) > cutoff_fc & res.RUV.symbol$padj < cutoff_adjP,]

#ggthemr_reset()
p.volcano <- EnhancedVolcano(res.RUV.symbol,
                             lab = res.RUV.symbol$gene_name,
                             x = "log2FoldChange",
                             y = "padj",
                             xlab = bquote(~Log[2]~ "fold change"),
                             ylab = bquote(~-Log[10]~adjusted~italic(P)),
                             pCutoff = cutoff_adjP,
                             FCcutoff = cutoff_fc,
                             #xlim=c(-6,6),
                             transcriptLabSize = 3.0,
                             colAlpha = 1,
                             #title = paste0("Method:DESeq2",'\n','cutoff: ','log2FC, ',cutoff_fc,' adjP, ',cutoff_adjP),
                             legend=c("NS","Log2 FC","Adjusted p-value",
                                      "Adjusted p-value & Log2 FC"),
                             legendPosition = "bottom",
                             legendLabSize = 12,
                             legendIconSize = 3.0)+
  geom_text(aes(label=sample2, x= max(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.symbol$padj))/3),size=5,color='#d95f02',alpha=0.5)+
  geom_text(aes(label=sample1, x= min(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.symbol$padj))/3),size=5,color='#d95f02',alpha=0.5)
pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_volcano.pdf'),width=8,height=8, useDingbats = F)
print(grid.arrange(p.volcano, top=textGrob(paste0("Method: DESeq2+RUV, up_DEGs: ",nrow(sig_RUV.up),", down_DEGs:",nrow(sig_RUV.down),'\n','cutoff:  ','|log2FC|>',cutoff_fc,', adjP<',cutoff_adjP), gp=gpar(fontsize=15,fontface="bold"))))
dev.off()

#2*max(-log10(res.RUV.symbol$padj))/3
#2*max(-log10(res.RUV.symbol$padj))/3

## enrichment
DESeq_RUV_result_enrich_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/enrichment/")
dir.create(DESeq_RUV_result_enrich_folder, showWarnings = FALSE)
DEGs_RUV_sig_all <- res.RUV.sig$entrez
DEGs_RUV_sig_up <- res.RUV.sig[res.RUV.sig$log2FoldChange>0,]$entrez
DEGs_RUV_sig_down <- res.RUV.sig[res.RUV.sig$log2FoldChange<0,]$entrez
DEGs_RUV_df <- data.frame(Entrez=res.RUV.sig$entrez,FC=res.RUV.sig$log2FoldChange)
DEGs_RUV_df$group <- ifelse(DEGs_RUV_df$FC>0,"upregulated","downregulated")

# enrich KEGG
DESeq_RUV_result_enrich_KEGG_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/enrichment/KEGG")
dir.create(DESeq_RUV_result_enrich_KEGG_folder , showWarnings = FALSE)
DESeq_RUV_enrich_KEGG <- compareCluster(Entrez~group, data=DEGs_RUV_df, fun="enrichKEGG",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_RUV_enrich_KEGG)),as.data.frame(DESeq_RUV_enrich_KEGG)),
            file=paste0(DESeq_RUV_result_enrich_KEGG_folder,"/",sample1,'_',sample2,'_enrich_KEGG.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

pdf(file=paste0(DESeq_RUV_result_enrich_KEGG_folder,"/",sample1,'_',sample2,'_enrich_KEGG_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_RUV_enrich_KEGG,showCategory = 30)+
  labs(title='KEGG enrichment anlaysis',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_RUV_enrich_KEGG))))
dev.off()

# enrich GO
DESeq_RUV_result_enrich_GO_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/enrichment/GO")
dir.create(DESeq_RUV_result_enrich_GO_folder , showWarnings = FALSE)

DESeq_RUV_enrich_GO_BP <- compareCluster(Entrez~group, data=DEGs_RUV_df, 'org.Hs.eg.db',ont = "BP",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_RUV_enrich_GO_BP)),as.data.frame(DESeq_RUV_enrich_GO_BP)),
            file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_BP.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_BP_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_RUV_enrich_GO_BP,showCategory = 30)+
  labs(title='GO enrichment anlaysis, BP',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_RUV_enrich_GO_BP))))
dev.off()

DESeq_RUV_enrich_GO_MF <- compareCluster(Entrez~group, data=DEGs_RUV_df, 'org.Hs.eg.db',ont = "MF",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_RUV_enrich_GO_MF)),as.data.frame(DESeq_RUV_enrich_GO_MF)),
            file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_MF.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_MF_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_RUV_enrich_GO_MF,showCategory = 30)+
  labs(title='GO enrichment anlaysis, MF',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_RUV_enrich_GO_MF))))
dev.off()

DESeq_RUV_enrich_GO_CC <- compareCluster(Entrez~group, data=DEGs_RUV_df, 'org.Hs.eg.db',ont = "CC",fun="enrichGO",qvalueCutoff = 0.05,pvalueCutoff = 0.05)
write.table(data.frame("number"=rownames(as.data.frame(DESeq_RUV_enrich_GO_CC)),as.data.frame(DESeq_RUV_enrich_GO_CC)),
            file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_CC.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
pdf(file=paste0(DESeq_RUV_result_enrich_GO_folder,"/",sample1,'_',sample2,'_enrich_GO_CC_top30.pdf'),width=10,height=10, useDingbats = F)
dotplot(DESeq_RUV_enrich_GO_CC,showCategory = 30)+
  labs(title='GO enrichment anlaysis, CC',
       subtitle = paste0('total number:',nrow(as.data.frame(DESeq_RUV_enrich_GO_CC))))
dev.off()

## GAGE
DESeq_RUV_result_gage_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/gage/")
dir.create(DESeq_RUV_result_gage_folder, showWarnings = FALSE)
DESeq_RUV_foldchanges = res.RUV.Ordered$log2FoldChange
names(DESeq_RUV_foldchanges) = res.RUV.Ordered$entrez

# GAGE KEGG
DESeq_RUV_result_gage_KEGG_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/gage/KEGG/")
dir.create(DESeq_RUV_result_gage_KEGG_folder, showWarnings = FALSE)

DESeq_RUV_keggres = gage(DESeq_RUV_foldchanges, gsets=kegg.sets.hs, same.dir=FALSE)
write.table(data.frame("keggid"=rownames(DESeq_RUV_keggres$greater),DESeq_RUV_keggres$greater),
            file=paste0(DESeq_RUV_result_gage_KEGG_folder,"/",sample1,'_',sample2,'_gage_KEGG.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_RUV_keggres_select <- DESeq_RUV_keggres$greater[, "q.val"] < 0.05 & !is.na(DESeq_RUV_keggres$greater[, "q.val"])
DESeq_RUV_keggres_sigids <- rownames(DESeq_RUV_keggres$greater)[DESeq_RUV_keggres_select]
DESeq_RUV_keggres_sigids <- substr(DESeq_RUV_keggres_sigids, 1, 8)

setwd(DESeq_RUV_result_gage_KEGG_folder)
dir.create(paste('/Users/wen.zhong/Work/BloodAtlas/result',Sys.Date(),paste0(sample1,'_',sample2),'DESeq_RUV/gage/KEGG/','tmp',sep='/'), showWarnings = FALSE)
for(pwid in DESeq_RUV_keggres_sigids){
  pathview(gene.data=DESeq_RUV_foldchanges, 
           pathway.id=pwid, 
           kegg.native = TRUE,
           out.suffix=paste0('DESeq_',sample1,'_',sample2),
           kegg.dir = paste('/Users/wen.zhong/Work/BloodAtlas/result',Sys.Date(),paste0(sample1,'_',sample2),'DESeq_RUV/gage/KEGG/','tmp',sep='/'),
           species="hsa")}
setwd('/Users/wen.zhong/Work/BloodAtlas/')

# GAGE GO
DESeq_RUV_result_gage_GO_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/gage/GO/")
dir.create(DESeq_RUV_result_gage_GO_folder, showWarnings = FALSE)

DESeq_RUV_GO_MF <- gage(DESeq_RUV_foldchanges, gsets = go.sets.hs[go.subs.hs$MF],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_MF$greater),DESeq_RUV_GO_MF$greater),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_MF_greater.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_MF$less),DESeq_RUV_GO_MF$less),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_MF_less.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_RUV_GO_BP <- gage(DESeq_RUV_foldchanges, gsets = go.sets.hs[go.subs.hs$BP],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_BP$greater),DESeq_RUV_GO_BP$greater),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_BP_greater.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_BP$less),DESeq_RUV_GO_BP$less),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_BP_less.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

DESeq_RUV_GO_CC <- gage(DESeq_RUV_foldchanges, gsets = go.sets.hs[go.subs.hs$CC],same.dir=TRUE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_CC$greater),DESeq_RUV_GO_CC$greater),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_CC_greater.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
write.table(data.frame("goid"=rownames(DESeq_RUV_GO_CC$less),DESeq_RUV_GO_CC$less),
            file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO_CC_less.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)

#cutoff_GO_adjp <- 0.1
#cutoff_GO_p <- 0.01

DESeq_RUV_GO_select_BP_great <- as.data.frame(subset(DESeq_RUV_GO_BP$greater,DESeq_RUV_GO_BP$greater[,'p.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_BP_great$type <- rep('up',nrow(DESeq_RUV_GO_select_BP_great))
DESeq_RUV_GO_select_BP_great$category <- rep('BP',nrow(DESeq_RUV_GO_select_BP_great))
DESeq_RUV_GO_select_BP_less <- as.data.frame(subset(DESeq_RUV_GO_BP$less,DESeq_RUV_GO_BP$less[,'p.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_BP_less$type <- rep('down',nrow(DESeq_RUV_GO_select_BP_less))
DESeq_RUV_GO_select_BP_less$category <- rep('BP',nrow(DESeq_RUV_GO_select_BP_less))

DESeq_RUV_GO_select_MF_great <- as.data.frame(subset(DESeq_RUV_GO_MF$greater,DESeq_RUV_GO_MF$greater[,'p.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_MF_great$type <- rep('up',nrow(DESeq_RUV_GO_select_MF_great))
DESeq_RUV_GO_select_MF_great$category <- rep('MF',nrow(DESeq_RUV_GO_select_MF_great))
DESeq_RUV_GO_select_MF_less <- as.data.frame(subset(DESeq_RUV_GO_MF$less,DESeq_RUV_GO_MF$less[,'p.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_MF_less$type <- rep('down',nrow(DESeq_RUV_GO_select_MF_less))
DESeq_RUV_GO_select_MF_less$category <- rep('MF',nrow(DESeq_RUV_GO_select_MF_less))

DESeq_RUV_GO_select_CC_great <- as.data.frame(subset(DESeq_RUV_GO_CC$greater,DESeq_RUV_GO_CC$greater[,'p.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_CC_great$type <- rep('up',nrow(DESeq_RUV_GO_select_CC_great))
DESeq_RUV_GO_select_CC_great$category <- rep('CC',nrow(DESeq_RUV_GO_select_CC_great))
DESeq_RUV_GO_select_CC_less <- as.data.frame(subset(DESeq_RUV_GO_CC$less,DESeq_RUV_GO_CC$less[,'q.val'] < cutoff_GO_p))
DESeq_RUV_GO_select_CC_less$type <- rep('down',nrow(DESeq_RUV_GO_select_CC_less))
DESeq_RUV_GO_select_CC_less$category <- rep('CC',nrow(DESeq_RUV_GO_select_CC_less))
DESeq_RUV_GO_select <- rbind(DESeq_RUV_GO_select_BP_great,
                             DESeq_RUV_GO_select_BP_less,
                             DESeq_RUV_GO_select_MF_great,
                             DESeq_RUV_GO_select_MF_less,
                             DESeq_RUV_GO_select_CC_great,
                             DESeq_RUV_GO_select_CC_less)
DESeq_RUV_GO_select$GO <- str_split_fixed(rownames(DESeq_RUV_GO_select),' ',2)[,2]

pdf(file=paste0(DESeq_RUV_result_gage_GO_folder,"/",sample1,'_',sample2,'_gage_GO.pdf'),width=10,height=10, useDingbats = F)
# ggthemr('light')
print(ggplot(DESeq_RUV_GO_select)+
        geom_point(aes(x=stat.mean, y=-log10(p.val),size=-log10(p.val), fill=category, color=category))+
        geom_text_repel(mapping = aes(x=stat.mean, y=-log10(p.val),label=ifelse(q.val<0.05,GO,''),color=category),size=3)+
        geom_vline(xintercept = 0, color='black',linetype="dotted",size=1)+
        ggtitle(paste0('GO: DESeq+RUV','\n','cutoff: p<',cutoff_adjP,'(dots), adjp<0.05(text)'))+
        geom_text(aes(label=sample2, x= max(stat.mean)/2,y=5*max(-log10(p.val))/6),size=5,color='#984ea3')+
        geom_text(aes(label=sample1, x= min(stat.mean)/2,y=5*max(-log10(p.val))/6),size=5,color='#984ea3')+
        theme(plot.title = element_text(face = "bold",
                                        size = rel(1.5), hjust = 0.5),
              text = element_text(),
              #panel.grid.major = element_blank(),
              #panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              #panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_blank(),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text.x = element_text(angle=45,hjust=1),
              axis.line = element_line(colour="black",size=0.5),
              axis.ticks = element_line(),
              panel.grid.major = element_line(color='grey80',linetype = "dashed"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size=10),
              legend.key.size= unit(0.7, "cm"),
              legend.title = element_text(face="italic",size=12),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")))
#ggthemr_reset()
dev.off()

## compare RUV and without RUV
# venn.diagram(x=list(DESeq2=rownames(res.sig),
#                   DESeq2_RUV=rownames(res.RUV.sig)),
#              file=paste("./result/plot/DESeq_RUV_compare/",paste(sample1,sample2,sep='_'),'.tiff',sep=''),
#              fill=colors)
# 

sum.out <- data.frame(sample1=sample1,
                      sample2=sample2,
                      DESeq_padj_0.05=nrow(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]),
                      DESeq_RUV_padj_0.05=nrow(res.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
                      overlap_padj_0.05=length(intersect(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
                                                         res.RUV.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
                      DESeq_padj_0.01=nrow(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]),
                      DESeq_RUV_padj_0.01=nrow(res.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
                      overlap_padj_0.01=length(intersect(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
                                                         res.RUV.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
                      DESeq_padj_0.005=nrow(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]),
                      DESeq_RUV_padj_0.005=nrow(res.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
                      overlap_padj_0.005=length(intersect(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
                      DESeq_padj_0.001=nrow(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]),
                      DESeq_RUV_padj_0.001=nrow(res.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
                      overlap_padj_0.001=length(intersect(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)))
write.table(sum.out,
            file=paste0(sample1_sample2_folder,"/",sample1,'_',sample2,'_DEGs_summary.txt'),
            sep='\t',
            row.names=FALSE,
            quote=FALSE)
