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

setwd('/Users/max.karlsson/Documents/Scilifelab/Projects/HPA-classification/')
setwd('/Users/wen.zhong/Work/BloodAtlas/')
result_folder <- paste0("./result/", Sys.Date())
dir.create(result_folder, showWarnings = FALSE)

gene.info <- read.delim("../database/ensembl92/geneinfo_92.txt",header = TRUE, sep='\t',check.names = T)
#protein.class <- read.table(paste(workpath,'RNAseq','annotation','gene.classes.txt',sep='/'),header = TRUE, sep='\t')
#rownames(protein.class) <- protein.class$rna.genes

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

## filter PBMCs
blood.Set <- blood.Set[,!blood.Set$Cat2 == 'total_PBMC']
blood.tpm.Set <- blood.tpm.Set[,!blood.tpm.Set$Cat2 == 'total_PBMC']

## sample select
for(i in unique(pData(blood.Set)$Cat2)){
  combns <- combn(unique(pData(blood.Set)[pData(blood.Set)$Cat2==i,]$Cell), 2)
  N <- seq_len(ncol(combns))
  for(cols in 1:ncol(combns)){
    pair <- combns[, cols]
    sample1 <- as.character(pair[1])
    sample2 <- as.character(pair[2])  

    
    sample1_sample2_folder <- paste0(result_folder, "/", sample1,'_',sample2)
    dir.create(sample1_sample2_folder, showWarnings = FALSE)
    
    blood.Set.select <- blood.Set[,blood.Set$Cell == sample1|blood.Set$Cell == sample2]
    #blood.Set.select <- blood.Set[,blood.Set$Cat2 == sample1|blood.Set$Cat2 == sample2]
    filter <- apply(counts(blood.Set.select),1,function(x) length(x[x>5])>=2)
    genes.filter <- counts(blood.Set.select)[filter,]
    pheno.select <- pData(blood.Set.select)
    pheno.select$Cell <- factor(pheno.select$Cell)
    #pheno.select$Cat2 <- factor(pheno.select$Cat2)
    blood.Set.select <- newSeqExpressionSet(genes.filter,
                                            phenoData = pheno.select)
    idxbycell <- order(blood.Set.select$Cell)
    #idxbytype<- order(blood.Set.select$Cat2)
    blood.Set.select <- blood.Set.select[,idxbycell]
    #blood.Set.select <- blood.Set.select[,idxbytype]
    
    blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Cell == sample1|blood.tpm.Set$Cell == sample2]
    #blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Cat2 == sample1|blood.tpm.Set$Cat2 == sample2]
    pheno.tpm.select <- pData(blood.tpm.Set.select)
    pheno.tpm.select$Cell <- factor(pheno.tpm.select$Cell)
    #pheno.tpm.select$Cat2 <- factor(pheno.tpm.select$Cat2)
    blood.tpm.Set.select <- newSeqExpressionSet(counts(blood.tpm.Set.select),
                                                phenoData = pheno.tpm.select)
    
    ## PCA
    colors <- c("#f1a340","#998ec3")
    x <- blood.Set.select$Cell
    #x <- blood.Set.select$Cat2
    #plotRLE(blood.Set.select,outline=FALSE, ylim=c(-4, 4),col=colors[x])
    #plotPCA(blood.Set.select, col=colors[x], cex=1.2)
    
    ## RUV
    design <- model.matrix(~Cell, data=pData(blood.Set.select))
    #design <- model.matrix(~Cat2, data=pData(blood.Set.select))
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
    
    ## DESeq2 + RUV
    dds.RUV <- DESeqDataSetFromMatrix(countData = counts(blood.Set.select.RUV),
                                      colData = pData(blood.Set.select.RUV),
    #                                  design = ~ W_1+Cat2)
                                      design = ~ W_1+Cell)
    dds.RUV <- DESeq(dds.RUV)
    vsd.RUV <- vst(dds.RUV, blind=FALSE)
    res.RUV <- results(dds.RUV, contrast=c("Cell",sample2,sample1))
    #res.RUV <- results(dds.RUV, contrast=c("Cat2", as.character(sample2), as.character(sample1)))
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
    sig_RUV.up <- res.RUV.sig[res.RUV.sig$log2FoldChange>0,]
    sig_RUV.down <- res.RUV.sig[res.RUV.sig$log2FoldChange<0,]
    
    
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
    #DESeq_RUV_basemean <- sapply(levels(dds.RUV$Cat2), function(lvl) rowMeans(counts(dds.RUV,normalized=TRUE)[,dds.RUV$Cat2 == lvl]))
    colnames(DESeq_RUV_basemean) <- c('sample2','sample1') 
    DESeq_RUV_basemean <- merge(DESeq_RUV_basemean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
    DESeq_RUV_basemean$sig <- ifelse(DESeq_RUV_basemean$Row.names %in% res.RUV.sig$Row.names,'DEG','no')
    #pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_basemean.pdf'),width=7,height=7, useDingbats = F)
    g1 <- ggplot(DESeq_RUV_basemean)+
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
                  strip.text = element_text(face="bold"),
                  aspect.ratio = 1)
    
    
    
    DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cell), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cell == lvl]))
    #DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cat2), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cat2 == lvl]))
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
    
    
    pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_pca.pdf'),width=7,height=7, useDingbats = F)
    print(pcaplot(vsd.RUV,intgroup=c('Cell'),ntop=1000,
    #print(pcaplot(vsd.RUV,intgroup=c('Cat2'),ntop=1000,
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
    vsd.RUV.counts <- as.data.frame(merge(vsd.RUV.counts, des.RUV, by.x="Sample", by.y="row.names"))
    
    pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_boxplot.pdf'),width=5,height=5, useDingbats = F)
    print(ggplot(vsd.RUV.counts, aes_string(x = "Sample", y = "vstExpression")) + 
            geom_boxplot(aes_string(color= "Cat2", fill = "Cat2"), 
                         #geom_boxplot(aes_string(color= "Cell", fill = "Cell"), 
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
    
    DESeq_RUV.sig.genes <- res.RUV.symbol[abs(res.RUV.symbol$log2FoldChange) > cutoff_fc & res.RUV.symbol$padj < cutoff_adjP,]
    
    
    
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
                                 axisLabSize = 10,
                                 gridlines.major = F,
                                 gridlines.minor = F,
                                 transcriptLabSize = 3.0,
                                 colAlpha = 1,
                                 col = c("grey30", "forestgreen", "royalblue", "red2"),
                                 title = paste0(sample1," vs ",sample2, "\n","up_DEGs: ",nrow(sig_RUV.up),", down_DEGs:",nrow(sig_RUV.down)),
                                 #titleLabSize = 10,
                                 legend=c("NS","Log2 FC","Adjusted p-value",
                                          "Adjusted p-value & Log2 FC"),
                                 legendPosition = "right",
                                 legendLabSize = 10,
                                 legendIconSize = 3.0)+
      # theme_option_2+
      theme(panel.border=element_blank(),
            plot.margin=unit(c(0,0,0,0), "lines"),
            plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            # legend.position = 'right',
            #axis.line = element_blank(),
            aspect.ratio = 0.8)
    #geom_text(label=sample2, x= max(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.RUV.symbol$padj)[!is.infinite(-log10(res.RUV.symbol$padj))])/3,size=2,color='black',alpha=0.5)+
    #geom_text(label=sample1, x= min(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.RUV.symbol$padj)[!is.infinite(-log10(res.RUV.symbol$padj))])/3,size=2,color='black',alpha=0.5)
    pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_volcano_basemean.pdf'),width=23,height=8, useDingbats = F)
    grid.arrange(g1, p.volcano, ncol = 2)
    dev.off()
    
    ## compare RUV and without RUV
    # venn.diagram(x=list(DESeq2=rownames(res.sig),
    #                   DESeq2_RUV=rownames(res.RUV.sig)),
    #              file=paste("./result/plot/DESeq_RUV_compare/",paste(sample1,sample2,sep='_'),'.tiff',sep=''),
    #              fill=colors)
    # 
    
    # sum.out <- data.frame(sample1=sample1,
    #                       sample2=sample2,
    #                       DESeq_padj_0.05=nrow(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]),
    #                       DESeq_RUV_padj_0.05=nrow(res.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
    #                       overlap_padj_0.05=length(intersect(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
    #                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
    #                       DESeq_padj_0.01=nrow(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]),
    #                       DESeq_RUV_padj_0.01=nrow(res.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
    #                       overlap_padj_0.01=length(intersect(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
    #                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
    #                       DESeq_padj_0.005=nrow(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]),
    #                       DESeq_RUV_padj_0.005=nrow(res.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
    #                       overlap_padj_0.005=length(intersect(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
    #                                                           res.RUV.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
    #                       DESeq_padj_0.001=nrow(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]),
    #                       DESeq_RUV_padj_0.001=nrow(res.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
    #                       overlap_padj_0.001=length(intersect(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
    #                                                           res.RUV.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)))
    # write.table(sum.out,
    #             file=paste0(sample1_sample2_folder,"/",sample1,'_',sample2,'_DEGs_summary.txt'),
    #             sep='\t',
    #             row.names=FALSE,
    #             quote=FALSE)
  }
}



# 
# 
# 
# combns <- combn(unique(pData(blood.Set)$Cell), 2)
# N <- seq_len(ncol(combns))
# for(cols in 1:ncol(combns)){
#   pair <- combns[, cols]
#   sample1 <- as.character(pair[1])
#   sample2 <- as.character(pair[2])  
#   
#   # combn(levels(pdata$Cat2), 2, function(x){
#   #   sample1 <- x[1]
#   #   sample2 <- x[2]
#   # sample1 <- 'tregs'
#   # sample2 <- 'memory_cd8'
#   #sample1 <- 'B_cells'
#   #sample2 <- 'T_cells'
#   
#   sample1_sample2_folder <- paste0(result_folder, "/", sample1,'_',sample2)
#   dir.create(sample1_sample2_folder, showWarnings = FALSE)
#   
#   #blood.Set.select <- blood.Set[,blood.Set$Cell == sample1|blood.Set$Cell == sample2]
#   blood.Set.select <- blood.Set[,blood.Set$Cat2 == sample1|blood.Set$Cat2 == sample2]
#   filter <- apply(counts(blood.Set.select),1,function(x) length(x[x>5])>=2)
#   genes.filter <- counts(blood.Set.select)[filter,]
#   pheno.select <- pData(blood.Set.select)
#   #pheno.select$Cell <- factor(pheno.select$Cell)
#   pheno.select$Cat2 <- factor(pheno.select$Cat2)
#   blood.Set.select <- newSeqExpressionSet(genes.filter,
#                                           phenoData = pheno.select)
#   #idxbycell <- order(blood.Set.select$Cell)
#   idxbytype<- order(blood.Set.select$Cat2)
#   #blood.Set.select <- blood.Set.select[,idxbycell]
#   blood.Set.select <- blood.Set.select[,idxbytype]
#   
#   #blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Cell == sample1|blood.tpm.Set$Cell == sample2]
#   blood.tpm.Set.select <- blood.tpm.Set[,blood.tpm.Set$Cat2 == sample1|blood.tpm.Set$Cat2 == sample2]
#   pheno.tpm.select <- pData(blood.tpm.Set.select)
#   #pheno.tpm.select$Cell <- factor(pheno.tpm.select$Cell)
#   pheno.tpm.select$Cat2 <- factor(pheno.tpm.select$Cat2)
#   blood.tpm.Set.select <- newSeqExpressionSet(counts(blood.tpm.Set.select),
#                                               phenoData = pheno.tpm.select)
#   
#   ## PCA
#   colors <- c("#f1a340","#998ec3")
#   #x <- blood.Set.select$Cell
#   x <- blood.Set.select$Cat2
#   #plotRLE(blood.Set.select,outline=FALSE, ylim=c(-4, 4),col=colors[x])
#   #plotPCA(blood.Set.select, col=colors[x], cex=1.2)
#   
#   ## RUV
#   #design <- model.matrix(~Cell, data=pData(blood.Set.select))
#   design <- model.matrix(~Cat2, data=pData(blood.Set.select))
#   y <- DGEList(counts=counts(blood.Set.select), group=x)
#   y <- calcNormFactors(y, method="upperquartile")
#   y <- estimateGLMCommonDisp(y, design)
#   y <- estimateGLMTagwiseDisp(y, design)
#   fit <- glmFit(y, design)
#   lrt <- glmLRT(fit, coef=2)
#   top <- topTags(lrt, n=nrow(blood.Set.select))$table
#   #empirical <- rownames(top[top$FDR>0.5,])
#   empirical <- rownames(blood.Set.select)[which(!(rownames(blood.Set.select) %in% rownames(top)[1:5000]))]
#   
#   blood.Set.select.RUV <- RUVg(blood.Set.select, empirical, k=1)
#   #plotRLE(blood.Set.select.RUV, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#   #plotPCA(blood.Set.select.RUV, col=colors[x], cex=1.2)
#   
#   ## set cutoff
#   cutoff_fc <- 1
#   cutoff_adjP <- 0.01
#   
#   ## DESeq2 + RUV
#   dds.RUV <- DESeqDataSetFromMatrix(countData = counts(blood.Set.select.RUV),
#                                     colData = pData(blood.Set.select.RUV),
#                                     design = ~ W_1+Cat2)
#                                     #design = ~ W_1+Cell)
#   dds.RUV <- DESeq(dds.RUV)
#   vsd.RUV <- vst(dds.RUV, blind=FALSE)
#   #res.RUV <- results(dds.RUV, contrast=c("Cell",sample2,sample1))
#   res.RUV <- results(dds.RUV, contrast=c("Cat2", as.character(sample2), as.character(sample1)))
#   res.RUV <- na.omit(res.RUV)
#   res.RUV.symbol <- merge(data.frame(res.RUV), gene.info, by.x="row.names", by.y="ensg_id",all.x=TRUE)
#   res.RUV.Ordered <- res.RUV.symbol[order(res.RUV.symbol$pvalue),]
#   res.RUV.Ordered$entrez <- mapIds(org.Hs.eg.db,
#                                    keys=res.RUV.Ordered$Row.names, 
#                                    column="ENTREZID",
#                                    keytype="ENSEMBL",
#                                    multiVals="first")
#   res.RUV.sig <- res.RUV.Ordered[res.RUV.Ordered$padj<cutoff_adjP & abs(res.RUV.Ordered$log2FoldChange)>cutoff_fc,]
#   sig_DEGs_RUV <- res.RUV.sig[order(abs(res.RUV.sig$log2FoldChange),decreasing = TRUE),]$Row.names
#   sig_RUV.up <- res.RUV.sig[res.RUV.sig$log2FoldChange>0,]
#   sig_RUV.down <- res.RUV.sig[res.RUV.sig$log2FoldChange<0,]
#   
#   
#   DESeq_RUV_result_folder <- paste0(sample1_sample2_folder,"/DESeq_RUV/")
#   dir.create(DESeq_RUV_result_folder, showWarnings = FALSE)
#   write.table(data.frame("ensgid"=rownames(res.RUV.Ordered),res.RUV.Ordered),
#               file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_pvalue.txt'),
#               sep='\t',
#               row.names=FALSE,
#               quote=FALSE)
#   write.table(data.frame("ensgid"=rownames(res.RUV.sig),res.RUV.sig),
#               file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_sig_pvalue.txt'),
#               sep='\t',
#               row.names=FALSE,
#               quote=FALSE)
#   
#   #DESeq_RUV_basemean <- sapply(levels(dds.RUV$Cell), function(lvl) rowMeans(counts(dds.RUV,normalized=TRUE)[,dds.RUV$Cell == lvl]))
#   DESeq_RUV_basemean <- sapply(levels(dds.RUV$Cat2), function(lvl) rowMeans(counts(dds.RUV,normalized=TRUE)[,dds.RUV$Cat2 == lvl]))
#   colnames(DESeq_RUV_basemean) <- c('sample2','sample1') 
#   DESeq_RUV_basemean <- merge(DESeq_RUV_basemean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
#   DESeq_RUV_basemean$sig <- ifelse(DESeq_RUV_basemean$Row.names %in% res.RUV.sig$Row.names,'DEG','no')
#   pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_basemean.pdf'),width=7,height=7, useDingbats = F)
#   print(ggplot(DESeq_RUV_basemean)+
#           geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
#           geom_text(data=DESeq_RUV_basemean[DESeq_RUV_basemean$Row.names %in% sig_DEGs_RUV,],
#                     aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
#                     size=2.8, check_overlap = TRUE)+
#           xlab(paste0('log2(',sample1,'+1)'))+
#           ylab(paste0('log2(',sample2,'+1)'))+
#           labs(title = 'BaseMean: DESeq+RUV normalized counts',
#                subtitle = paste0('up_DEGs:',nrow(sig_RUV.up),'  ','down_DEGs:',nrow(sig_RUV.down)))+
#           scale_color_manual(values=c("#e41a1c","#999999"))+
#           theme(plot.title = element_text(face = "bold",
#                                           size = rel(1.5), hjust = 0.5),
#                 plot.subtitle = element_text(face = "bold",
#                                              size = rel(1.2), hjust = 0.5),
#                 text = element_text(),
#                 panel.background = element_blank(),
#                 plot.background = element_rect(colour = NA),
#                 panel.border = element_blank(),
#                 axis.title = element_text(face = "bold",size = rel(1)),
#                 axis.title.y = element_text(angle=90,vjust =2),
#                 axis.title.x = element_text(vjust = -0.2),
#                 axis.text = element_text(),
#                 axis.line = element_line(colour="black",size=0.5),
#                 axis.ticks = element_line(),
#                 panel.grid.major = element_line(color='grey80',linetype = "dashed"),
#                 panel.grid.minor = element_blank(),
#                 legend.key = element_rect(colour = NA),
#                 legend.position = "bottom",
#                 legend.direction = "horizontal",
#                 legend.text = element_text(size=10),
#                 legend.key.size= unit(0.8, "cm"),
#                 legend.title = element_text(face="italic",size=12),
#                 plot.margin=unit(c(10,5,5,5),"mm"),
#                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                 strip.text = element_text(face="bold")))
#   dev.off()
#   
# 
#   
#   #DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cell), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cell == lvl]))
#   DESeq_RUV_tpm_mean <- sapply(levels(blood.tpm.Set.select$Cat2), function(lvl) rowMeans(counts(blood.tpm.Set.select)[,blood.tpm.Set.select$Cat2 == lvl]))
#   colnames(DESeq_RUV_tpm_mean) <- c('sample2','sample1') 
#   DESeq_RUV_tpm_mean <- merge(DESeq_RUV_tpm_mean, gene.info, by.x="row.names", by.y="ensg_id", all.x=TRUE)
#   DESeq_RUV_tpm_mean$sig <- ifelse(DESeq_RUV_tpm_mean$Row.names %in% res.RUV.sig$Row.names,'DEG','no')
#   pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_tpm_mean.pdf'),width=7,height=7, useDingbats = F)
#   print(ggplot(DESeq_RUV_tpm_mean)+
#           geom_point(aes(x=log2(sample1+1),y=log2(sample2+1),color=sig),size=1.5)+
#           geom_text(data=DESeq_RUV_tpm_mean[DESeq_RUV_tpm_mean$Row.names %in% sig_DEGs_RUV,],
#                     aes(x=log2(sample1+1),y=log2(sample2+1),label=gene_name),
#                     size=2.8, check_overlap = TRUE)+
#           xlab(paste0('log2(',sample1,'+1)'))+
#           ylab(paste0('log2(',sample2,'+1)'))+
#           labs(title = 'TPM: DESeq+RUV',
#                subtitle = paste0('up_DEGs:',nrow(sig_RUV.up),'  ','down_DEGs:',nrow(sig_RUV.down)))+
#           scale_color_manual(values=c("#e41a1c","#999999"))+
#           theme(plot.title = element_text(face = "bold",
#                                           size = rel(1.5), hjust = 0.5),
#                 plot.subtitle = element_text(face = "bold",
#                                              size = rel(1.2), hjust = 0.5),
#                 text = element_text(),
#                 panel.background = element_blank(),
#                 plot.background = element_rect(colour = NA),
#                 panel.border = element_blank(),
#                 axis.title = element_text(face = "bold",size = rel(1)),
#                 axis.title.y = element_text(angle=90,vjust =2),
#                 axis.title.x = element_text(vjust = -0.2),
#                 axis.text = element_text(),
#                 axis.line = element_line(colour="black",size=0.5),
#                 axis.ticks = element_line(),
#                 panel.grid.major = element_line(color='grey80',linetype = "dashed"),
#                 panel.grid.minor = element_blank(),
#                 legend.key = element_rect(colour = NA),
#                 legend.position = "bottom",
#                 legend.direction = "horizontal",
#                 legend.text = element_text(size=10),
#                 legend.key.size= unit(0.8, "cm"),
#                 legend.title = element_text(face="italic",size=12),
#                 plot.margin=unit(c(10,5,5,5),"mm"),
#                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                 strip.text = element_text(face="bold")))
#   dev.off()
#   
#   
#   pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_pca.pdf'),width=7,height=7, useDingbats = F)
#   #print(pcaplot(vsd.RUV,intgroup=c('Cell'),ntop=1000,
#   print(pcaplot(vsd.RUV,intgroup=c('Cat2'),ntop=1000,
#                 text_labels = FALSE,
#                 pcX=1, pcY=2,title='PCA: DESeq + RUV', ellipse = TRUE)+
#           geom_text_repel(aes(label=paste(vsd.RUV$Cell,vsd.RUV$Subject,sep='.')))+
#           coord_fixed(ratio=0.8)+
#           theme(plot.title = element_text(face = "bold",
#                                           size = rel(1.5), hjust = 0.5),
#                 text = element_text(),
#                 #panel.grid.major = element_blank(),
#                 #panel.grid.minor = element_blank(),
#                 panel.background = element_blank(),
#                 #panel.background = element_rect(colour = NA),
#                 plot.background = element_rect(colour = NA),
#                 panel.border = element_blank(),
#                 axis.title = element_text(face = "bold",size = rel(1)),
#                 axis.title.y = element_text(angle=90,vjust =2),
#                 axis.title.x = element_text(vjust = -0.2),
#                 axis.text = element_text(),
#                 axis.line = element_line(colour="black",size=0.5),
#                 axis.ticks = element_line(),
#                 panel.grid.major = element_line(color='grey80',linetype = "dashed"),
#                 panel.grid.minor = element_blank(),
#                 legend.key = element_rect(colour = NA),
#                 legend.position = "bottom",
#                 legend.direction = "horizontal",
#                 legend.text = element_text(size=10),
#                 legend.key.size= unit(0.8, "cm"),
#                 legend.title = element_text(face="italic",size=12),
#                 plot.margin=unit(c(10,5,5,5),"mm"),
#                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                 strip.text = element_text(face="bold")))
#   dev.off()
#   
#   vsd.RUV.counts <- tidyr::gather(as.data.frame(assay(vsd.RUV)))
#   names(vsd.RUV.counts) <- c("Sample", "vstExpression")
#   vsd.RUV.counts$Sample <- factor(vsd.RUV.counts$Sample, levels = colnames(vsd.RUV))
#   des.RUV <- colData(vsd.RUV)
#   vsd.RUV.counts <- as.data.frame(merge(vsd.RUV.counts, des.RUV, by.x="Sample", by.y="row.names"))
#   
#   pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_vsd_boxplot.pdf'),width=5,height=5, useDingbats = F)
#   print(ggplot(vsd.RUV.counts, aes_string(x = "Sample", y = "vstExpression")) + 
#           geom_boxplot(aes_string(color= "Cat2", fill = "Cat2"), 
#           #geom_boxplot(aes_string(color= "Cell", fill = "Cell"), 
#                        alpha = 0.5)+
#           ggtitle('Boxplot: DESeq + RUV')+
#           scale_fill_manual(values=c("#fc8d59", "#99d594"))+
#           scale_color_manual(values=c("#fc8d59", "#99d594"))+
#           theme(plot.title = element_text(face = "bold",
#                                           size = rel(1.5), hjust = 0.5),
#                 text = element_text(),
#                 #panel.grid.major = element_blank(),
#                 #panel.grid.minor = element_blank(),
#                 panel.background = element_blank(),
#                 #panel.background = element_rect(colour = NA),
#                 plot.background = element_rect(colour = NA),
#                 panel.border = element_blank(),
#                 axis.title = element_text(face = "bold",size = rel(1)),
#                 axis.title.y = element_text(angle=90,vjust =2),
#                 axis.title.x = element_text(vjust = -0.2),
#                 axis.text.x = element_text(angle=45,hjust=1),
#                 axis.line = element_line(colour="black",size=0.5),
#                 axis.ticks = element_line(),
#                 panel.grid.major = element_line(color='grey80',linetype = "dashed"),
#                 panel.grid.minor = element_blank(),
#                 legend.key = element_rect(colour = NA),
#                 legend.position = "bottom",
#                 legend.direction = "horizontal",
#                 legend.text = element_text(size=10),
#                 legend.key.size= unit(0.7, "cm"),
#                 legend.title = element_text(face="italic",size=12),
#                 plot.margin=unit(c(10,5,5,5),"mm"),
#                 strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
#                 strip.text = element_text(face="bold")))
#   dev.off()
#   
#   DESeq_RUV.sig.genes <- res.RUV.symbol[abs(res.RUV.symbol$log2FoldChange) > cutoff_fc & res.RUV.symbol$padj < cutoff_adjP,]
#   
#   #ggthemr_reset()
#   p.volcano <- EnhancedVolcano(res.RUV.symbol,
#                                lab = res.RUV.symbol$gene_name,
#                                x = "log2FoldChange",
#                                y = "padj",
#                                xlab = bquote(~Log[2]~ "fold change"),
#                                ylab = bquote(~-Log[10]~adjusted~italic(P)),
#                                pCutoff = cutoff_adjP,
#                                FCcutoff = cutoff_fc,
#                                #xlim=c(-6,6),
#                                axisLabSize = 10,
#                                gridlines.major = F,
#                                gridlines.minor = F,
#                                transcriptLabSize = 3.0,
#                                colAlpha = 1,
#                                col = c("grey30", "forestgreen", "royalblue", "red2"),
#                                title = paste0(sample1," vs ",sample2, "\n","up_DEGs: ",nrow(sig_RUV.up),", down_DEGs:",nrow(sig_RUV.down)),
#                                #titleLabSize = 10,
#                                legend=c("NS","Log2 FC","Adjusted p-value",
#                                         "Adjusted p-value & Log2 FC"),
#                                legendPosition = "right",
#                                legendLabSize = 10,
#                                legendIconSize = 3.0)+
#    # theme_option_2+
#     theme(panel.border=element_blank(),
#           plot.margin=unit(c(0,0,0,0), "lines"),
#           plot.title = element_text(face = "bold",
#                                     size = rel(1.2), hjust = 0.5),
#          # legend.position = 'right',
#           #axis.line = element_blank(),
#           aspect.ratio = 0.8)
#     #geom_text(label=sample2, x= max(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.RUV.symbol$padj)[!is.infinite(-log10(res.RUV.symbol$padj))])/3,size=2,color='black',alpha=0.5)+
#     #geom_text(label=sample1, x= min(res.RUV.symbol$log2FoldChange)/2,y=2*max(-log10(res.RUV.symbol$padj)[!is.infinite(-log10(res.RUV.symbol$padj))])/3,size=2,color='black',alpha=0.5)
#   pdf(file=paste0(DESeq_RUV_result_folder,"/",sample1,'_',sample2,'_DESeq_RUV_volcano.pdf'),width=8,height=8, useDingbats = F)
#   print(p.volcano)
#   dev.off()
#   
#   ## compare RUV and without RUV
#   # venn.diagram(x=list(DESeq2=rownames(res.sig),
#   #                   DESeq2_RUV=rownames(res.RUV.sig)),
#   #              file=paste("./result/plot/DESeq_RUV_compare/",paste(sample1,sample2,sep='_'),'.tiff',sep=''),
#   #              fill=colors)
#   # 
#   
#   # sum.out <- data.frame(sample1=sample1,
#   #                       sample2=sample2,
#   #                       DESeq_padj_0.05=nrow(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]),
#   #                       DESeq_RUV_padj_0.05=nrow(res.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
#   #                       overlap_padj_0.05=length(intersect(res.Ordered[res.Ordered$padj<0.05 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
#   #                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.05 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
#   #                       DESeq_padj_0.01=nrow(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]),
#   #                       DESeq_RUV_padj_0.01=nrow(res.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
#   #                       overlap_padj_0.01=length(intersect(res.Ordered[res.Ordered$padj<0.01 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
#   #                                                          res.RUV.Ordered[res.RUV.Ordered$padj<0.01 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
#   #                       DESeq_padj_0.005=nrow(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]),
#   #                       DESeq_RUV_padj_0.005=nrow(res.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
#   #                       overlap_padj_0.005=length(intersect(res.Ordered[res.Ordered$padj<0.005 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
#   #                                                           res.RUV.Ordered[res.RUV.Ordered$padj<0.005 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)),
#   #                       DESeq_padj_0.001=nrow(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]),
#   #                       DESeq_RUV_padj_0.001=nrow(res.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]),
#   #                       overlap_padj_0.001=length(intersect(res.Ordered[res.Ordered$padj<0.001 & abs(res.Ordered$log2FoldChange)>1,]$Row.names,
#   #                                                           res.RUV.Ordered[res.RUV.Ordered$padj<0.001 & abs(res.RUV.Ordered$log2FoldChange)>1,]$Row.names)))
#   # write.table(sum.out,
#   #             file=paste0(sample1_sample2_folder,"/",sample1,'_',sample2,'_DEGs_summary.txt'),
#   #             sep='\t',
#   #             row.names=FALSE,
#   #             quote=FALSE)
# }
