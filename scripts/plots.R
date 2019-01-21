library('pcaMethods')
library('ggrepel')
library('gridExtra')
library('ClassDiscovery')
library('gplots')


make_gene_expression_barplot <- function(atlas.max.tb, maxEx_columns, content_column, content_color) {
  
  plot.data <- 
    atlas.max.tb %>% 
    select("content_column" = content_column, ensg_id, maxEx_columns) %>%
    gather(key = "Type", value = "Expression", -(1:2)) 
  
  genes <- unique(plot.data$ensg_id)
  for(i in seq(1, length(genes), 4)) {
    j = i + 3
    if(j > length(genes)) j = length(genes)
    
    
    plot.data %>%
      filter(ensg_id %in% genes[i:j]) %>%
      ggplot(aes(content_column, Expression, fill = content_column)) + 
      geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
      
      geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
      simple_theme+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
      scale_fill_manual(values = content_color)+
      facet_grid(ensg_id ~ Type, scales = "free") 
    ggsave(paste(result_folder, paste0(prefix, "gene expression barplot", i, "-", j, ".png"), sep = "/"), width = 20, height = 10)
  }
}

make_tissue_distributions_plot <- function(atlas.tb, Ex_column, content_column, und.lim, do.tissues = "all", outpath, prefix) {
  
  atlas.tb %>%
    mutate(tissue = eval(parse(text = content_column)),
           expression = eval(parse(text = Ex_column))) %>%
    filter(ifelse(rep(do.tissues[1] == "all", nrow(.)), T, tissue %in% do.tissues) & expression > 0) %>%
    mutate(tissue = factor(tissue, levels = with({group_by(., tissue) %>% 
        summarise(methods = paste(unique(method), collapse = " "))},
        tissue[order(methods)]))) %>%
        {ggplot(., aes(tissue, expression, fill = method, color = method))+
            stat_summary(fun.y = "median", fun.args = c("na.rm" = T), geom = "line", aes(group = method), size = 2, alpha = 0.2)+
            stat_summary(fun.y = "min", geom = "line", aes(tissue, eval(parse(text = Ex_column)), group = 1), inherit.aes = F, size = 1, alpha = 0.1)+
            stat_summary(fun.y = "max", geom = "line", aes(tissue, eval(parse(text = Ex_column)), group = 1), inherit.aes = F, size = 1, alpha = 0.1)+
            stat_summary(fun.y = "min", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            stat_summary(fun.y = "max", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            
            
            geom_violin(draw_quantiles = 0.5, alpha = 0.2, position = "identity")+
            geom_hline(yintercept = und.lim, linetype = "dashed")+
            simple_theme+
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
            theme(axis.title = element_blank())+
            scale_fill_manual(values = dataset.colors)+
            scale_color_manual(values = dataset.colors)} 
  
  ggsave(filename = paste(outpath, paste0(prefix, " tissue distributions.png"), sep = "/"), width = 10, height = 5, dpi = 300)
}


### 1. PCA plot
pca.cal <- function(tb.wide){
  pca <- 
    tb.wide %>%
    #filter row with NA
    {names <- rownames(.); 
    filt <- !apply(., 1, function(x){any(is.na(x))}); 
    filter(., filt) %>% mutate(ensg_id = names[filt])} %>%
    column_to_rownames("ensg_id") %>% 
    as.matrix() %>%
    `+`(1) %>%
    log10() %>%
    t() %>%
    pcaMethods::pca(center = T, nPcs = ncol(tb.wide))
  
  PCA.scores <- pcaMethods::scores(pca)
  PCA.loadings <- pcaMethods::loadings(pca)
  return(list(PCA.scores, PCA.loadings))
} 

make_PCA_plots <- function(scores, loadings, groups, groups.color, ellipse = T, outpath, prefix) {
  loadings %>% 
    ggplot(aes(PC1, PC2, label = labels))+
    geom_text(show.legend = F)+
    simple_theme+
    ggtitle("Loadings")
  ggsave(paste(outpath, paste0(prefix, '_PCA_loadings.pdf'),sep='/'), width=10, height=10)
  
  scores %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(sample = names,
                                                 groups = groups[match(names, names(groups))])} %>%
                                                 {means <- group_by(., groups) %>%
                                                   dplyr::summarise(mean.PC1 = mean(PC1),
                                                                    mean.PC2 = mean(PC2))
                                                 left_join(., means, by = "groups") %>%
                                                 {ggplot(., aes(PC1, PC2, color = groups, label = groups))+
                                                     geom_text_repel(data = dplyr::select(., mean.PC1, mean.PC2, groups) %>%
                                                                       unique(),
                                                                     aes(mean.PC1, mean.PC2), show.legend = F, size = 4)+
                                                     geom_point(show.legend = F)+
                                                     geom_segment(aes(xend = mean.PC1, yend = mean.PC2), show.legend = F)+
                                                     simple_theme+
                                                     scale_color_manual(values = groups.color)+
                                                     ggtitle("Scores")}} 
  
  ggsave(paste(outpath, paste0(prefix, '_PCA_scores.pdf'),sep='/'), width=10, height=10)

}


### 2. clustering plot
make_clustering_plot <- function(tb.wide, colors, outpath, prefix){
  consensus.cor.spearman<-cor(tb.wide, method="spearman", use="pairwise.complete.obs")
  d <- 1 - consensus.cor.spearman
  hcl <- hclust(as.dist(d), "average") 
  consensus.clust.order<-hcl$order
  consensus.colors <- colors[colnames(tb.wide)]
  
  pdf(file = paste(outpath, paste0(prefix, '_spearman_clustering.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  par(mar=c(10,4,4,2))
  print(plotColoredClusters(hcl, hcl$labels, cols=consensus.colors, cex=1,cex.lab=0.1, xlab="", sub=""))
  dev.off()
  
  my_palette <- colorRampPalette(c("white", "yellow", "red"))(n = 100)
  pdf(file = paste(outpath, paste0(prefix, '_spearman_heatmap.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  heatmap.2(consensus.cor.spearman, Rowv=as.dendrogram(hcl), Colv=rev(as.dendrogram(hcl)), 
            trace="none", na.col="gray" ,margins=c(8,8), cexCol=1, cexRow=1, col=my_palette, 
            RowSideColors = consensus.colors, ColSideColors = consensus.colors)
  dev.off()
}


# ------ umap, PCA, tsne plots -----------
make_umap_plot <- function(eset,outpath,prefix){
  require(umap)
  require(Rtsne)
  exprs <- t(exprs(eset))
  pdata <- pData(eset)
  embedding <- umap(exprs)
  umap.plot.matrix <-
    embedding$layout %>%
    as.tibble() %>%
    mutate(tissue = as.factor(pdata$tissue),
           sample = as.factor(pdata$sample))
  
  pdf(file = paste(outpath, paste0(prefix, '_umap.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(ggplot(umap.plot.matrix, aes(V1, V2, label = sample)) + 
    geom_point(aes(color = tissue),size=4, alpha=0.6)+
    geom_path(aes(color = tissue), size=0.5, linetype='dotted')+
    #  geom_line(aes(group = subj),size = 0.1)+
    geom_text_repel(aes(label=sample, color=tissue),size=3)+ #aes(color = Age_PN)
    #scale_color_viridis(option="D", direction = -1)+
    #stat_ellipse(aes(x=V1, y=V2,group = tissue),type = "t", linetype = "dotted", show.legend = F) +
    xlab('UMAP1')+
    ylab('UMAP2')+
    labs(title= paste0('umap: ',prefix))+
    theme_option_2)
  dev.off()

  exprs.merge <- merge(exprs,pdata,by.x="row.names",by.y="row.names",all.x=TRUE)
  rownames(exprs.merge) <- exprs.merge$Row.names
  exprs.merge.matrix <- exprs.merge[,-c(1,(ncol(exprs.merge)-1):ncol(exprs.merge))]
  tsne_out <- Rtsne(exprs.merge.matrix, perplexity=20)
  tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], tissue = exprs.merge$tissue, sample = exprs.merge$sample)
  
  pdf(file = paste(outpath, paste0(prefix, '_tsne.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(ggplot(tsne_plot,aes(x=x,y=y,label=sample)) + 
    geom_point(aes(x=x, y=y, color=tissue),size=4, alpha=0.6)+
    geom_path(aes(color = tissue), size=0.5, linetype='dotted')+
    geom_text_repel(aes(label=sample, color=tissue),size=3)+
    xlab('tSNE1')+
    ylab('tSNE2')+
    labs(title= paste0('tSNE: ',prefix))+
    theme_option_2)
  dev.off()
}


### 3. elevated bar plot
make_elevated_bar_plot <- function(elevated.summary.table, translate_categories = c("Tissue" = "Tissue"), outpath, prefix){
  pdf(file = paste(outpath, paste0(prefix, '_elevated_bar.pdf'),sep='/'), width=8, height=5, useDingbats = F)
  print(elevated.summary.table %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(tissue = names)} %>%
    gather(key = "Classification", value = "Number of genes", -tissue) %>%
    mutate(tissue = factor(tissue, levels = rev(unique(tissue[order(mapply(tissue, FUN = function(x) sum(`Number of genes`[tissue == x & Classification %in% c("Tissue enriched","Celltype enriched",
                                                                                                                                                               "Group enriched","Tissue enhanced","Celltype enhanced")])))]))),
           Classification = factor(Classification, levels = c("Not detected in any tissues","Not detected in any celltypes","Not detected in this tissue","Not detected in this celltype","Mixed in this tissue", "Mixed in this celltype","Expressed in all tissues","Expressed in all celltypes","Tissue enhanced", "Celltype enhanced","Group enriched","Tissue enriched", "Celltype enriched")),
           Classification = gsub(pattern = translate_categories, replacement = names(translate_categories), Classification)) %>%
    filter(Classification %in% c("Tissue enriched","Celltype enriched",
                                 "Group enriched","Tissue enhanced","Celltype enhanced")) %>%
    ggplot(aes(tissue, `Number of genes`, fill = Classification))+
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "",values = cat2.cols)+
    simple_theme+
    xlab("")+
    ylab("Number of genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = c(0.8, 0.8)))
  dev.off()
}

make_all_genes_bar_plot <- function(atlas.cat, translate_categories = c("Tissue" = "Tissue"), outpath, prefix){
  
  tissues <- 
    atlas.cat$`tissues over lim` %>%
    unique() %>%
    paste(collapse = ", ") %>%
    strsplit(split = ", ") %>%
    {.[[1]]} %>%
    unique() %>%
    {.[.[] != ""]}
  
  atlas.cat %>%
    {tibble(tissue = tissues, 
            `not expressed in any tissue` = filter(., `tissues over lim` == "") %>% nrow(),
            `not expressed in this tissue` = sapply(tissues, FUN = function(x) filter(., !grepl(x, `tissues over lim`)) %>% nrow()),
            `tissue enhanced` = sapply(tissues, 
                                       FUN = function(x) filter(., grepl(x, `enriched tissues`) & elevated.category == "tissue enhanced") %>% nrow()),
            `low tissue specificity` = sapply(tissues, 
                                              FUN = function(x) filter(., grepl(x, `tissues over lim`) & elevated.category == "low tissue specificity") %>% nrow()),
            `group enriched` = sapply(tissues, 
                                      FUN = function(x) filter(., grepl(x, `enriched tissues`) & elevated.category == "group enriched") %>% nrow()),
            `tissue enriched` = sapply(tissues, 
                                       FUN = function(x) filter(., grepl(x, `enriched tissues`) & elevated.category == "tissue enriched") %>% nrow()))} %>%
    mutate(`not expressed in this tissue` = `not expressed in this tissue` - `not expressed in any tissue`) %>%
    gather(key = "Classification", value = "Number of genes", -tissue) %>%
    mutate(tissue = factor(tissue, levels = rev(unique(tissue[order(mapply(tissue, FUN = function(x) sum(`Number of genes`[tissue == x & Classification %in% c("tissue enriched","celltype enriched",
                                                                                                                                                               "group enriched","tissue enhanced",
                                                                                                                                                               "celltype enhanced")])))]))),
           Classification = gsub(pattern = translate_categories, replacement = names(translate_categories), Classification),
           Classification = factor(Classification, levels = c("not expressed in any tissue","not expressed in any celltype",
                                                              "not expressed in this tissue","not expressed in this celltype",
                                                              "low tissue specificity","tissue enhanced", "celltype enhanced",
                                                              "group enriched","tissue enriched", "celltype enriched"))) %>%
    ggplot(aes(tissue, `Number of genes`, fill = Classification))+
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "",values = c(elevated.cat.cols, 
                                           "not expressed in any tissue" = "gray",
                                           "not expressed in any celltype" = "gray",
                                           "not expressed in this tissue" = "lightgray",
                                           "not expressed in this celltype" = "lightgray"))+
    simple_theme+
    xlab("")+
    ylab("Number of genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = c(0.8, 0.8))
  ggsave(paste(outpath, paste0(prefix, '_all_genes_bar.pdf'),sep='/'), width=8, height=5)
}

### 4. tissue distribution
make_tissue_distribution_plot <- function(tb.atlas, expr_column, und.lim = 1, do.tissues = "all", outpath, prefix) {
  pdf(file = paste(outpath, paste0(prefix, '_tissue_distribution.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(
  tb.atlas %>%
    mutate(tissue = factor(content_name, levels = with({group_by(., content_name) %>% 
        summarise(methods = paste(unique(method), collapse = " "))},content_name[order(methods)]))) %>%
        {ggplot(., aes(tissue, eval(parse(text = expr_column)), fill = method, color = method))+
            stat_summary(fun.y = "median", fun.args = c("na.rm" = T), geom = "line", aes(group = method), size = 2, alpha = 0.2)+
            stat_summary(fun.y = "min", geom = "line", aes(tissue, eval(parse(text = expr_column)), group = 1), inherit.aes = F, size = 1, alpha = 0.1)+
            stat_summary(fun.y = "max", geom = "line", aes(tissue, eval(parse(text = expr_column)), group = 1), inherit.aes = F, size = 1, alpha = 0.1)+
            stat_summary(fun.y = "min", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            stat_summary(fun.y = "max", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            
            
            geom_violin(draw_quantiles = 0.5, alpha = 0.2, position = "identity")+
            geom_hline(yintercept = und.lim, linetype = "dashed")+
            simple_theme+
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
            theme(axis.title = element_blank())+
            scale_fill_manual(values = dataset.colors)+
            scale_color_manual(values = dataset.colors)})
  dev.off()
}

## 5. specificity distribution
make_specificity_distribution_plot <- function(atlas.cat, type = "Tissue", outpath, prefix) {
  cat <- 
    atlas.cat %>% 
    group_by(elevated.category,express.category.2) %>% 
    summarize(num.genes=n()) %>% 
    ungroup() %>%
    mutate(elevated.category = factor(elevated.category, levels = rev(names(elevated.cat.cols))),
           express.category.2 = factor(express.category.2,levels = rev(names(expressed.cat.cols))))
  
  pdf(file = paste(outpath, paste0(prefix, '_Bar_specificity_distr_1.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(
    ggplot(cat, aes(x=express.category.2, y=num.genes, fill=elevated.category)) + 
    geom_bar(stat = "identity")+ 
    ggtitle(paste(type, "distribution category"))+ 
    scale_fill_manual(values=elevated.cat.cols, name=paste(type, "specificity category"))+
    ylab("Number of genes")+
    xlab("")+theme_light() +
    simple_theme+
    theme(axis.text=element_text(size=14),axis.text.x = element_text(angle=50, hjust=1),axis.title=element_text(size=14)))
  dev.off()
  
  pdf(file = paste(outpath, paste0(prefix, '_Bar_specificity_distr_2.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(
    ggplot(cat, aes(x=elevated.category, y=num.genes, fill=express.category.2)) + 
    geom_bar(stat = "identity") +
    ggtitle(paste(type, "specificity category")) + 
    scale_fill_manual(values=expressed.cat.cols, name=paste(type, "distribution category"))+
    ylab("Number of genes")+
    xlab("")+
    simple_theme+
    theme(axis.text=element_text(size=14),axis.text.x = element_text(angle=50, hjust=1),axis.title=element_text(size=14)))
  dev.off()
}

## 6. chord plot
chord_classification <- function(from, to, sizes, grid.col, groups = rep(1, length(from)), plot.order = c(unique(from), unique(to)), line_expansion = 10000, size_labels = F){
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, 
                 order = plot.order)
  
  if(size_labels) {
    for(i in 1:nrow(chord)) {
      value <- chord$value[i]
      if(is.null(value)) value <- chord$value1[1]
      x1 <- chord$x1[i] - value / 2
      x2 <- chord$x2[i] - value / 2
      to_ <- chord$cn[i]
      from_ <- chord$rn[i]
      circos.text(x = x1, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = from_, niceFacing = T)
      circos.text(x = x2, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = to_, niceFacing = T)
    }
  }
  
  
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)
    
    width <- strwidth(sector.name)*line_expansion
    
    circos.segments(x0 = mean(xlim), x1 = mean(xlim), 
                    y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment, 
                    sector.name)
    
    circos.segments(x0 = mean(xlim) - width/2, x1 = mean(xlim) + width/2, 
                    y0 = mean(ylim) - 0.2 + adjustment, y1 = mean(ylim) - 0.2 + adjustment, 
                    sector.name) 
    
    circos.text(mean(xlim), mean(ylim) + adjustment, sector.name, niceFacing = TRUE, facing = "bending")
  }, bg.border = NA)
  
  circos.clear()
}

make_classification_chord_plot <- function(atlas.cat, outpath, prefix) {
  pdf(file = paste(outpath, paste0(prefix, '_chord.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  atlas.cat %>%
    group_by(elevated.category,express.category.2) %>% 
    summarize(num.genes=n()) %$%
    chord_classification(from = elevated.category, 
                         to = express.category.2, 
                         sizes = num.genes, 
                         grid.col = c(elevated.cat.cols, expressed.cat.cols),#c(cat.cols, "not expressed" = "gray", "low tissue specificity" = "#4daf4a", "expressed in all" = "#377eb8", "expressed in some" = "#377eb8", "expressed in single" = "#377eb8", "expressed in many" = "#377eb8"),
                         groups = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                         plot.order = c(c("not detected", "low tissue specificity","tissue enhanced", "group enriched", "tissue enriched"), 
                                        c("not expressed","expressed in single", "expressed in some","expressed in many", "expressed in all")),
                         size_labels = T)
  dev.off()
}

## 7. swarm plot
make_swarm_expression_plot <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, outpath, prefix) {
  plot.data <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column))) %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    #Plot 1 % highest
    filter(expression >= quantile(expression, probs = 0.99)) %>%
    group_by(Grouping) %>%
    # Write gene name if highest 0.5 % per tissue and/or highest 1 % in total
    mutate(highest = expression >= quantile(expression, probs = 0.995) | rank(expression) >= (length(expression) - 8)) %>%
    ungroup() %>%
    mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98))) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`), by = "ensg_id") %>%
    left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
    left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))

  
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = elevated.category)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 2)+
      simple_theme+
      scale_color_manual(values = elevated.cat.cols)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_abundance_jitter_1.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 2)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_score_jitter_2.pdf'),sep='/'), width=15, height=10)
  
  plot.data <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column))) %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
    mutate(score = as.numeric(`tissue/group specific score`)) %>%
    # tissue enriched
    filter(elevated.category == "tissue enriched" & Grouping == `enriched tissues`) %>%
    group_by(Grouping) %>%
    # Write gene name if highest 2 % per tissue and/or highest 1 % in total
    mutate(highest = expression >= quantile(expression, probs = 0.99) | rank(expression) >= (length(expression) - 10),
           highest_score = score >= quantile(score, probs = 0.99) | rank(score) >= (length(score) - 10)) %>%
    ungroup() %>%
    mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98)),
           highest_score = ifelse(highest_score, T, score >= quantile(score, probs = 0.98))) %>%
    
    left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
    left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%

    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single)) 
  
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 2)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_tissue_enriched_jitter.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
  {ggplot(., aes(Grouping, score, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest_score,], size = 2)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_score_enriched_jitter.pdf'),sep='/'), width=15, height=10)
  
  plot.data <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column))) %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
    mutate(score = as.numeric(`tissue/group specific score`)) %>%
    # tissue enriched
    filter(elevated.category %in% c("tissue enriched", "tissue enhanced", "group enriched")) %>%
    filter(expression >= quantile(expression, probs = 0.99)) %>%
    group_by(Grouping) %>%
    #Plot 1 % highest
    #filter(expression >= quantile(expression, probs = 0.9)) %>%
    # Write gene name if highest 2 % per tissue and/or highest 1 % in total
    mutate(highest = expression >= quantile(expression, probs = 0.99) | rank(expression) >= (length(expression) - 10),
           highest_score = score >= quantile(score, probs = 0.99) | rank(score) >= (length(score) - 10)) %>%
    ungroup() %>%
    mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98)),
           highest_score = ifelse(highest_score, T, score >= quantile(score, probs = 0.98))) %>%
    
    left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
    left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
    
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single)) 
  
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 2)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_tissue_elevated_jitter.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
  {ggplot(., aes(Grouping, score, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest_score,], size = 2)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
  ggsave(paste(outpath, paste0(prefix, '_high_score_elevated_jitter.pdf'),sep='/'), width=15, height=10)
}

## 8. Bland-Altman plot

make_bland_altman_plot <- function(x, y, fill, fillname = "", title = "", Points = F, alpha = 0.1, 
                                   xl = "Mean of NX and X", yl = "NX - X"){
  
      
      
  
  if(Points) {
    p <- 
      tibble(x, y, fill) %>%
      {ggplot(., aes(x, y, color = fill)) + 
          geom_point(shape = 16, alpha = alpha) + 
          {ungroup(.) %>%
              summarise(SD = sd(y, na.rm = T),
                        Median = median(y, na.rm = T)) %$%
              geom_hline(yintercept = c(Median - 2 * SD, Median, Median + 2 * SD, 0), 
                         linetype = c("dotted", "dashed", "dotted", "solid"), 
                         color = c("black", "black", "black", "lightgray"))}} + 
      scale_color_discrete(name = fillname)
  } else {
    p <- 
      tibble(x, y, fill) %>%
      {ggplot(., aes(x, y, fill = fill, alpha = ..count..)) + 
          geom_hex(bins = 100, color = NA) +
          {ungroup(.) %>%
              summarise(SD = sd(y, na.rm = T),
                        Median = median(y, na.rm = T)) %$%
              geom_hline(yintercept = c(Median - 2 * SD, Median, Median + 2 * SD, 0), 
                         linetype = c("dotted", "dashed", "dotted", "solid"), 
                         color = c("black", "black", "black", "lightgray"))}} + 
      scale_fill_discrete(name = fillname)
  }
  
  p + 
    simple_theme + 
    xlab(xl) + 
    ylab(yl) + 
    ggtitle(title) 
    
    
}


make_chord_group_enriched <- function(elevated.table, grid.col, tissue_hierarcy, palet = colorRamps::blue2red, 
                                                reverse = F, outpath, prefix){
  circos.par(RESET = T)
  pdf(file = paste(outpath, paste0(prefix, '_group_enriched_hierarchy_chord.pdf'),sep='/'), width=15, height=10, useDingbats = F)
  
  content_gene_class <- 
    elevated.table %>%
    as_tibble(., rownames = "ensg_id") %>%
    gather(key = "content", "classification", -ensg_id) %>%
    # Only include group enriched
    filter(classification %in% c(3))
  
  mat <- 
    content_gene_class %>%
    {mapply(unique(.$content), 
            FUN = function(cell1) mapply(unique(.$content), 
                                         FUN = function(cell2) length(intersect(filter(., content == cell1)$ensg_id, 
                                                                                filter(., content == cell2)$ensg_id))))} 
  df <- 
    mat%>%
    as_tibble(., rownames = "from") %>%
    gather(key = "to", value = "number of genes", -from) %>%
    filter(to != from) %>%
    mutate(transfer = mapply(.$from, .$to, FUN = function(a,b) paste(sort(c(a,b)), collapse = " to "))) %>%
    # Take only first occurence of transfer
    {.[match(unique(.$transfer), .$transfer), ]}
  
  require(circlize)
  
  total.number.of.genes <- sapply(unique(c(df$from,df$to)), FUN = function(x) sum(df$`number of genes`[grep(x, df$transfer)]))
  
  classes <- unique(c(df$from, df$to))
  tissues <- tissue_hierarcy$content[match(classes, tissue_hierarcy$content)]
  track.levels <- names(tissue_hierarcy)
  tissue_hierarcy <- 
    tissue_hierarcy %>%
    filter(content %in% classes)
  
  chord <- 
    df %>% 
    mutate(from.sum = total.number.of.genes[match(from, names(total.number.of.genes))],
           to.sum = total.number.of.genes[match(to, names(total.number.of.genes))],
           frac.mean = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) mean(c(n/a, n/b))),
           frac.min = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) min(c(n/a, n/b))),
           frac.max = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) max(c(n/a, n/b)))) %>%
    {chordDiagram(.[,1:3],
                  grid.col = grid.col,
                  col = palet(100)[cut(.$frac.max, breaks = 100)],
                  directional = 0,link.largest.ontop = T,
                  annotationTrack="grid",
                  annotationTrackHeight = 0.02,
                  preAllocateTracks = eval(parse(text = paste0("list(", 
                                                               paste0(rep("list(track.height = 0.05, 
                                                                          track.margin = c(0.0035, 0.01))", 
                                                                          length(track.levels)), collapse = ", "), ")"))),
                  order = {if(length(track.levels) == 1) {
                    classes
                    } else {
                      classes[with(tissue_hierarcy[match(tissues, tissue_hierarcy[, 1][[1]]),],
                                   eval(parse(text = paste0("order(", paste(track.levels[-1], collapse = ", "), ")"))))]
                    }
                    })}
  
  
  
  sectors <-
    tissue_hierarcy %>%
    mutate(main.track.handle = content) %>%
    gather(key = "level", value = "handle", -main.track.handle) %>%
    mutate(color = grid.col[match(handle, names(grid.col))],
           track = ifelse(rep(reverse, nrow(.)), match(level, rev(track.levels)), match(level, track.levels))) %>%
    group_by(handle, color) %>%
    summarise(track = list(unique(track)),
              main.track.handles = list(unique(main.track.handle))) 
  
  for(tissue in tissue_hierarcy$content) {
    highlight.sector(tissue, 
                     track.index = 1:(length(track.levels)), 
                     col = NA, lwd = 2,
                     border = grid.col[match(tissue, names(grid.col))])
  }
  
  for(i in 1:nrow(sectors)) {
    highlight.sector(sectors$main.track.handles[i][[1]], 
                     track.index = sectors$track[i][[1]], 
                     col = sectors$color[i], 
                     text = sectors$handle[i], text.col = "white", facing = "bending", niceFacing = T, cex = 0.8)
  }
  
  dev.off()
  circos.clear()
  
  
}

chord_vertical_text <- function(from, to, sizes, grid.col, plot.group, plot.order){
  
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(plot.group)-1)) {
    if(plot.group[i] == plot.group[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(plot.group[length(plot.group)] == plot.group[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par, canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, order = plot.order)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    
    
    
    circos.text(mean(xlim), min(ylim), sector.name, adj = 0, niceFacing = TRUE, facing = "clockwise", cex = 0.9)
    
  }, bg.border = NA)
  
  circos.clear()
  
}


# ------ classification comparison plot

chord_classification_clockwise <- function(from, to, sizes, grid.col, plot.group, plot.order, line_expansion = 10000){
  
  require(circlize) 
  
  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)
  
  
  tb <- 
    tibble(from, to, sizes)
  
  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(plot.group)-1)) {
    if(plot.group[i] == plot.group[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }
  
  if(plot.group[length(plot.group)] == plot.group[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }
  
  circos.par(gap.after = gap.after.par)
  
  chord <-
    tb %>% 
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, order = plot.order)
  
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")
    
    adjustment <- (sector.index %% 5) * 0.2 - 0.4
    
    width <- strwidth(sector.name)*line_expansion
    
    circos.segments(x0 = mean(xlim), x1 = mean(xlim), 
                    y0 = min(ylim), y1 = mean(ylim)-0.12 + adjustment, 
                    sector.name)
    
    circos.segments(x0 = mean(xlim) - width/2, x1 = mean(xlim) + width/2, 
                    y0 = mean(ylim) - 0.12 + adjustment, y1 = mean(ylim) - 0.12 + adjustment, 
                    sector.name) 
    
    circos.text(mean(xlim), mean(ylim) + adjustment, sector.name, niceFacing = TRUE, facing = "bending", cex = 0.6)
    #circos.text(mean(xlim), ylim[1], sector.name, niceFacing = TRUE, facing = "clockwise",  adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  
  circos.clear()
  
}

make_class_comparison_chord <- function(cat1, cat2, line_expansion = 10000,
                                        outpath, prefix,
                                        excat_suffix = c(" celltypes", " tissues"),
                                        elcat_suffix = c(" celltypes", " tissues"),
                                        excat_prefix = c("", ""),
                                        elcat_prefix = c("", ""),
                                        excat_gsub_pattern = c(" tissue", ""),
                                        excat_gsub_replace = c(" celltype", ""),
                                        elcat_gsub_pattern = c(" tissue", ""),
                                        elcat_gsub_replace = c(" celltype", ""),
                                        excat_cats = c("expressed in all", "expressed in many", "expressed in some", "expressed in single", "not expressed"),
                                        elcat_cats = c("tissue enriched", "group enriched", "tissue enhanced", "low tissue specificity", "not detected")) {
  
  plot.order <- c(paste0(excat_prefix[1], gsub(excat_gsub_pattern[1], excat_gsub_replace[1], excat_cats), excat_suffix[1]),
                  paste0(elcat_prefix[1], gsub(elcat_gsub_pattern[1], elcat_gsub_replace[1], elcat_cats), elcat_suffix[1]),
                  rev(paste0(elcat_prefix[2], gsub(elcat_gsub_pattern[2], elcat_gsub_replace[2], elcat_cats), elcat_suffix[2])),
                  rev(paste0(excat_prefix[2], gsub(excat_gsub_pattern[2], excat_gsub_replace[2], excat_cats), excat_suffix[2])))
  chord_col <- setNames(c(expressed.cat.cols[match(excat_cats, names(expressed.cat.cols))], 
                          elevated.cat.cols[match(elcat_cats, names(elevated.cat.cols))],
                          rev(elevated.cat.cols[match(elcat_cats, names(elevated.cat.cols))]),
                          rev(expressed.cat.cols[match(excat_cats, names(expressed.cat.cols))])), 
                        plot.order)
  
  
  
  df <- 
    left_join(cat1, cat2, by = "ensg_id") %>% 
    rename(excat1 = express.category.2.x, excat2 = express.category.2.y, 
           elcat1 = elevated.category.x, elcat2 = elevated.category.y) %>%
    select(excat1, excat2, elcat1, elcat2) %>%
    mutate(excat1 = paste0(excat_prefix[1], gsub(excat_gsub_pattern[1], excat_gsub_replace[1], excat1), excat_suffix[1], "excat1"),
           excat2 = paste0(excat_prefix[2], gsub(excat_gsub_pattern[2], excat_gsub_replace[2], excat2), excat_suffix[2], "excat2"),
           elcat1 = paste0(elcat_prefix[1], gsub(elcat_gsub_pattern[1], elcat_gsub_replace[1], elcat1), elcat_suffix[1], "elcat1"),
           elcat2 = paste0(elcat_prefix[2], gsub(elcat_gsub_pattern[2], elcat_gsub_replace[2], elcat2), elcat_suffix[2], "elcat2")) %>% 
           {first = T
           for(col1 in names(.)){
             for(col2 in names(.)){
               tb.temp <- 
                 group_by(., eval(parse(text = col1)), eval(parse(text = col2))) %>%
                 dplyr::summarise(sizes = n()) %>%
                 rename("from" = 1,
                        "to" = 2)
               
               if(first) {
                 tb <- tb.temp
                 first = F
               } else{
                 tb <- 
                   rbind(tb, tb.temp)
               }
               
             }
           }
           tb
           } %>%
    filter(from != to) %>%
    mutate(transfer = mapply(from, to, 
                             FUN = function(a, b) {
                               paste(case_when(grepl("excat1", a) & grepl("excat2", b) ~ "excat",
                                               grepl("excat2", a) & grepl("excat1", b) ~ "excat rev",
                                               grepl("elcat1", a) & grepl("elcat2", b) ~ "elcat",
                                               grepl("elcat2", a) & grepl("elcat1", b) ~ "elcat rev",
                                               grepl("excat1", a) & grepl("elcat1", b) ~ "exel within 1",
                                               grepl("elcat1", a) & grepl("excat1", b) ~ "exel within 1 rev",
                                               grepl("excat2", a) & grepl("elcat2", b) ~ "exel within 2",
                                               grepl("elcat2", a) & grepl("excat2", b) ~ "exel within 2 rev",
                                               grepl("excat1", a) & grepl("elcat2", b) ~ "exel without",
                                               grepl("elcat2", a) & grepl("excat1", b) ~ "exel without rev",
                                               grepl("elcat1", a) & grepl("excat2", b) ~ "exel without",
                                               grepl("excat2", a) & grepl("elcat1", b) ~ "exel without rev"))})) %>%
    ungroup() %>%
    mutate(from = gsub("e(x|l)cat(1|2)$", "", from),
           to = gsub("e(x|l)cat(1|2)$", "", to)) %>%
    filter(!grepl("rev$", transfer))
  
  
  
  plot.color <- setNames(c(expressed.cat.cols,
                           rev(expressed.cat.cols),               
                           elevated.cat.cols[!grepl("celltype", names(elevated.cat.cols))],               
                           rev(elevated.cat.cols[!grepl("tissue", names(elevated.cat.cols))])),
                         plot.order)
  
  plot.group <- c(rep(1, length(expressed.cat.cols)),
                  rep(2, length(expressed.cat.cols)),
                  rep(3, length(elevated.cat.cols[!grepl("celltype", names(elevated.cat.cols))])),
                  rep(4, length(elevated.cat.cols[!grepl("tissue", names(elevated.cat.cols))])))
  
  for(filt in c("excat elcat", "exel within 1 exel within 2", "exel without")) {
    pdf(paste(outpath, paste0(prefix, " classification comparison ", filt, ".pdf"), sep = "/"), 
        width = 10, height = 10, useDingbats = F)
    df %>%
      filter(sapply(transfer, FUN = function(x) grepl(x, filt))) %$%
      chord_classification_clockwise(from, to, sizes, chord_col, plot.group, plot.order, line_expansion = 10000)
    dev.off()
  }
  
}





##
make_elevated_organ_total_chord <- function(cat1, cat2, line_expansion = 10000, 
                                            grid.col, elevated_cats = c(2,3,4), 
                                            direction = 1, cat1_name, cat2_name,
                                            outpath, prefix) {
  joined_cats <-
    left_join(cat1, cat2, by = "ensg_id") %>%
    filter(case_when(direction == 1 ~ category.x %in% elevated_cats, 
                     direction == 2 ~ category.y %in% elevated_cats)) %>%
    select(`enriched tissues.x`, `enriched tissues.y`, express.category.2.x, express.category.2.y, elevated.category.x, elevated.category.y)

  cat1_tissues <- unique(strsplit(paste(joined_cats$`enriched tissues.x`, collapse = ", "), split = ", ")[[1]])
  cat2_tissues <- unique(strsplit(paste(joined_cats$`enriched tissues.y`, collapse = ", "), split = ", ")[[1]])

  plot.data <- 
    crossing(cat1_tissues, cat2_tissues) %>%
    mutate(sizes = mapply(cat1_tissues, cat2_tissues, FUN = function(a, b) {
      filter(joined_cats, `enriched tissues.x` == a & `enriched tissues.y` == b) %>%
        nrow()
    })) %>%
    filter(sizes > 0) %>%
    mutate(cat2_tissues = ifelse(cat2_tissues == "", ifelse(direction == 1, 
                                                            paste("Not tissue elevated in", cat2_name), 
                                                            paste("Not tissue elevated in", cat1_name)), cat2_tissues),
           cat1_tissues = ifelse(cat1_tissues == "", ifelse(direction == 1, 
                                                            paste("Not tissue elevated in", cat2_name), 
                                                            paste("Not tissue elevated in", cat1_name)), cat1_tissues))
  
    
  cat1_factors <- unique(plot.data$cat1_tissues)
  cat2_factors <- unique(plot.data$cat2_tissues)
  
  
  number_cat1 <- 
    plot.data %>%
    group_by(cat1_tissues) %>%
    summarise(n = sum(sizes)) %>%
    filter(n > 0)
  
  number_cat2 <- 
    plot.data %>%
    group_by(cat2_tissues) %>%
    summarise(n = sum(sizes)) %>%
    filter(n > 0) 
  
  grid.col.cat1.name <- paste("Not tissue elevated in", cat1_name)
  grid.col.cat2.name <- paste("Not tissue elevated in", cat2_name)
  
  pdf(paste(outpath, paste0(prefix, " Elevated genes comparison local global.pdf"), sep = "/"), 
      width = 10, height = 10, useDingbats = F)
  plot.data %$%
    chord_vertical_text(from = cat1_tissues, to = cat2_tissues, sizes = sizes, c(grid.col, 
                                                                                 setNames(c("black", "black"), 
                                                                                          c(grid.col.cat1.name, grid.col.cat2.name))), 
                         plot.group = c(rep(1, length(cat1_factors)), rep(2, length(cat2_factors))), 
                         plot.order = c(number_cat1$cat1_tissues[order(number_cat1$n)], number_cat2$cat2_tissues[order(number_cat2$n)]))
  dev.off()
}
##
# Plot group enrich chord but as a heatmap

make_heatmap_group_enriched <- function(elevated.table, outpath, prefix) {
  
  group_enriched_number <- 
    elevated.table %>%
    as_tibble(., rownames = "ensg_id") %>%
    gather(key = "content", "classification", -ensg_id) %>%
    # Only include group enriched
    filter(classification %in% c(3)) %>%
    {mapply(unique(.$content), 
            FUN = function(cell1) mapply(unique(.$content), 
                                         FUN = function(cell2) length(intersect(filter(., content == cell1)$ensg_id, 
                                                                                filter(., content == cell2)$ensg_id))))} %>%
    as_tibble(., rownames = "from") %>%
    gather(key = "to", value = "number of genes", -from) %>%
    left_join(ungroup(.) %>%
                group_by(from) %>%
                summarise(total = sum(`number of genes`),
                          Min = min(`number of genes`),
                          Max = max(`number of genes`)),
              by = "from") %>%
    left_join(ungroup(.) %>%
                group_by(to) %>%
                summarise(total = sum(`number of genes`),
                          Min = min(`number of genes`),
                          Max = max(`number of genes`)),
              by = "to") %>%
    mutate(Min = mapply(Min.x, Min.y, FUN = function(a, b) min(c(a, b))),
           Max = mapply(Max.x, Max.y, FUN = function(a, b) max(c(a, b))),
           n = (`number of genes` - Min)/ (Max - Min),
           from.sum = filter(., to == from)$`number of genes`[match(from, filter(., to == from)$from)],
           to.sum = filter(., to == from)$`number of genes`[match(to, filter(., to == from)$from)],
           frac.max = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) max(c(n/a, n/b)))) 
  
  GEN_dendrogram <- 
    group_enriched_number %>%
    select(from, to, frac.max) %>%
    spread(key = to, value = frac.max) %>%
    column_to_rownames("from") %>%
    as.matrix() %>%
    cor(method="spearman", use="pairwise.complete.obs")  %>%
    {1 - .} %>%
    as.dist() %>%
    hclust("average")  %>%
    as.dendrogram() 
  
  
  group_enriched_number %>%
    mutate(from = factor(from, levels = unique(from)[order.dendrogram(GEN_dendrogram)]),
           to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
    ggplot(aes(from, to, fill = frac.max)) +
    geom_tile(alpha = 0.8) + 
    simple_theme + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), axis.title = element_blank()) + 
    scale_fill_gradientn(colors = c("yellow", "orangered", "#800026")) 
  
  ggsave(paste(outpath, paste0(prefix, '_group_enriched_heatmap.pdf'),sep='/'), width=15, height=10)
}


# Plot group enrich chord but as a heatmap

make_heatmap_group_enriched_circle <- function(elevated.table, outpath, prefix) {
  
  group_enriched_number <- 
    elevated.table %>%
    as_tibble(., rownames = "ensg_id") %>%
    gather(key = "content", "classification", -ensg_id) %>%
    # Only include group enriched
    filter(classification %in% c(3)) %>%
    {mapply(unique(.$content), 
            FUN = function(cell1) mapply(unique(.$content), 
                                         FUN = function(cell2) length(intersect(filter(., content == cell1)$ensg_id, 
                                                                                filter(., content == cell2)$ensg_id))))} %>%
    as_tibble(., rownames = "from") %>%
    gather(key = "to", value = "number of genes", -from) %>%
    left_join(ungroup(.) %>%
                group_by(from) %>%
                summarise(total = sum(`number of genes`),
                          Min = min(`number of genes`),
                          Max = max(`number of genes`)),
              by = "from") %>%
    left_join(ungroup(.) %>%
                group_by(to) %>%
                summarise(total = sum(`number of genes`),
                          Min = min(`number of genes`),
                          Max = max(`number of genes`)),
              by = "to") %>%
    mutate(Min = mapply(Min.x, Min.y, FUN = function(a, b) min(c(a, b))),
           Max = mapply(Max.x, Max.y, FUN = function(a, b) max(c(a, b))),
           n = (`number of genes` - Min)/ (Max - Min),
           from.sum = filter(., to == from)$`number of genes`[match(from, filter(., to == from)$from)],
           to.sum = filter(., to == from)$`number of genes`[match(to, filter(., to == from)$from)],
           frac.max = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) max(c(n/a, n/b)))) 
  
  GEN_dendrogram <- 
    group_enriched_number %>%
    select(from, to, frac.max) %>%
    spread(key = to, value = frac.max) %>%
    column_to_rownames("from") %>%
    as.matrix() %>%
    cor(method="spearman", use="pairwise.complete.obs")  %>%
    {1 - .} %>%
    as.dist() %>%
    hclust("average")  %>%
    as.dendrogram()
  
  dend_segments <- 
    ggdendro::dendro_data(GEN_dendrogram)$segments
    
  group_enriched_number %>%
    mutate(from = factor(from, levels = unique(from)[order.dendrogram(GEN_dendrogram)]),
           to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
    {nfactors  <- length(levels(.$from));
    n_sectors <- nfactors + 3
    text_angle <- function(x) - ((x + 3) * 360/n_sectors - 360/(n_sectors*2))
    ggplot(.) +
        geom_tile(aes(from, to, fill = frac.max), alpha = 0.8) + 
        geom_segment(data = dend_segments, aes(x = x, y = -10*yend, xend = xend, yend = -10*y))+
        annotate(geom = "text", x = 1:nfactors, y = n_sectors, label = levels(.$from), 
                 angle = text_angle(1:nfactors),
                 size = 3)+
        annotate(geom = "text", x = -2.5, y =1:nfactors, label = levels(.$from), hjust = 0,
                 size = 3)+
        annotate("point", x = -2.5, y = n_sectors, color = NA)+
        #coord_fixed() +
        coord_polar()+
        #simple_theme + 
        theme_minimal()+
        theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
        scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
  
  ggsave(paste(outpath, paste0(prefix, '_group_enriched_heatmap_circle.pdf'),sep='/'), width=15, height=10)
}



make_heatmap_group_and_enhanced_expression_levels_circle <- 
  function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix) {
    
    group_enriched_genes <- 
      elevated.table %>%
      as_tibble(., rownames = "ensg_id") %>%
      gather(key = "content", "classification", -ensg_id) %>%
      # Only include group enriched
      filter(classification %in% c(3, 4))
    
    group_enriched_expression <- 
      group_enriched_genes %>%
      {mapply(unique(.$content), 
              FUN = function(cell1) mapply(unique(.$content), 
                                           FUN = function(cell2) {
                                             genes <- 
                                               filter(., content == cell1)$ensg_id
                                             all.atlas.max.tb %>%
                                               filter(eval(parse(text = tissue_column)) == cell2 & 
                                                        ensg_id %in% genes) %$%
                                               median(eval(parse(text = maxEx_column)), na.rm = T)
                                           }))} %>%
      as_tibble(., rownames = "from") %>%
      gather(key = "to", value = "median expression", -from) 
    
    GEN_dendrogram <- 
      group_enriched_expression %>%
      select(from, to, `median expression`) %>%
      spread(key = to, value = `median expression`) %>%
      column_to_rownames("from") %>%
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram()
    
    dend_segments <- 
      ggdendro::dendro_data(GEN_dendrogram)$segments
    
    group_enriched_expression %>%
      mutate(from = factor(from, levels = unique(from)[order.dendrogram(GEN_dendrogram)]),
             to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
             {nfactors  <- length(levels(.$from));
             n_sectors <- nfactors + 3
             text_angle <- function(x) - ((x + 3) * 360/n_sectors - 360/(n_sectors*2))
             ggplot(.) +
               geom_tile(aes(from, to, fill = `median expression`), alpha = 0.8) + 
               geom_segment(data = dend_segments, aes(x = x, y = -10*yend, xend = xend, yend = -10*y))+
               annotate(geom = "text", x = 1:nfactors, y = n_sectors, label = levels(.$from), 
                        angle = text_angle(1:nfactors),
                        size = 3)+
               annotate(geom = "text", x = -2.5, y =1:nfactors, label = levels(.$from), hjust = 0,
                        size = 3)+
               annotate("point", x = -2.5, y = n_sectors, color = NA)+
               #coord_fixed() +
               coord_polar()+
               #simple_theme + 
               theme_minimal()+
               theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
               scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
    
    ggsave(paste(outpath, paste0(prefix, '_group_and_enhanced_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
  }

make_heatmap_group_enriched_expression_levels_circle <- 
  function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix, y_dendrogram = F) {
    
    group_enriched_genes <- 
      elevated.table %>%
      as_tibble(., rownames = "ensg_id") %>%
      gather(key = "content", "classification", -ensg_id) %>%
      # Only include group enriched
      filter(classification %in% c(3))
    
    group_enriched_expression <- 
      group_enriched_genes %>%
      right_join(filter(all.atlas.max.tb, ensg_id %in% group_enriched_genes$ensg_id), by = c("ensg_id", "content" = tissue_column)) %$%
      tibble(from = ensg_id,
             to = content,
             expression = log10(eval(parse(text = maxEx_column)) + 1)) 
      
    
    GEN_dendrogram <- 
      group_enriched_expression %>%
      select(from, to, expression) %>%
      spread(key = to, value = expression) %>%
      column_to_rownames("from") %>%
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram()
    
    
    
    dend_segments <- 
      ggdendro::dendro_data(GEN_dendrogram)$segments
    
    
    
    if(y_dendrogram) {
      Genes_dendrogram <- 
        group_enriched_expression %>%
        select(from, to, expression) %>%
        spread(key = from, value = expression) %>%
        column_to_rownames("to") %>%
        as.matrix() %>%
        cor(method="spearman", use="pairwise.complete.obs")  %>%
        {1 - .} %>%
        as.dist() %>%
        hclust("average")  %>%
        as.dendrogram()
      genes_dend_segments <- 
        ggdendro::dendro_data(Genes_dendrogram)$segments
    }
    
    x_margin = 2.5
    
    g <-
      group_enriched_expression %>%
      #filter(from %in% .$from[1:100]) %>%
      mutate(from = factor(from, levels = unique(from)[order.dendrogram(Genes_dendrogram)]),
             to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
             {nfactors  <- length(levels(.$to));
             n_sectors <- round((1 + 1/6)*nfactors)
             y_sectors <- length(levels(.$from))
             

             text_angle <- function(x) - (round((x + (1/6)*nfactors)) * 360/n_sectors - 360/(n_sectors*2))
             ggplot(.) +
               geom_tile(aes(to, from, fill = expression), alpha = 0.8) + 
               geom_segment(data = dend_segments, aes(x = x, y = -0.5*y_sectors*yend, 
                                                      xend = xend, yend = -0.5*y_sectors*y))+
               
               annotate(geom = "text", x = 1:nfactors, y = 1.15*y_sectors, label = levels(.$to),
                        angle = text_angle(1:nfactors),
                        size = 3)+
               annotate("point", x = -x_margin, y = n_sectors, color = NA)+
               #coord_fixed() +
               coord_polar()+
               #simple_theme + 
               theme_minimal()+
               theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
               scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
    
    if(y_dendrogram){
      genes_dend_segments_scaled <-
      {max_y = max(c(genes_dend_segments$y,
                     genes_dend_segments$yend))
      min_y = min(c(genes_dend_segments$y,
                    genes_dend_segments$yend))
      genes_dend_segments %>%
        mutate(y = ((y - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin),
               yend = ((yend - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin))
      }
      g <- 
        g + 
        geom_segment(data = genes_dend_segments_scaled, aes(x = y, y = x,
                                                            xend = yend, yend = xend))
      
    }
    
    
    ggsave(plot = g, paste(outpath, paste0(prefix, '_group_enriched_gene_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
  }

make_heatmap_all_elevated_expression_levels_circle <- 
  function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix, y_dendrogram = F) {
    
    group_enriched_genes <- 
      elevated.table %>%
      as_tibble(., rownames = "ensg_id") %>%
      gather(key = "content", "classification", -ensg_id) %>%
      # Only include group enriched
      filter(classification %in% c(2, 3, 4))
    
    group_enriched_expression <- 
      group_enriched_genes %>%
      right_join(filter(all.atlas.max.tb, ensg_id %in% group_enriched_genes$ensg_id), by = c("ensg_id", "content" = tissue_column)) %$%
      tibble(from = ensg_id,
             to = content,
             expression = log10(eval(parse(text = maxEx_column)) + 1)) 
    
    
    GEN_dendrogram <- 
      group_enriched_expression %>%
      select(from, to, expression) %>%
      spread(key = to, value = expression) %>%
      column_to_rownames("from") %>%
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram()
    
    
    
    dend_segments <- 
      ggdendro::dendro_data(GEN_dendrogram)$segments
    
    
    
    Genes_dendrogram <- 
      group_enriched_expression %>%
      select(from, to, expression) %>%
      spread(key = from, value = expression) %>%
      column_to_rownames("to") %>%
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram()
    genes_dend_segments <- 
      ggdendro::dendro_data(Genes_dendrogram)$segments
    
    
    x_margin = 2.5
    
    g <-
      group_enriched_expression %>%
      #filter(from %in% .$from[1:100]) %>%
      mutate(from = factor(from, levels = unique(from)[order.dendrogram(Genes_dendrogram)]),
             to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
             {nfactors  <- length(levels(.$to));
             n_sectors <- round((1 + 1/6)*nfactors)
             y_sectors <- length(levels(.$from))
             
             
             text_angle <- function(x) - (round((x + (1/6)*nfactors)) * 360/n_sectors - 360/(n_sectors*2))
             ggplot(.) +
               geom_tile(aes(to, from, fill = expression), alpha = 0.8) + 
               geom_segment(data = dend_segments, aes(x = x, y = -0.5*y_sectors*yend, 
                                                      xend = xend, yend = -0.5*y_sectors*y))+
               
               annotate(geom = "text", x = 1:nfactors, y = 1.15*y_sectors, label = levels(.$to),
                        angle = text_angle(1:nfactors),
                        size = 3)+
               annotate("point", x = -x_margin, y = n_sectors, color = NA)+
               #coord_fixed() +
               coord_polar()+
               #simple_theme + 
               theme_minimal()+
               theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
               scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
    
    if(y_dendrogram){
      genes_dend_segments_scaled <-
      {max_y = max(c(genes_dend_segments$y,
                     genes_dend_segments$yend))
      min_y = min(c(genes_dend_segments$y,
                    genes_dend_segments$yend))
      genes_dend_segments %>%
        mutate(y = ((y - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin),
               yend = ((yend - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin))
      }
      g <- 
        g + 
        geom_segment(data = genes_dend_segments_scaled, aes(x = y, y = x,
                                                            xend = yend, yend = xend))
      
    }
    
    
    ggsave(plot = g, paste(outpath, paste0(prefix, '_all_elevated_gene_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
  }

make_spearman_dendrogram <- function(var1, var2, value){
  data.frame(var1, var2, value) %>%
    spread(key = var1, value = value) %>%
    column_to_rownames("var2") %>% 
    as.matrix() %>%
    cor(method="spearman", use="pairwise.complete.obs")  %>%
    {1 - .} %>%
    as.dist() %>%
    hclust("average")  %>%
    as.dendrogram() 
}

get_dendrogram_segments <- function(dendrogram) {
  dendrogram %>%
    ggdendro::dendro_data() %$%
    segments
}

range_scale <- function(x, xmax, xmin, span, dodge = 0) { 
  ((x - xmin)/(xmax - xmin))*
    (span[2] - span[1]) + span[1] + 
    dodge
}

ggdendroheat <- function(x, y, value, show.legend = T, xdendrogram = T, ydendrogram = T,
                         x_margin = 0.2,
                         y_margin = 0.2){
  
  x_dendrogram <- 
    make_spearman_dendrogram(x, y, value)
  x_dendrogram_segments <- get_dendrogram_segments(x_dendrogram)
  
  y_dendrogram <- 
    make_spearman_dendrogram(y, x, value)
  y_dendrogram_segments <- get_dendrogram_segments(y_dendrogram)
  
  g <- 
    tibble(x, y, value) %>%
    mutate(x = factor(x, levels = labels(x_dendrogram)),
           y = factor(y, levels = labels(y_dendrogram))) %>%
    ggplot() +
    geom_tile(aes(x, y, fill = value), show.legend = show.legend) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8), 
          panel.grid = element_blank())
  
  xp <- 0
  yp <- 0
  
  
  x_factors <- length(labels(x_dendrogram))
  y_factors <- length(labels(y_dendrogram))
  
  if(xdendrogram) {
    g <- g + 
      geom_segment(data = x_dendrogram_segments, 
                   aes(x = x, 
                       xend = xend,
                       y = range_scale(yend, max(y, yend), min(y, yend), 
                                       span = c(y_factors, y_factors/(1 - y_margin)), 
                                       dodge = 0.5), 
                       yend = range_scale(y, max(y, yend), min(y, yend), 
                                          span = c(y_factors, y_factors/(1 - y_margin)), 
                                          dodge = 0.5)))
    xp <- c(x_factors + x_margin)
  }
  
  if(ydendrogram) {
    g <- g + 
      geom_segment(data = y_dendrogram_segments, 
                   aes(x = range_scale(yend, max(y, yend), min(y, yend), 
                                       span = c(x_factors, x_factors/(1 - x_margin)), 
                                       dodge = 0.5), 
                       xend = range_scale(y, max(y, yend), min(y, yend), 
                                          span = c(x_factors, x_factors/(1 - x_margin)), 
                                          dodge = 0.5),
                       y = x, 
                       yend = xend)) 
    yp <- c(y_factors + y_margin)
  }
  
  g + annotate("point", x = xp, y = yp, color = NA)
  
  
}

make_heatmap_expression_levels <- function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, enrichment = c(3), outpath, prefix, y_dendrogram = F) {
  genes <- 
    elevated.table %>%
    as_tibble(., rownames = "ensg_id") %>%
    gather(key = "content", "classification", -ensg_id) %>%
    filter(classification %in% enrichment)
  
  genes_expression <- 
    genes %>%
    right_join(filter(all.atlas.max.tb, ensg_id %in% genes$ensg_id), by = c("ensg_id", "content" = tissue_column)) %$%
    tibble(from = ensg_id,
           to = content,
           expression = log10(eval(parse(text = maxEx_column)) + 1)) 
  
  g <- 
    genes_expression %$%
    ggdendroheat(from, to, expression)+
    scale_fill_gradientn(colors = c("yellow", "orangered", "#800026")) + 
    theme(axis.text.x = element_blank(), axis.title.y = element_blank()) + 
    xlab("genes")
  
  ggsave(plot = g, paste(outpath, paste0(prefix, '__gene_expression_heatmap.pdf'),sep='/'), width=15, height=10)
}

make_heatmap_median_expression_levels <- 
  function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, enrichment = c(3), outpath, prefix, y_dendrogram = F) {
    genes <- 
      elevated.table %>%
      as_tibble(., rownames = "ensg_id") %>%
      gather(key = "content", "classification", -ensg_id) %>%
      filter(classification %in% enrichment)
    
    genes_expression <- 
      genes %>%
      {mapply(unique(.$content), 
              FUN = function(cell1) mapply(unique(.$content), 
                                           FUN = function(cell2) {
                                             genes <- 
                                               filter(., content == cell1)$ensg_id
                                             all.atlas.max.tb %>%
                                               filter(eval(parse(text = tissue_column)) == cell2 & 
                                                        ensg_id %in% genes) %$%
                                               median(eval(parse(text = maxEx_column)), na.rm = T)
                                           }))} %>%
      as_tibble(., rownames = "from") %>%
      gather(key = "to", value = "median expression", -from)
    
    g <- 
      genes_expression %$%
      ggdendroheat(from, to, `median expression`)+
      scale_fill_gradientn(colors = c("yellow", "orangered", "#800026")) + 
      theme(axis.title = element_blank())
    
    ggsave(plot = g, paste(outpath, paste0(prefix, '__median_expression_heatmap.pdf'),sep='/'), width=15, height=10)
    
    }



make_number_detected_genes_barplot <- function(all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix) {
  all.atlas.max.tb %>%
    filter(eval(parse(text = maxEx_column)) >= 1) %>%
    mutate(content = eval(parse(text = tissue_column))) %>%
    group_by(content) %>%
    summarise(`number expressed` = length(ensg_id)) %>%  
    mutate(content = factor(content, levels = content[order(`number expressed`, decreasing = T)])) %>%
    ggplot(aes(content, `number expressed`)) +
    geom_bar(stat = "identity", fill = "gray80", color = "black") +
    
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8), 
          axis.title.x = element_blank())+
    ylab("Number of expressed genes")
  
  ggsave(paste(outpath, paste0(prefix, '__number_detected_genes_barplot.pdf'),sep='/'), width=6, height=4)
}

# Plot method spearman cluster

make_spearman_method_dendrogram <- function(all.atlas.tb, Ex_column, content_column, named_color_replacement, outpath, prefix) {
  
  
  dendr <- 
    all.atlas.tb %>%
    mutate(content_method = paste(eval(parse(text = content_column)), method),
           ex = eval(parse(text = Ex_column)),
           ex = log(ex + 1)) %$%
    make_spearman_dendrogram(content_method, ensg_id, ex)
  
  dendr %>%
    get_dendrogram_segments() %>%
    {ggplot(., aes(x, y, xend = xend, yend = yend)) +
        geom_segment() +
        geom_text(data = filter(., yend == 0), aes(y = yend - 0.002, 
                                                   label = labels(dendr), 
                                                   color = trimws(str_extract(labels(dendr), paste(names(named_color_replacement), collapse = "|")))), 
                  angle = 90, hjust = 1, vjust = 0.3, size = 3, show.legend = F) + 
        scale_y_continuous(expand = c(0.05,0.2)) +
        scale_color_manual(values = named_color_replacement) +
        theme_minimal() + 
        theme(axis.text = element_blank(), axis.title = element_blank(), 
              panel.grid = element_blank())}
    
  
   ggsave(paste(outpath, paste0(prefix, '_spearman_method_cluster.pdf'),sep='/'), width=16, height=8)
}
