library('pcaMethods')
library('ggrepel')
library('gridExtra')
library('ClassDiscovery')
library('gplots')

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
  pdf(file = paste(outpath, paste0(prefix, '_PCA_loadings.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(loadings %>% 
  {ggplot(., aes(PC1, PC2, label = labels))+
      geom_text(show.legend = F)+
      simple_theme+
      ggtitle("Loadings")})
  dev.off()
  
  pdf(file = paste(outpath, paste0(prefix, '_PCA_scores.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(
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
                                                     #scale_color_manual(values = groups.color)+
                                                     ggtitle("Scores")}})
  dev.off()
    # scores %>%
    # {names <- rownames(.); as.tibble(.) %>% mutate(sample = names,
    #                                                groups = groups[match(names, names(groups))])} %>%
    #                                                {means <- group_by(., groups) %>%
    #                                                  dplyr::summarise(mean.PC1 = mean(PC1),
    #                                                                   mean.PC2 = mean(PC2))
    #                                                left_join(., means, by = "groups") %>%
    #                                                {ggplot(., aes(PC1, PC2, color = groups, label = groups))+
    # 
    #                                                    geom_point(show.legend = F)+
    #                                                    geom_segment(aes(xend = mean.PC1, yend = mean.PC2), show.legend = F)+
    #                                                    simple_theme+
    #                                                    #scale_color_manual(values = groups.color)+
    #                                                    ggtitle("Scores")}}
  #grid.arrange(g1, g2, g3, ncol = 3)

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


### 3. elevated bar plot
make_elevated_bar_plot <- function(elevated.summary.table, outpath, prefix){
  pdf(file = paste(outpath, paste0(prefix, '_elevated_bar.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  print(elevated.summary.table %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(tissue = names)} %>%
    gather(key = "Classification", value = "Number of genes", -tissue) %>%
    mutate(tissue = factor(tissue, levels = rev(unique(tissue[order(mapply(tissue, FUN = function(x) sum(`Number of genes`[tissue == x & Classification %in% c("Tissue enriched","Celltype enriched",
                                                                                                                                                               "Group enriched","Tissue enhanced","Celltype enhanced")])))]))),
           Classification = factor(Classification, levels = c("Not detected in any tissues","Not detected in any celltypes","Not detected in this tissue","Not detected in this celltype","Mixed in this tissue", "Mixed in this celltype","Expressed in all tissues","Expressed in all celltypes","Tissue enhanced", "Celltype enhanced","Group enriched","Tissue enriched", "Celltype enriched"))) %>%
    filter(Classification %in% c("Tissue enriched","Celltype enriched",
                                 "Group enriched","Tissue enhanced","Celltype enhanced")) %>%
    ggplot(aes(tissue, `Number of genes`, fill = Classification))+
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "",values = cat2.cols)+
    simple_theme+
    xlab("")+
    ylab("Number of genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = c(0.6, 0.8)))
  dev.off()
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
chord_classification <- function(from, to, sizes, grid.col, groups = rep(1, length(from)), plot.order = c(unique(from), unique(to)), line_expansion = 10000){
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
                 directional = 1,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05, 
                 preAllocateTracks = 1, order = plot.order)
  
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
                         grid.col = c(cat.cols, "not expressed" = "gray", "low tissue specificity" = "#4daf4a", "expressed in all" = "#377eb8", "expressed in some" = "#377eb8", "expressed in single" = "#377eb8", "expressed in many" = "#377eb8"),
                         groups = c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
                         plot.order = c("not detected", "low tissue specificity","tissue enhanced", "group enriched", "tissue enriched", 
                                        "not expressed","expressed in single", "expressed in some","expressed in many", "expressed in all"))
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
    # Write gene name if highest 2 % per tissue and/or highest 1 % in total
    mutate(highest = expression >= quantile(expression, probs = 0.99)) %>%
    ungroup() %>%
    mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98))) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category), by = "ensg_id") %>%
    left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
    left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))

  pdf(file = paste(outpath, paste0(prefix, '_high_abundance_jitter_1.pdf'),sep='/'), width=15, height=10, useDingbats = F)
  print(
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = elevated.category)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 1.5)+
      simple_theme+
      scale_color_manual(values = elevated.cat.cols)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))})
  dev.off()
  
  pdf(file = paste(outpath, paste0(prefix, '_high_abundance_jitter_2.pdf'),sep='/'), width=15, height=10, useDingbats = F)
  print(
  plot.data %>%
  {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
      geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
      geom_text_repel(data = .[.$highest,], size = 1.5)+
      simple_theme+
      scale_color_manual(values = protein.class.palette)+
      scale_y_log10()+
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))})
  dev.off()
}
