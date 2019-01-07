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
  loadings %>% 
    {ggplot(., aes(PC1, PC2, label = labels))+
      geom_text(show.legend = F)+
      simple_theme+
      ggtitle("Loadings")} %>%
    print()
  dev.off()
  
  pdf(file = paste(outpath, paste0(prefix, '_PCA_scores.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  
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
                                                     ggtitle("Scores")}} %>%
    print()
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
    chordDiagram(grid.col = grid.col,
                 col = palet(100)[cut(.$frac.max, breaks = 100)],
                 directional = 0,link.largest.ontop = T,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.02,
                 preAllocateTracks = eval(parse(text = paste0("list(", 
                                                              paste0(rep("list(track.height = 0.05, 
                                                                         track.margin = c(0.0035, 0.01))", 
                                                                         length(track.levels)), collapse = ", "), ")"))),
                 order = classes[with(tissue_hierarcy[match(tissues, tissue_hierarcy[, 1][[1]]),],
                                      eval(parse(text = paste0("order(", paste(track.levels[-1], collapse = ", "), ")"))))])
  
  
  
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
    chordDiagram(grid.col = chord_col,
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
  
  for(filt in c("excat elcat", "exel within 1 exel within 2", "exel without")) {
    pdf(paste(outpath, paste0(prefix, " classification comparison ", filt, ".pdf"), sep = "/"), 
        width = 10, height = 10, useDingbats = F)
    df %>%
      filter(sapply(transfer, FUN = function(x) grepl(x, filt))) %$%
      chord_classification_clockwise(from, to, sizes, chord_col, plot.group, plot.order, line_expansion = 10000)
    dev.off()
  }
  
}



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


# For enriched genes NX between tissues 

# atlas_max <- all.atlas.max
# cats <- all.atlas.category
# elevated.table <- all.atlas.elevated.table
# 
# left_join(atlas_max, cats, by = "ensg_id") %>%
#   select(ensg_id, elevated.category, limma_gene_dstmm.zero.impute.expression_maxEx, consensus_content_name, `enriched tissues`) %>%
#   filter(elevated.category %in% c("group enriched", "tissue enriched")) %>%
#   mutate(NX = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
#   filter(!is.na(NX)) %>%
#   mutate(content = ifelse(is.na(content), `enriched tissues`, content)) %>%
#   left_join({group_by(., content) %>%
#       summarise(mean_NX = mean(NX, na.rm = T),
#                 sd_NX = sd(NX))}, by = c("content")) %>%
#   mutate(scaled_NX = (NX - mean_NX) / sd_NX) %>%
#   # left_join({ungroup(.) %>%
#   #     group_by(`enriched tissues`) %>%
#   #     summarise(mean_NX_tiss = mean(NX),
#   #               sd_NX_tiss = sd(NX))}, by = c("enriched tissues")) %>%
#   # mutate(scaled_NX = (scaled_NX - mean_NX_tiss) / sd_NX_tiss) %>%
#   ungroup() %>%
#   group_by(., consensus_content_name, content) %>%
#   summarise(number = length(content),
#             mean_scaled_NX = mean(scaled_NX, na.rm = T)) %>%
#   ggplot(aes(consensus_content_name, content, fill = number)) +
#   geom_tile() + 
#   simple_theme + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3)) + 
#   scale_fill_viridis_c()