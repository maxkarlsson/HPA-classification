library('pcaMethods')
library('ggrepel')
library('gridExtra')
library('ClassDiscovery')
library('gplots')


make_general_environment_plots <- function(outpath) {
  
  plot.data <- 
    all.atlas.max %>%
    select(ensg_id, 
           consensus_content_name, 
           norm_exp = max_norm_exp) %>%
    mutate(atlas = "tissue") %>%
    # rbind(all.atlas.max %>%
    #         select(ensg_id, 
    #                consensus_content_name, 
    #                norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
    #         mutate(atlas = "tissue"),
    #       brain.atlas.max %>%
    #         select(ensg_id, 
    #                consensus_content_name = subgroup, 
    #                norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
    #         mutate(atlas = "brain")) %>%
    rbind(blood.atlas.max %>%
            select(ensg_id, 
                   consensus_content_name, 
                   norm_exp = max_norm_exp) %>%
            mutate(atlas = "blood cells")) %>% 
    ungroup() %>%
    rbind(select(cell.lines.atlas, 1:3) %>%
            rename(consensus_content_name = celline_name) %>%
            mutate(atlas = "celline")) %>%
    filter(norm_exp >= 1) %>%
    group_by(consensus_content_name, atlas) %>%
    summarise(n_detected_genes = length(ensg_id)) %>% 
    ungroup() %>%
    mutate(consensus_content_name = factor(consensus_content_name, levels = unique(consensus_content_name[order(atlas, n_detected_genes)]))) 
  
  plot.data %>%
    ggplot(aes(consensus_content_name, n_detected_genes, fill = atlas)) +
    geom_bar(stat = "identity") +
    simple_theme + 
    theme(axis.text.x = element_blank())
  ggsave(paste(outpath,'detected_genes_bar.pdf',sep='/'), width=8, height=8)
  
  plot.data %>%
    ggplot(aes(atlas, n_detected_genes, fill = atlas, color = atlas)) +
    geom_boxplot(alpha = 0.2) + 
    ggbeeswarm::geom_beeswarm() +
    simple_theme
  ggsave(paste(outpath,'detected_genes_boxplot.pdf',sep='/'), width=8, height=8)
  
  plot.data <- 
    all.atlas.max %>%
    select(ensg_id, 
           consensus_content_name, 
           norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
    filter(!consensus_content_name %in% c("brain", "blood", "intestine", "lymphoid system")) %>%
    mutate(atlas = "tissues") %>%
    ungroup() %>%
    rbind(all.atlas %>%
            select(ensg_id, 
                   content_name,
                   consensus_content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression) %>%
            filter(consensus_content_name %in% c("brain", "intestine", "lymphoid system")) %>%
            mutate(atlas = consensus_content_name) %>%
            select(ensg_id,
                   consensus_content_name = content_name,
                   norm_exp,
                   atlas)) %>%
    rbind(blood.atlas.max %>%
            select(ensg_id, 
                   consensus_content_name = content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
            mutate(atlas = "blood cells") %>%
            ungroup()) %>% 
    ungroup() %>%
    rbind(select(cell.lines.atlas, 1:3) %>%
            rename(consensus_content_name = celline_name) %>%
            mutate(atlas = "celline")) %>%
    filter(norm_exp >= 1) %>%
    group_by(consensus_content_name, atlas) %>%
    summarise(n_detected_genes = length(unique(ensg_id))) %>% 
    ungroup()
  
  plot.data %>%
    ggplot(aes(atlas, n_detected_genes, fill = atlas, color = atlas)) +
    geom_boxplot(alpha = 0.2) + 
    ggbeeswarm::geom_beeswarm() +
    simple_theme
  ggsave(paste(outpath,'detected_genes_separate_boxplot.pdf',sep='/'), width=8, height=8)
  
  
  plot.data <- 
    all.atlas.max %>%
    select(ensg_id, 
           consensus_content_name, 
           norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
    filter(!consensus_content_name %in% c("brain", "blood", "intestine", "lymphoid system")) %>%
    mutate(atlas = "tissues") %>%
    ungroup() %>%
    rbind(all.atlas %>%
            select(ensg_id, 
                   content_name,
                   consensus_content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression) %>%
            filter(consensus_content_name %in% c("brain", "intestine", "lymphoid system")) %>%
            mutate(atlas = consensus_content_name) %>%
            select(ensg_id,
                   consensus_content_name = content_name,
                   norm_exp,
                   atlas)) %>%
    rbind(blood.atlas.max %>%
            select(ensg_id, 
                   consensus_content_name = content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
            mutate(atlas = "blood cells") %>%
            ungroup()) %>% 
    ungroup() %>%
    rbind(select(cell.lines.atlas, 1:3) %>%
            rename(consensus_content_name = celline_name) %>%
            mutate(atlas = "celline")) %>%
    filter(norm_exp >= 1) %>%
    group_by(atlas) %>%
    summarise(total_detected_genes = length(unique(ensg_id))) %>% 
    ungroup()
  
  plot.data %>%
    ggplot(aes(atlas, total_detected_genes, fill = atlas, color = atlas)) +
    geom_bar(stat = "identity") + 
    simple_theme
  ggsave(paste(outpath,'total_detected_genes_bar.pdf',sep='/'), width=8, height=8)
  
  
  
  ###
  plot.data <- 
    all.atlas.max %>%
    select(ensg_id, 
           consensus_content_name, 
           norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
    filter(!consensus_content_name %in% c("brain", "blood", "intestine", "lymphoid system")) %>%
    mutate(atlas = "single tissues") %>%
    ungroup() %>%
    rbind(all.atlas.max %>%
            select(ensg_id, 
                   consensus_content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
            filter(consensus_content_name %in% c("brain", "blood", "intestine", "lymphoid system")) %>%
            mutate(atlas = "multiple tissues") %>%
            ungroup()) %>%
    rbind(blood.atlas.max %>%
            select(ensg_id, 
                   consensus_content_name = content_name, 
                   norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
            mutate(atlas = "single blood cells") %>%
            ungroup()) %>% 
    ungroup() %>%
    rbind(select(cell.lines.atlas, 1:3) %>%
            rename(consensus_content_name = celline_name) %>%
            mutate(atlas = "cell line")) %>%
    filter(norm_exp >= 1) %>%
    group_by(consensus_content_name, atlas) %>%
    summarise(n_detected_genes = length(unique(ensg_id))) %>% 
    ungroup()
  
  plot.data %>%
    mutate(atlas = factor(atlas, levels = c("single blood cells", "cell line", "single tissues", "multiple tissues"))) %>%
    {ggplot(., aes(atlas, n_detected_genes, fill = atlas, color = atlas)) +
        geom_boxplot(alpha = 0.2) + 
        ggbeeswarm::geom_beeswarm() +
        simple_theme + 
        geom_text(data = subset(., atlas == "multiple tissues" | consensus_content_name == "testis"), 
                  aes(atlas, n_detected_genes, color = atlas, label = consensus_content_name), hjust = 1, nudge_x = -0.1)}
  ggsave(paste(outpath,'detected_genes_separate_boxplot2.pdf',sep='/'), width=8, height=8)
}


#####

make_plots <- function(atlas, atlas.max, atlas.cat, Ex_column, maxEx_column, content_column, consensus_content_column, content_hierarchy = NULL, 
                       content_colors, plots = "all", plot.atlas = c("tissue", "blood", "brain"), plot.order,
                       subatlas_unit = "celltype", global.atlas.category = NULL,
                       sample.atlas = NULL, FACS_markers = NULL, sample_content_column = NULL, sample_Ex_column = NULL, 
                       outpath, prefix) {
  
  atlas.max.wide <- generate_wide(atlas.max, 
                                  ensg_column = 'ensg_id', 
                                  group_column = consensus_content_column, 
                                  max_column = maxEx_column)
  
  atlas.elevated.table <- calc_elevated.table(tb.wide = atlas.max.wide, 
                                              atlas.categories = atlas.cat)
  
  atlas.elevated.summary.table <- calc_elevated.summary.table(atlas.elevated.table)
  
  ## ----- Spearman method cluster ------
  
  if("spearman dendrogram" %in% plots | "all" %in% plots) {
    
    make_spearman_method_dendrogram(all.atlas.tb = atlas, 
                                    Ex_column = Ex_column, 
                                    content_column = content_column, 
                                    named_color_replacement = dataset.colors, 
                                    outpath = outpath, 
                                    prefix = paste0(prefix, "_method_color"))
    
    make_spearman_method_dendrogram(all.atlas.tb = atlas, 
                                    Ex_column = Ex_column, 
                                    content_column = content_column, 
                                    named_color_replacement = tissue.colors, 
                                    outpath = outpath, 
                                    prefix = paste0(prefix, "_tissue_color"))
  }
  
  
  
  
  ## ----- tissue distribution of normalized values -----
  
  if("tissue distribution" %in% plots | "all" %in% plots) {
    make_tissue_distribution_plot(tb.atlas = atlas, 
                                  expr_column = Ex_column,
                                  outpath = outpath,
                                  prefix = prefix)
  }
  
  ## ----- PCA and clustering plots -----
  if("PCA" %in% plots | "all" %in% plots) {
    atlas.max.pca.values <- pca.cal(atlas.max.wide)
    scores <- atlas.max.pca.values[[1]]
    loadings <- 
      atlas.max.pca.values[[2]] %>%
      as.tibble(rownames = "ensg_id") %>%
      mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
    
    make_PCA_plots(scores = scores,
                   loadings = loadings,
                   groups = setNames(rownames(atlas.max.pca.values[[1]]), 
                                     rownames(atlas.max.pca.values[[1]])),
                   groups.color = content_colors,
                   outpath = outpath,
                   prefix = prefix)
  }
  
  if("cluster" %in% plots | "all" %in% plots) {
    make_clustering_plot(tb.wide = atlas.max.wide, 
                         colors = content_colors, 
                         outpath = outpath,
                         prefix = prefix)
  }
  
  ## ----- tissue elevated plot -----
  if("elevated bar" %in% plots | "all" %in% plots) {
    make_elevated_bar_plot(elevated.summary.table = atlas.elevated.summary.table, 
                           outpath = outpath,
                           prefix = prefix)
  }
  
  ## ----- specificity distribution ----
  if("specificity distribution" %in% plots | "all" %in% plots) {
    make_specificity_distribution_plot(atlas.cat = atlas.cat, 
                                     type = "Tissue",
                                     outpath = outpath,
                                     prefix = prefix)
  }
  
  ## ----- chord diagrams ----
  ## chord plot
  if("class chord" %in% plots | "all" %in% plots) {
    make_classification_chord_plot(atlas.cat = atlas.cat,
                                   outpath = outpath,
                                   prefix = prefix)
  }
  
  
  # group enriched chord diagram
  if("group chord" %in% plots | "all" %in% plots) {
    make_chord_group_enriched(atlas.elevated.table, 
                              grid.col = content_colors, 
                              tissue_hierarcy = content_hierarchy,
                              palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
                              outpath = outpath, 
                              prefix = prefix)
  }
  
  ## ----- heatmaps ----
  if("heatmaps" %in% plots | "all" %in% plots) {
    make_heatmap_group_enriched(atlas.elevated.table, 
                                outpath = outpath,
                                prefix = prefix)
    
    make_expression_heatmaps(atlas.max.tb = atlas.max, 
                             atlas.cat = atlas.cat, 
                             maxEx_column = maxEx_column, 
                             tissue_column = consensus_content_column, 
                             ensemblanno.table = ensemblanno.table,
                             proteinclass.table = proteinclass.table, 
                             proteinclass.table_ensg_id_column = "rna.genes", 
                             proteinclass.table_class_column = "proteinclass.vec.single", 
                             outpath = outpath, 
                             prefix = prefix)
    
    # make_expression_heatmaps(atlas.max.tb = atlas.max, 
    #                          atlas.cat = atlas.cat, 
    #                          maxEx_column = maxEx_column, 
    #                          tissue_column = consensus_content_column, 
    #                          ensemblanno.table = ensemblanno.table,
    #                          proteinclass.table = proteinclass.table, 
    #                          proteinclass.table_ensg_id_column = "rna.genes", 
    #                          proteinclass.table_class_column = "proteinclass.vec.single", 
    #                          outpath = outpath, 
    #                          range_scale_x = T,
    #                          prefix = paste(prefix, "range scaled", sep = "_"))
    
    
    
    # make_heatmap_median_expression_levels(elevated.table = atlas.elevated.table,
    #                                       all.atlas.max.tb = atlas.max, 
    #                                       maxEx_column = maxEx_column,
    #                                       tissue_column = content_column,
    #                                       enrichment = c(3),
    #                                       outpath = outpath,
    #                                       prefix = paste(prefix, "group_enriched", sep = "_"))
    # 
    # make_heatmap_median_expression_levels(elevated.table = atlas.elevated.table,
    #                                       all.atlas.max.tb = atlas.max, 
    #                                       maxEx_column = maxEx_column,
    #                                       tissue_column = content_column,
    #                                       enrichment = c(2, 3, 4),
    #                                       outpath = outpath,
    #                                       prefix = paste(prefix, "all_elevated", sep = "_"))
    
    # make_heatmap_expression_levels(elevated.table = atlas.elevated.table,
    #                                all.atlas.max.tb = atlas.max, 
    #                                maxEx_column = maxEx_column,
    #                                tissue_column = content_column,
    #                                enrichment = c(3),
    #                                outpath = outpath,
    #                                prefix = paste(prefix, "group_enriched", sep = "_"))
    # make_heatmap_expression_levels(elevated.table = atlas.elevated.table,
    #                                all.atlas.max.tb = atlas.max, 
    #                                maxEx_column = maxEx_column,
    #                                tissue_column = content_column,
    #                                enrichment = c(2, 3, 4),
    #                                outpath = outpath,
    #                                prefix = paste(prefix, "all_elevated", sep = "_"))
    # make_heatmap_expression_levels(elevated.table = atlas.elevated.table,
    #                                all.atlas.max.tb = atlas.max, 
    #                                maxEx_column = maxEx_column,
    #                                tissue_column = content_column,
    #                                enrichment = c(2),
    #                                outpath = outpath,
    #                                prefix = paste(prefix, "tissue_enriched", sep = "_"))
  }
  
  ## ----- swarm plots ----
  ## swarm plot
  if("swarm expression" %in% plots | "all" %in% plots) {
    # make_swarm_expression_plot(atlas.max = atlas.max, 
    #                            atlas.cat = atlas.cat, 
    #                            maxEx_column = maxEx_column, 
    #                            tissue_column = consensus_content_column,
    #                            outpath = outpath,
    #                            prefix = prefix)
    # 
    # make_swarm_expression_plot(atlas.max = atlas.max, 
    #                            atlas.cat = atlas.cat, 
    #                            maxEx_column = maxEx_column, 
    #                            tissue_column = consensus_content_column, 
    #                            func = swarm_plot_points,
    #                            outpath = outpath,
    #                            prefix = paste0(prefix, "_points"))
    
    make_swarm_expression_plot(atlas.max = atlas.max, 
                               atlas.cat = atlas.cat, 
                               maxEx_column = maxEx_column, 
                               tissue_column = consensus_content_column, 
                               func = swarm_plot_points_top5,
                               outpath = outpath,
                               prefix = paste0(prefix, "_points_top5"))
    
  }
  
  
  ## ----- bar and pie plots ----
  # Number of expressed genes
  if("number detected bar" %in% plots | "all" %in% plots) {
    make_number_detected_genes_barplot(all.atlas.max.tb = atlas.max, 
                                       maxEx_column = maxEx_column,
                                       tissue_column = consensus_content_column,
                                       outpath = outpath,
                                       prefix = prefix)
  }
  
  # Total elevated expression fraction
  # if("NX fraction bar" %in% plots | "all" %in% plots) {
  #   make_elevated_NX_fraction_barplots(atlas.max = atlas.max, 
  #                                      atlas.cat = atlas.cat, 
  #                                      maxEx_column = maxEx_column,
  #                                      tissue_column = consensus_content_column,
  #                                      outpath = outpath, 
  #                                      prefix = prefix)
  # }
  
  if("classification pie" %in% plots | "all" %in% plots) {
    make_classification_pie_chart(atlas.cat = atlas.cat, 
                                  outpath = outpath, 
                                  prefix = prefix)
  }
  
  
  
  # TPM & NX for 100 random genes
  if("TPM NX example genes bar" %in% plots | "all" %in% plots) {
    atlas.max %>%
      filter(ensg_id %in% unique(ensg_id)[1:100]) %>%
      make_gene_expression_barplot(maxEx_columns = c("PTPM" = "max_ptpm", 
                                                     "NX" = "max_nx"),
                                   content_column = consensus_content_column, 
                                   content_color = content_colors)
  }
  
  ## ----- elevated score plots ----
  if("score plots" %in% plots | "all" %in% plots) {
    make_score_expression_scatter(atlas.max.tb = atlas.max, 
                                  atlas.cat = atlas.cat, 
                                  maxEx_column = maxEx_column, 
                                  tissue_column = consensus_content_column, 
                                  ensemblanno.table = ensemblanno.table,
                                  plot.order = plot.order,
                                  outpath = outpath, 
                                  prefix = prefix)
  }
  
  ## ----- Sum TPM and similar plots -----
  if("sum TPM" %in% plots | "all" %in% plots) {
    make_sum_TPM_plot(atlas = atlas, 
                      atlas.cat = atlas.cat, 
                      tissue_column = tissue_column, 
                      method_column = method_column, 
                      outpath = outpath, 
                      prefix = prefix)
    make_sum_TPM_gene_bar_plot(atlas = atlas, 
                               atlas.cat = atlas.cat, 
                               tissue_column = tissue_column, 
                               method_column = method_column, 
                               outpath = outpath, 
                               prefix = prefix)
  }
  
  ## ----- UMAP and tSNE -----
  if("UMAP tSNE" %in% plots | "all" %in% plots) {
    
    make_umap_plot_2(tb = atlas, 
                     ensg_column = "ensg_id", 
                     sample_column = content_column, 
                     group_column = NULL, 
                     Ex_column = sample_Ex_column, 
                     plot_colors = content_colors,
                     mean_expression_filter = 0.1,
                     nudge_x = 1.5, 
                     seed = 42, 
                     tSNE_perplexity = 3,
                     outpath = outpath,
                     prefix = prefix)
  }
  
  # =========== *Subatlas*    ===========
  if("brain" %in% plot.atlas | "blood" %in% plot.atlas) {
    
    if("class comparison chord" %in% plots | "all" %in% plots) {
      # Categories between blood and all atlas
      make_class_comparison_chord(cat1 = atlas.cat, 
                                  cat2 = all.atlas.category,
                                  outpath = outpath, 
                                  prefix = prefix)
    }
    
    if("class organ chord" %in% plots | "all" %in% plots) {
      #Comparison of elevated genes to tissue atlas
      make_elevated_organ_total_chord(cat1 = atlas.cat, 
                                      cat2 = all.atlas.category, 
                                      grid.col = c(content_colors, tissue.colors), 
                                      elevated_cats = c(2,3,4), 
                                      direction = 1, 
                                      cat1_name = subatlas_unit, 
                                      cat2_name = "tissues",
                                      outpath = outpath, 
                                      prefix = paste(prefix, "_tissue_elevated", sep = "_"))
      
      
      make_elevated_organ_total_chord(cat1 = atlas.cat, 
                                      cat2 = all.atlas.category, 
                                      grid.col = c(content_colors, tissue.colors), 
                                      elevated_cats = c(2), 
                                      direction = 1, 
                                      cat1_name = subatlas_unit, 
                                      cat2_name = "tissues",
                                      outpath = outpath, 
                                      prefix = paste(prefix, "_tissue_enriched", sep = "_"))
      
    }
    
  }
  
  if("umap tsne" %in% plots | "all" %in% plots) {
    
    # make_umap_plot(eset,outpath,prefix)
    
  }
  
  
  
  # =========== *Brain altas* =========== 
  
  if("brain" %in% plot.atlas) {
    brain.atlas.max.wide_all_regions <- generate_wide(brain.atlas.max_all_regions, ensg_column='ensg_id',
                                                      group_column='content_name',
                                                      max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
    
    if("spearman dendrogram" %in% plots | "all" %in% plots) {
      
      
      cell.colors <- with(brainregions.table, setNames(subgroup.color, tissue.type))
      make_clustering_plot(tb.wide = brain.atlas.max.wide_all_regions, 
                           colors = cell.colors, 
                           outpath = outpath,
                           prefix = 'brain_all_cells')
    }
    
    if("double donut chord" %in% plots | "all" %in% plots) {
      
      
      make_double_donut_chord(cat1 = all.atlas.category, 
                              cat2 = brain.atlas.category, 
                              regional_class_dictionary = c("low regional specificity",
                                                            "not detected in brain regions",
                                                            "regionally elevated"), 
                              grid.col = c("brain elevated" = "seagreen1", 
                                           "elevated in other tissue" = "skyblue3", 
                                           "low tissue specificity" = "grey40", 
                                           "low regional specificity" = "grey40", 
                                           "regionally elevated" = "seagreen1", 
                                           "not detected in brain regions" = "gray"), 
                              global_name = "brain", 
                              prefix = prefix, 
                              outpath = outpath)
      
    }
    
    
  }
  
  # =========== *Blood altas* =========== 
  
  if("blood" %in% plot.atlas) {
    
    if("blood class tissue expression" %in% plots | "all" %in% plots) {
      
      make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
                               atlas.cat = blood.atlas.category, 
                               maxEx_column = maxEx_column, 
                               tissue_column = "consensus_content_name", 
                               ensemblanno.table = ensemblanno.table,
                               proteinclass.table = proteinclass.table, 
                               proteinclass.table_ensg_id_column = "rna.genes", 
                               proteinclass.table_class_column = "proteinclass.vec.single", 
                               outpath = outpath, 
                               prefix = "blood atlas cat on all atlas")
      
      make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
                               atlas.cat = blood.atlas.category, 
                               maxEx_column = maxEx_column, 
                               tissue_column = "consensus_content_name", 
                               ensemblanno.table = ensemblanno.table,
                               proteinclass.table = proteinclass.table, 
                               proteinclass.table_ensg_id_column = "rna.genes", 
                               proteinclass.table_class_column = "proteinclass.vec.single", 
                               outpath = outpath, 
                               range_scale_x = T,
                               prefix = "blood atlas cat on all atlas range scaled")
      
      make_immunodeficiency_expression_heatmaps(atlas.max.tb = atlas.max,
                                                atlas.cat = atlas.cat,
                                                immunodeficiency.table = immunodeficiency.table,
                                                maxEx_column = maxEx_column,
                                                tissue_column = content_column,
                                                ensemblanno.table = ensemblanno.table,
                                                proteinclass.table = proteinclass.table,
                                                proteinclass.table_ensg_id_column = "rna.genes",
                                                proteinclass.table_class_column = "proteinclass.vec.single",
                                                outpath = outpath,
                                                range_scale_x = F,
                                                prefix = prefix)
      
      make_immunodeficiency_expression_heatmaps(atlas.max.tb = atlas.max,
                                                atlas.cat = atlas.cat,
                                                immunodeficiency.table = immunodeficiency.table,
                                                maxEx_column = maxEx_column,
                                                tissue_column = content_column,
                                                ensemblanno.table = ensemblanno.table,
                                                proteinclass.table = proteinclass.table,
                                                proteinclass.table_ensg_id_column = "rna.genes",
                                                proteinclass.table_class_column = "proteinclass.vec.single",
                                                outpath = outpath,
                                                range_scale_x = F,
                                                prefix = paste(prefix, "range_scaled"))

      
    }
    if("double donut chord" %in% plots | "all" %in% plots) {
      
      
      make_double_donut_chord(cat1 = all.atlas.category, 
                              cat2 = blood.atlas.category, 
                              regional_class_dictionary = c("low blood cell specificity",
                                                            "not detected in blood cells",
                                                            "blood cell elevated"), 
                              grid.col = c("blood elevated" = "dark red", 
                                           "elevated in other tissue" = "skyblue3", 
                                           "low tissue specificity" = "grey40", 
                                           "low blood cell specificity" = "grey40", 
                                           "blood cell elevated" = "darkred", 
                                           "not detected in blood cells" = "gray"), 
                              global_name = "blood", 
                              prefix = prefix, 
                              outpath = outpath)
      
    }
    
    if("sample FACS boxplot" %in% plots | "all" %in% plots) {
      
      make_FACS_boxplot(sample_atlas = sample.atlas, 
                        FACS_markers = FACS_markers, 
                        plot.order = plot.order, 
                        Ex_column = sample_Ex_column, 
                        content_hierarchy = content_hierarchy, 
                        content_column = sample_content_column, 
                        plot_colors = content_colors, 
                        x_title = "TMM normalized PTPM")
      
    }
    
    if("sample UMAP tSNE" %in% plots | "all" %in% plots) {
      
      make_umap_plot_2(tb = sample.atlas, 
                       ensg_column = "ensg_id", 
                       sample_column = "tissue_sample", 
                       group_column = sample_content_column, 
                       Ex_column = sample_Ex_column, 
                       plot_colors = content_colors,
                       mean_expression_filter = 0.1,
                       nudge_x = 1.5, 
                       seed = 42, 
                       tSNE_perplexity = 6,
                       outpath = outpath,
                       prefix = paste(prefix, "samples", sep = "_"))
      
    }
    
     
  }
  
  
  
  
}


#####

make_blood_atlas_paper_plots <- function(atlas, atlas.max, atlas.cat, Ex_column, maxEx_column, content_column, 
                                         consensus_content_column, content_hierarchy = NULL, 
                                         content_colors, plot.order, global.atlas.category = NULL,
                                         sample.atlas = NULL, FACS_markers = NULL, sample_content_column = NULL, sample_Ex_column = NULL, 
                                         plot_theme, 
                                         outpath, prefix) {
  
  atlas.max.wide <- generate_wide(atlas.max, 
                                  ensg_column = 'ensg_id', 
                                  group_column = consensus_content_column, 
                                  max_column = maxEx_column)
  
  atlas.elevated.table <- calc_elevated.table(tb.wide = atlas.max.wide, 
                                              atlas.categories = atlas.cat)
  
  atlas.elevated.summary.table <- calc_elevated.summary.table(atlas.elevated.table)
  
  # 1b
  
  # *** This plot could be smaller but the legend is too big right now
  
  atlas.max %>%
    rename(consensus_content_column = consensus_content_column) %>%
    left_join(ensemblanno.table, by = "ensg_id") %>%
    left_join(content_hierarchy, by = c("consensus_content_column" = "content")) %>%
    filter(gene_name %in% c("CD40", "CD8A", "CTLA4")) %>%
    mutate(consensus_content_column = factor(consensus_content_column, levels = plot.order)) %>% 
    ggplot(aes(consensus_content_column, nx, fill = content_l1)) + 
    geom_bar(stat = "identity") + 
    facet_wrap(~ gene_name, nrow = 3, scales = "free_y") + 
    plot_theme+
    theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top") + 
    ylab("NX") + 
    xlab("") + 
    scale_fill_manual(values = content_colors, name = "")
  
  ggsave(paste(outpath, "Fig 1B.pdf", sep = "/"), width = 5, height = 8) 
  
  # 1c
  
  
}

#####
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

make_umap_plot_2 <- function(tb, ensg_column, sample_column, group_column = NULL, Ex_column, plot_colors, 
                             mean_expression_filter = 0, nudge_x = 1, seed = 42, tSNE_perplexity = 6, outpath, prefix) {
  require(Biobase)
  require(umap)
  require(Rtsne)
  
  genes_passing_filters <-
    tb %>%
    rename(ensg_column = ensg_column,
           Ex_column = Ex_column) %>%
    group_by(ensg_column) %>%
    summarise(mean_expression = mean(Ex_column, na.rm = T)) %>%
    filter(mean_expression >= mean_expression_filter) %$%
    ensg_column
  
  if(is.null(group_column)) {
    tb_ <- 
      tb %>%
      rename(ensg_column = ensg_column,
             Ex_column = Ex_column, 
             sample_column = sample_column) %>%
      filter(ensg_column %in% genes_passing_filters)
  } else {
    tb_ <- 
      tb %>%
      rename(ensg_column = ensg_column,
             Ex_column = Ex_column, 
             group_column = group_column,
             sample_column = sample_column) %>%
      filter(ensg_column %in% genes_passing_filters)
  }
  
  
  exprs <- 
    tb_ %>%
    dplyr::select(sample_column, ensg_column, Ex_column) %>%
    spread(key = ensg_column, value = Ex_column) %>%
    column_to_rownames("sample_column") 
  
  if(is.null(group_column)) {
    pdata <- 
      tb_ %>%
      select(sample_column = sample_column) %>%
      distinct() %>%
      {.[order(match(.$sample_column, rownames(exprs))),]} 
  } else {
    pdata <- 
      tb_ %>%
      select(sample_column = sample_column,
             group_column = group_column) %>%
      distinct() %>%
      {.[order(match(.$sample_column, rownames(exprs))),]} 
  }
  

  stopifnot(all(rownames(exprs) == pdata$sample_column))
  
  # UMAP
  set.seed(seed)
  embedding <- umap(as.matrix(exprs))
  
  
  if(is.null(group_column)) {
    umap.plot.matrix <-
      embedding$layout %>%
      as.tibble() %>%
      mutate(sample_column = as.factor(pdata$sample_column))
  } else {
    umap.plot.matrix <-
      embedding$layout %>%
      as.tibble() %>%
      mutate(group_column = as.factor(pdata$group_column),
             sample_column = as.factor(pdata$sample_column))
  }
  
  
  if(is.null(group_column)) {
    umap.plot.matrix %>%
      {ggplot(., aes(V1, V2, label = sample_column, fill = sample_column)) + 
          geom_point(size = 2.5, 
                     shape = 21,
                     show.legend = F)+
          
          geom_text_repel(aes(color = sample_column),
                          size=4, 
                          show.legend = F)+ 
          xlab('UMAP1')+
          ylab('UMAP2')+
          #labs(title= paste0('umap: ', prefix))+
          simple_theme+ 
          scale_color_manual(values = plot_colors)+
          scale_fill_manual(values = plot_colors)}
  } else {
    umap.plot.matrix %>%
      group_by(., group_column) %>% 
      mutate(mean_V1 = mean(V1), mean_V2 = mean(V2)) %>%
      {ggplot(., aes(V1, V2, label = group_column, fill = group_column)) + 
          geom_point(size = 2.5, 
                     shape = 21,
                     show.legend = F)+
          
          geom_segment(aes(xend = mean_V1, yend = mean_V2, color = group_column),
                       size=0.5,
                       #linetype='dotted',
                       show.legend = F)+
          group_by(., group_column) %>% 
          summarise(V1 = mean(V1), V2 = mean(V2)) %>%
          {geom_label_repel(data = .,
                            size=4, 
                            nudge_x = ifelse(.$V1 > 0, -nudge_x, nudge_x),
                            show.legend = F)}+ 
          xlab('UMAP1')+
          ylab('UMAP2')+
          #labs(title= paste0('umap: ', prefix))+
          simple_theme+ 
          scale_color_manual(values = plot_colors)+
          scale_fill_manual(values = plot_colors)}
  }
  
  
  
  
  ggsave(paste(outpath, paste0(prefix, '_umap.pdf'),sep='/'), width=10, height=10)
  
  
  # tSNE
  exprs_pdata <- 
    pdata %>%
    left_join(exprs %>%
                rownames_to_column("sample_column"), 
              by = "sample_column")
  
  set.seed(seed)
  tsne_out <- Rtsne(exprs, perplexity = tSNE_perplexity, check_duplicates = F)
  
  if(is.null(group_column)) {
    tsne_out$Y %>%
      as.tibble() %>% 
      bind_cols(pdata) %>%
      {ggplot(., aes(V1, V2, label = sample_column, fill = sample_column)) + 
          geom_point(size = 2.5, 
                     shape = 21,
                     show.legend = F)+
          
          geom_text_repel(aes(color = sample_column),
                          size=4, 
                          show.legend = F)+ 
          xlab('tSNE1')+
          ylab('tSNE2')+
          #labs(title= paste0('umap: ', prefix))+
          simple_theme+ 
          scale_color_manual(values = plot_colors)+
          scale_fill_manual(values = plot_colors)}
    
    
  } else {
    tsne_out$Y %>%
      as.tibble() %>% 
      bind_cols(pdata) %>%
      group_by(., group_column) %>% 
      mutate(mean_V1 = mean(V1), mean_V2 = mean(V2)) %>%
      {ggplot(., aes(V1, V2, label = group_column, fill = group_column)) + 
          geom_point(size = 2.5, 
                     shape = 21,
                     show.legend = F)+
          
          geom_segment(aes(xend = mean_V1, yend = mean_V2, color = group_column),
                       size=0.5,
                       #linetype='dotted',
                       show.legend = F)+
          group_by(., group_column) %>% 
          summarise(V1 = mean(V1), V2 = mean(V2)) %>%
          {geom_label_repel(data = .,
                            size=4, 
                            nudge_x = ifelse(.$V1 > 0, -nudge_x, nudge_x),
                            show.legend = F)}+ 
          xlab('tSNE1')+
          ylab('tSNE2')+
          #labs(title= paste0('umap: ', prefix))+
          simple_theme+ 
          scale_color_manual(values = plot_colors)+
          scale_fill_manual(values = plot_colors)}
    
    
  }
  
  
  ggsave(paste(outpath, paste0(prefix, '_tsne.pdf'),sep='/'), width=10, height=10)
  
}


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
# make_swarm_expression_plot <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, outpath, prefix) {
#   plot.data <- 
#     atlas.max %>%
#     ungroup() %>%
#     mutate(expression = eval(parse(text = maxEx_column)),
#            Grouping = eval(parse(text = tissue_column))) %>%
#     dplyr::select(Grouping, ensg_id, expression) %>%
#     #Plot 1 % highest
#     filter(expression >= quantile(expression, probs = 0.99)) %>%
#     group_by(Grouping) %>%
#     # Write gene name if highest 0.5 % per tissue and/or highest 1 % in total
#     mutate(highest = expression >= quantile(expression, probs = 0.995) | rank(expression) >= (length(expression) - 8)) %>%
#     ungroup() %>%
#     mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98))) %>%
#     left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`), by = "ensg_id") %>%
#     left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
#     left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
#     mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))
# 
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, expression, label = gene_name, color = elevated.category)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = elevated.cat.cols)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_abundance_jitter_1.pdf'),sep='/'), width=15, height=10)
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = protein.class.palette)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_score_jitter_2.pdf'),sep='/'), width=15, height=10)
#   
#   plot.data <- 
#     atlas.max %>%
#     ungroup() %>%
#     mutate(expression = eval(parse(text = maxEx_column)),
#            Grouping = eval(parse(text = tissue_column))) %>%
#     dplyr::select(Grouping, ensg_id, expression) %>%
#     left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
#     mutate(score = as.numeric(`tissue/group specific score`)) %>%
#     # tissue enriched
#     filter(elevated.category == "tissue enriched" & Grouping == `enriched tissues`) %>%
#     group_by(Grouping) %>%
#     # Write gene name if highest 2 % per tissue and/or highest 1 % in total
#     mutate(highest = expression >= quantile(expression, probs = 0.99) | rank(expression) >= (length(expression) - 10),
#            highest_score = score >= quantile(score, probs = 0.99) | rank(score) >= (length(score) - 10)) %>%
#     ungroup() %>%
#     mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98)),
#            highest_score = ifelse(highest_score, T, score >= quantile(score, probs = 0.98))) %>%
#     
#     left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
#     left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
# 
#     mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single)) 
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = protein.class.palette)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_tissue_enriched_jitter.pdf'),sep='/'), width=15, height=10)
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, score, label = gene_name, color = gene_class)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest_score,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = protein.class.palette)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_score_enriched_jitter.pdf'),sep='/'), width=15, height=10)
#   
#   plot.data <- 
#     atlas.max %>%
#     ungroup() %>%
#     mutate(expression = eval(parse(text = maxEx_column)),
#            Grouping = eval(parse(text = tissue_column))) %>%
#     dplyr::select(Grouping, ensg_id, expression) %>%
#     left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
#     mutate(score = as.numeric(`tissue/group specific score`)) %>%
#     # tissue enriched
#     filter(elevated.category %in% c("tissue enriched", "tissue enhanced", "group enriched")) %>%
#     filter(expression >= quantile(expression, probs = 0.99)) %>%
#     group_by(Grouping) %>%
#     #Plot 1 % highest
#     #filter(expression >= quantile(expression, probs = 0.9)) %>%
#     # Write gene name if highest 2 % per tissue and/or highest 1 % in total
#     mutate(highest = expression >= quantile(expression, probs = 0.99) | rank(expression) >= (length(expression) - 10),
#            highest_score = score >= quantile(score, probs = 0.99) | rank(score) >= (length(score) - 10)) %>%
#     ungroup() %>%
#     mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98)),
#            highest_score = ifelse(highest_score, T, score >= quantile(score, probs = 0.98))) %>%
#     
#     left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
#     left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
#     
#     mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single)) 
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, expression, label = gene_name, color = gene_class)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = protein.class.palette)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_tissue_elevated_jitter.pdf'),sep='/'), width=15, height=10)
#   
#   plot.data %>%
#   {ggplot(., aes(Grouping, score, label = gene_name, color = gene_class)) +
#       geom_jitter(alpha = 0.4, size = 1, width = 0.1)+
#       geom_text_repel(data = .[.$highest_score,], size = 2)+
#       simple_theme+
#       scale_color_manual(values = protein.class.palette)+
#       scale_y_log10()+
#       theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))}
#   ggsave(paste(outpath, paste0(prefix, '_high_score_elevated_jitter.pdf'),sep='/'), width=15, height=10)
# }

swarm_circle_plot <- function(plot.df, color_palette, color_column, legend_name, y_column) {
  plot.df <- 
    plot.df %>%
    rename(y = y_column)
  
  sectors <- levels(plot.df$Grouping)
  n_sectors <- length(sectors)
  text_angle <- 
    seq(360 - 360 / (n_sectors * 2), by = - 360 / n_sectors, length.out = n_sectors) 
  text_angle_nice_flip <- 
    text_angle %>%
    {ifelse(. > 90 & . < 270, . + 180, .)}
  
  ymax <- max(plot.df$y)
  ymin <- min(plot.df$y)
  
  
  magnitude <- round(log10(c(ymin, ymax)))
  
  break_n <- 4
  breaks <- 
    10^seq(log10(ceiling(ymin / (10^magnitude[1])) * (10^magnitude[1])), 
           log10(ceiling(ymax / (10^magnitude[2])) * (10^magnitude[2])), length.out = break_n) %>%
    round(-floor(log10(.)))
  
  xaxis_y <- 1.5*max(breaks) - min(breaks)
  inner_y <- breaks[1]/2
    
    
  segment_color <- "gray45"
  
  plot.df %>% 
    {ggplot(., aes(sector_index_swarm, y, label = gene_name, color = eval(parse(text = color_column)))) +
        # Dividers
        annotate(geom = "segment", 
                 x = 0:n_sectors + 0.5, 
                 xend = 0:n_sectors + 0.5,
                 y = inner_y, 
                 yend = xaxis_y, 
                 color = segment_color)+
        # y-axis breaks
        annotate(geom = "segment", 
                 x = rep(0:n_sectors + 0.5, each = break_n) - c(rep(0, break_n), 0.03 * rep(break_n:1, n_sectors)),
                 xend = rep(0:n_sectors + 0.5, each = break_n) + c(0.03 * rep(break_n:1, n_sectors), rep(0, break_n)),
                 y = rep(breaks, n_sectors + 1), 
                 yend = rep(breaks, n_sectors + 1), 
                 color = segment_color)+
        # Inner circle
        annotate(geom = "segment", 
                 x = seq(0.5, n_sectors - 1 + 0.5), 
                 xend = seq(1.5, n_sectors + 0.5),
                 y = inner_y, 
                 yend = inner_y, 
                 color = segment_color)+
        # Outer circle
        annotate(geom = "segment", 
                 x = seq(0.5, n_sectors - 1 + 0.5), 
                 xend = seq(1.5, n_sectors + 0.5),
                 y = xaxis_y, 
                 yend = xaxis_y,
                 color = segment_color)+
        
        geom_text(aes(size = range_scale(y, span = c(1, 4))), alpha = 0.8)+
        
        scale_color_manual(values = color_palette, name = legend_name)+
        scale_y_log10(expand = c(0.1, 0), breaks = breaks)+
        scale_x_continuous(limits = c(0.5, n_sectors + 0.5), expand = c(0,0)) + 
        scale_size_continuous(guide = F)+
        ylab(y_column) + 
        
        coord_polar()+
        
        annotate(geom = "text", 
                 x = 1:n_sectors, 
                 y = xaxis_y, 
                 label = levels(.$Grouping),
                 angle = text_angle_nice_flip,
                 vjust = ifelse(text_angle > 90 & text_angle < 270, 1.5, -1),
                 size = 4)+
        
        theme(axis.text.x = element_blank(), 
              axis.title.x = element_blank(), 
              panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line.y = element_line())}
}

make_swarm_expression_circle_plot <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, plot.order = NULL, outpath, prefix) {
  plot.data.unfilt <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column)))
  if(is.null(plot.order)) {
    plot.data.unfilt <- 
      plot.data.unfilt %>%
      mutate(Grouping = as.factor(Grouping))
  } else {
    plot.data.unfilt <- 
      plot.data.unfilt %>%
      mutate(Grouping = factor(Grouping, levels = plot.order))
  }
  
  plot.data.unfilt <- 
    plot.data.unfilt %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
    left_join(dplyr::select(ensemblanno.table, ensg_id, gene_name, gene_description, ncbi_gene_summary, chr_name) , by = "ensg_id") %>%
    left_join(dplyr::select(proteinclass.table, rna.genes, proteinclass.vec.single), by = c("ensg_id"="rna.genes")) %>%
    mutate(score = `tissue/group specific score`) 
  
  
  #---- Elevated genes
  
  plot.data <- 
    plot.data.unfilt %>%
    filter(elevated.category %in% c("tissue enriched", "group enriched", "tissue enhanced") &
             mapply(paste0("(^|, )", Grouping, "(, |$)"), 
                    `enriched tissues`, FUN = function(x,y) grepl(x, y))) %>%
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))
  
  sectors <- levels(plot.data$Grouping)
  
  plot.data <- 
    plot.data %>%
    mutate(sector_index = match(Grouping, sectors),
           sector_index_swarm = sector_index + runif(length(sector_index), min=-0.35, max=0.35)) 
  
  # --All elevated genes
  # Expression
  plot.data %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_all_elevated_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_all_elevated_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  # Score
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_all_elevated_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_all_elevated_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  # --Top 20 % elevated genes
  
  plot.data <- 
    plot.data %>%
    group_by(Grouping) %>%
    filter(expression >= quantile(expression, probs = 0.8)) %>%
    ungroup()
  
  # Expression
  plot.data %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_elevated_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_elevated_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  # Score
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_top_elevated_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_top_elevated_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  #---- Tissue enriched genes
  
  plot.data <- 
    plot.data.unfilt %>%
    filter(elevated.category %in% c("tissue enriched") &
             Grouping == `enriched tissues`) %>%
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))
  
  sectors <- levels(plot.data$Grouping)
  
  plot.data <- 
    plot.data %>%
    mutate(sector_index = match(Grouping, sectors),
           sector_index_swarm = sector_index + runif(length(sector_index), min=-0.35, max=0.35)) 
  
  # Expression
  plot.data %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_tissue_enriched_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_tissue_enriched_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  # Score
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_tissue_enriched_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_tissue_enriched_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  
  #---- All genes
  
  plot.data <- 
    plot.data.unfilt %>%
    mutate(gene_class = ifelse(is.na(proteinclass.vec.single), "other", proteinclass.vec.single))
  
  sectors <- levels(plot.data$Grouping)
  
  plot.data <- 
    plot.data %>%
    mutate(sector_index = match(Grouping, sectors),
           sector_index_swarm = sector_index + runif(length(sector_index), min=-0.35, max=0.35)) 
  
  
  # --Top 1 % genes
  
  plot.data <- 
    plot.data %>%
    group_by(Grouping) %>%
    filter(expression >= quantile(expression, probs = 0.99)) %>%
    ungroup()
  
  # Expression
  plot.data %>%
    swarm_circle_plot(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_all_genes_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    swarm_circle_plot(color_palette = protein.class.palette, 
                      color_column = "gene_class",
                      legend_name = "Protein class", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_all_genes_category_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
}

swarm_plot <- function(plot.df, color_palette, color_column, legend_name, y_column) {
  plot.df %>%
    rename(y = y_column) %>%
    {ggplot(., aes(Grouping, y, label = display_name, color = eval(parse(text = color_column)))) +
    #geom_text_repel(direction = "x", segment.color = NA)+
        
    ggbeeswarm::geom_beeswarm(data = subset(., eval(parse(text = color_column)) == ""), color = "black") + 
    geom_text(data = subset(., eval(parse(text = color_column)) != ""), fontface = "bold")+
        
    scale_color_manual(values = color_palette, name = legend_name)+
    scale_y_log10()+
    ylab(y_column)+
    scale_size_continuous(guide = F)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
          axis.title.x = element_blank(), 
          panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.line.y = element_line(),
          axis.line.x = element_line())}
}

swarm_plot_points <- function(plot.df, color_palette, color_column, legend_name, y_column) {
  plot.df %>%
    rename(y = y_column) %>%
    {ggplot(., aes(Grouping, y, label = display_name, color = eval(parse(text = color_column)))) +
        ggbeeswarm::geom_beeswarm(cex = 0.5)+
        ggrepel::geom_text_repel(data = {
          
          genes <- 
            group_by(., Grouping) %>%
            mutate(rank = - (rank(y) - max(rank(y)) - 1)) %>%
            filter(rank <= 5) %$%
            display_name
          
          filter(., display_name %in% genes)},
          
                                 fontface = "bold")+
        
        scale_color_manual(values = color_palette, name = legend_name)+
        scale_y_log10()+
        ylab(y_column)+
        scale_size_continuous(guide = F)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
              axis.title.x = element_blank(), 
              panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line.y = element_line(),
              axis.line.x = element_line())}
}


swarm_plot_points_top5 <- function(plot.df, color_palette, color_column, legend_name, y_column) {
  plot.df %>%
    rename(y = y_column) %>%
    {ggplot(., aes(Grouping, y, label = display_name, color = eval(parse(text = color_column)))) +
        ggbeeswarm::geom_beeswarm(cex = 0.5)+
        ggrepel::geom_text_repel(data = {
          
          group_by(., Grouping) %>%
            mutate(rank = - (rank(y) - max(rank(y)) - 1)) %>%
            filter(rank <= 5) },
          
          fontface = "bold")+
        
        scale_color_manual(values = color_palette, name = legend_name)+
        scale_y_log10()+
        ylab(y_column)+
        scale_size_continuous(guide = F)+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
              axis.title.x = element_blank(), 
              panel.grid = element_blank(), 
              panel.background = element_blank(), 
              axis.line.y = element_line(),
              axis.line.x = element_line())}
}
make_swarm_expression_plot <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, plot.order = NULL, func = swarm_plot, outpath, prefix) {
  plot.data.unfilt <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column)))
  if(is.null(plot.order)) {
    plot.data.unfilt <- 
      plot.data.unfilt %>%
      mutate(Grouping = as.factor(Grouping))
  } else {
    plot.data.unfilt <- 
      plot.data.unfilt %>%
      mutate(Grouping = factor(Grouping, levels = plot.order))
  }
  
  plot.data.unfilt <- 
    plot.data.unfilt %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    left_join(dplyr::select(atlas.cat, ensg_id, express.category.2, elevated.category, `enriched tissues`, `tissue/group specific score`), by = "ensg_id") %>%
    left_join(proteinlocalization.table, by = "ensg_id") %>%
    mutate(score = `tissue/group specific score`,
           predicted_localization_class_simple = case_when(predicted_localization_class == "intracellular" ~ "",
                                                           predicted_localization_class == "intracellular and membrane isoforms" ~ "membrane",
                                                           predicted_localization_class == "intracellular and secreted isoforms" ~ "secreted",
                                                           predicted_localization_class == "intracellular, membrane, secreted isoforms" ~ "membrane and secreted isoforms",
                                                           T ~ predicted_localization_class),
           predicted_localization_class_simple = case_when(predicted_secreted ~ "secreted",
                                                           predicted_membrane ~ "membrane",
                                                           predicted_intracellular ~ "intracellular",
                                                           T ~ ""))
  
  
  #---- Elevated genes
  
  plot.data <- 
    plot.data.unfilt %>%
    filter(elevated.category %in% c("tissue enriched", "group enriched", "tissue enhanced") &
             mapply(paste0("(^|, )", Grouping, "(, |$)"), 
                    `enriched tissues`, FUN = function(x,y) grepl(x, y))) %>%
    {.[order(.$expression, decreasing = T),][1:1000,]}
  #####################################
  
  # --All elevated genes
  # Expression
  plot.data %>%
    func(color_palette = elevated.cat.cols, 
               color_column = "elevated.category",
               legend_name = "Elevated category", 
               y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_all_elevated_jitter.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    func(color_palette = protein.localization.palette, 
               color_column = "predicted_localization_class_simple",
               legend_name = "Protein location", 
               y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_all_elevated_category_jitter.pdf'),sep='/'), width=15, height=10)
  
  # Score
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    func(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_all_elevated_jitter.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    func(color_palette = protein.localization.palette, 
                      color_column = "predicted_localization_class_simple",
                      legend_name = "Protein location", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_all_elevated_category_jitter.pdf'),sep='/'), width=15, height=10)

  
  #---- Tissue enriched genes
  
  plot.data <- 
    plot.data.unfilt %>%
    filter(elevated.category %in% c("tissue enriched") &
             Grouping == `enriched tissues`) 
  
  # Expression
  plot.data %>%
    func(color_palette = protein.localization.palette, 
                      color_column = "predicted_localization_class_simple",
                      legend_name = "Protein location", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_tissue_enriched_category_jitter.pdf'),sep='/'), width=15, height=10)
  
  # Score
  plot.data %>%
    filter(score != "") %>%
    mutate(score = as.numeric(score)) %>%
    func(color_palette = protein.localization.palette, 
                      color_column = "predicted_localization_class_simple",
                      legend_name = "Protein location", 
                      y_column = "score")
  ggsave(paste(outpath, paste0(prefix, '_score_tissue_enriched_category_jitter.pdf'),sep='/'), width=15, height=10)
  
  
  #---- All genes

  # --Top 1000 genes
  
  plot.data <- 
    plot.data.unfilt %>%
    {.[order(.$expression, decreasing = T),][1:1000,]}
  
  # Expression
  plot.data %>%
    func(color_palette = elevated.cat.cols, 
                      color_column = "elevated.category",
                      legend_name = "Elevated category", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_all_genes_jitter_circle.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    func(color_palette = protein.localization.palette, 
                      color_column = "predicted_localization_class_simple",
                      legend_name = "Protein location", 
                      y_column = "expression")
  ggsave(paste(outpath, paste0(prefix, '_high_top_all_genes_category_jitter.pdf'),sep='/'), width=15, height=10)
  
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
  
  for(filt in c("excat", "elcat")) {
    pdf(paste(outpath, paste0(prefix, " classification ", filt, ".pdf"), sep = "/"), 
        width = 10, height = 10, useDingbats = F)
    df %>%
      filter(sapply(transfer, FUN = function(x) grepl(x, filt))) %$%
      chord_classification_clockwise(from, to, sizes, chord_col, plot.group[1:(length(plot.group)/2)], plot.order, line_expansion = 10000)
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

# make_heatmap_group_enriched_circle <- function(elevated.table, outpath, prefix) {
#   
#   group_enriched_number <- 
#     elevated.table %>%
#     as_tibble(., rownames = "ensg_id") %>%
#     gather(key = "content", "classification", -ensg_id) %>%
#     # Only include group enriched
#     filter(classification %in% c(3)) %>%
#     {mapply(unique(.$content), 
#             FUN = function(cell1) mapply(unique(.$content), 
#                                          FUN = function(cell2) length(intersect(filter(., content == cell1)$ensg_id, 
#                                                                                 filter(., content == cell2)$ensg_id))))} %>%
#     as_tibble(., rownames = "from") %>%
#     gather(key = "to", value = "number of genes", -from) %>%
#     left_join(ungroup(.) %>%
#                 group_by(from) %>%
#                 summarise(total = sum(`number of genes`),
#                           Min = min(`number of genes`),
#                           Max = max(`number of genes`)),
#               by = "from") %>%
#     left_join(ungroup(.) %>%
#                 group_by(to) %>%
#                 summarise(total = sum(`number of genes`),
#                           Min = min(`number of genes`),
#                           Max = max(`number of genes`)),
#               by = "to") %>%
#     mutate(Min = mapply(Min.x, Min.y, FUN = function(a, b) min(c(a, b))),
#            Max = mapply(Max.x, Max.y, FUN = function(a, b) max(c(a, b))),
#            n = (`number of genes` - Min)/ (Max - Min),
#            from.sum = filter(., to == from)$`number of genes`[match(from, filter(., to == from)$from)],
#            to.sum = filter(., to == from)$`number of genes`[match(to, filter(., to == from)$from)],
#            frac.max = mapply(`number of genes`, from.sum, to.sum, FUN = function(n,a,b) max(c(n/a, n/b)))) 
#   
#   GEN_dendrogram <- 
#     group_enriched_number %>%
#     select(from, to, frac.max) %>%
#     spread(key = to, value = frac.max) %>%
#     column_to_rownames("from") %>%
#     as.matrix() %>%
#     cor(method="spearman", use="pairwise.complete.obs")  %>%
#     {1 - .} %>%
#     as.dist() %>%
#     hclust("average")  %>%
#     as.dendrogram()
#   
#   dend_segments <- 
#     ggdendro::dendro_data(GEN_dendrogram)$segments
#     
#   group_enriched_number %>%
#     mutate(from = factor(from, levels = unique(from)[order.dendrogram(GEN_dendrogram)]),
#            to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
#     {nfactors  <- length(levels(.$from));
#     n_sectors <- nfactors + 3
#     text_angle <- function(x) - ((x + 3) * 360/n_sectors - 360/(n_sectors*2))
#     ggplot(.) +
#         geom_tile(aes(from, to, fill = frac.max), alpha = 0.8) + 
#         geom_segment(data = dend_segments, aes(x = x, y = -10*yend, xend = xend, yend = -10*y))+
#         annotate(geom = "text", x = 1:nfactors, y = n_sectors, label = levels(.$from), 
#                  angle = text_angle(1:nfactors),
#                  size = 3)+
#         annotate(geom = "text", x = -2.5, y =1:nfactors, label = levels(.$from), hjust = 0,
#                  size = 3)+
#         annotate("point", x = -2.5, y = n_sectors, color = NA)+
#         #coord_fixed() +
#         coord_polar()+
#         #simple_theme + 
#         theme_minimal()+
#         theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
#         scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
#   
#   ggsave(paste(outpath, paste0(prefix, '_group_enriched_heatmap_circle.pdf'),sep='/'), width=15, height=10)
# }



# make_heatmap_group_and_enhanced_expression_levels_circle <- 
#   function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix) {
#     
#     group_enriched_genes <- 
#       elevated.table %>%
#       as_tibble(., rownames = "ensg_id") %>%
#       gather(key = "content", "classification", -ensg_id) %>%
#       # Only include group enriched
#       filter(classification %in% c(3, 4))
#     
#     group_enriched_expression <- 
#       group_enriched_genes %>%
#       {mapply(unique(.$content), 
#               FUN = function(cell1) mapply(unique(.$content), 
#                                            FUN = function(cell2) {
#                                              genes <- 
#                                                filter(., content == cell1)$ensg_id
#                                              all.atlas.max.tb %>%
#                                                filter(eval(parse(text = tissue_column)) == cell2 & 
#                                                         ensg_id %in% genes) %$%
#                                                median(eval(parse(text = maxEx_column)), na.rm = T)
#                                            }))} %>%
#       as_tibble(., rownames = "from") %>%
#       gather(key = "to", value = "median expression", -from) 
#     
#     GEN_dendrogram <- 
#       group_enriched_expression %>%
#       select(from, to, `median expression`) %>%
#       spread(key = to, value = `median expression`) %>%
#       column_to_rownames("from") %>%
#       as.matrix() %>%
#       cor(method="spearman", use="pairwise.complete.obs")  %>%
#       {1 - .} %>%
#       as.dist() %>%
#       hclust("average")  %>%
#       as.dendrogram()
#     
#     dend_segments <- 
#       ggdendro::dendro_data(GEN_dendrogram)$segments
#     
#     group_enriched_expression %>%
#       mutate(from = factor(from, levels = unique(from)[order.dendrogram(GEN_dendrogram)]),
#              to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
#              {nfactors  <- length(levels(.$from));
#              n_sectors <- nfactors + 3
#              text_angle <- function(x) - ((x + 3) * 360/n_sectors - 360/(n_sectors*2))
#              ggplot(.) +
#                geom_tile(aes(from, to, fill = `median expression`), alpha = 0.8) + 
#                geom_segment(data = dend_segments, aes(x = x, y = -10*yend, xend = xend, yend = -10*y))+
#                annotate(geom = "text", x = 1:nfactors, y = n_sectors, label = levels(.$from), 
#                         angle = text_angle(1:nfactors),
#                         size = 3)+
#                annotate(geom = "text", x = -2.5, y =1:nfactors, label = levels(.$from), hjust = 0,
#                         size = 3)+
#                annotate("point", x = -2.5, y = n_sectors, color = NA)+
#                #coord_fixed() +
#                coord_polar()+
#                #simple_theme + 
#                theme_minimal()+
#                theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
#                scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
#     
#     ggsave(paste(outpath, paste0(prefix, '_group_and_enhanced_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
#   }

# make_heatmap_group_enriched_expression_levels_circle <- 
#   function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix, y_dendrogram = F) {
#     
#     group_enriched_genes <- 
#       elevated.table %>%
#       as_tibble(., rownames = "ensg_id") %>%
#       gather(key = "content", "classification", -ensg_id) %>%
#       # Only include group enriched
#       filter(classification %in% c(3))
#     
#     group_enriched_expression <- 
#       group_enriched_genes %>%
#       right_join(filter(all.atlas.max.tb, ensg_id %in% group_enriched_genes$ensg_id), by = c("ensg_id", "content" = tissue_column)) %$%
#       tibble(from = ensg_id,
#              to = content,
#              expression = log10(eval(parse(text = maxEx_column)) + 1)) 
#       
#     
#     GEN_dendrogram <- 
#       group_enriched_expression %>%
#       select(from, to, expression) %>%
#       spread(key = to, value = expression) %>%
#       column_to_rownames("from") %>%
#       as.matrix() %>%
#       cor(method="spearman", use="pairwise.complete.obs")  %>%
#       {1 - .} %>%
#       as.dist() %>%
#       hclust("average")  %>%
#       as.dendrogram()
#     
#     
#     
#     dend_segments <- 
#       ggdendro::dendro_data(GEN_dendrogram)$segments
#     
#     
#     
#     if(y_dendrogram) {
#       Genes_dendrogram <- 
#         group_enriched_expression %>%
#         select(from, to, expression) %>%
#         spread(key = from, value = expression) %>%
#         column_to_rownames("to") %>%
#         as.matrix() %>%
#         cor(method="spearman", use="pairwise.complete.obs")  %>%
#         {1 - .} %>%
#         as.dist() %>%
#         hclust("average")  %>%
#         as.dendrogram()
#       genes_dend_segments <- 
#         ggdendro::dendro_data(Genes_dendrogram)$segments
#     }
#     
#     x_margin = 2.5
#     
#     g <-
#       group_enriched_expression %>%
#       #filter(from %in% .$from[1:100]) %>%
#       mutate(from = factor(from, levels = unique(from)[order.dendrogram(Genes_dendrogram)]),
#              to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
#              {nfactors  <- length(levels(.$to));
#              n_sectors <- round((1 + 1/6)*nfactors)
#              y_sectors <- length(levels(.$from))
#              
# 
#              text_angle <- function(x) - (round((x + (1/6)*nfactors)) * 360/n_sectors - 360/(n_sectors*2))
#              ggplot(.) +
#                geom_tile(aes(to, from, fill = expression), alpha = 0.8) + 
#                geom_segment(data = dend_segments, aes(x = x, y = -0.5*y_sectors*yend, 
#                                                       xend = xend, yend = -0.5*y_sectors*y))+
#                
#                annotate(geom = "text", x = 1:nfactors, y = 1.15*y_sectors, label = levels(.$to),
#                         angle = text_angle(1:nfactors),
#                         size = 3)+
#                annotate("point", x = -x_margin, y = n_sectors, color = NA)+
#                #coord_fixed() +
#                coord_polar()+
#                #simple_theme + 
#                theme_minimal()+
#                theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
#                scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
#     
#     if(y_dendrogram){
#       genes_dend_segments_scaled <-
#       {max_y = max(c(genes_dend_segments$y,
#                      genes_dend_segments$yend))
#       min_y = min(c(genes_dend_segments$y,
#                     genes_dend_segments$yend))
#       genes_dend_segments %>%
#         mutate(y = ((y - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin),
#                yend = ((yend - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin))
#       }
#       g <- 
#         g + 
#         geom_segment(data = genes_dend_segments_scaled, aes(x = y, y = x,
#                                                             xend = yend, yend = xend))
#       
#     }
#     
#     
#     ggsave(plot = g, paste(outpath, paste0(prefix, '_group_enriched_gene_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
#   }

# make_heatmap_all_elevated_expression_levels_circle <- 
#   function(elevated.table, all.atlas.max.tb, maxEx_column, tissue_column, outpath, prefix, y_dendrogram = F) {
#     
#     group_enriched_genes <- 
#       elevated.table %>%
#       as_tibble(., rownames = "ensg_id") %>%
#       gather(key = "content", "classification", -ensg_id) %>%
#       # Only include group enriched
#       filter(classification %in% c(2, 3, 4))
#     
#     group_enriched_expression <- 
#       group_enriched_genes %>%
#       right_join(filter(all.atlas.max.tb, ensg_id %in% group_enriched_genes$ensg_id), by = c("ensg_id", "content" = tissue_column)) %$%
#       tibble(from = ensg_id,
#              to = content,
#              expression = log10(eval(parse(text = maxEx_column)) + 1)) 
#     
#     
#     GEN_dendrogram <- 
#       group_enriched_expression %>%
#       select(from, to, expression) %>%
#       spread(key = to, value = expression) %>%
#       column_to_rownames("from") %>%
#       as.matrix() %>%
#       cor(method="spearman", use="pairwise.complete.obs")  %>%
#       {1 - .} %>%
#       as.dist() %>%
#       hclust("average")  %>%
#       as.dendrogram()
#     
#     
#     
#     dend_segments <- 
#       ggdendro::dendro_data(GEN_dendrogram)$segments
#     
#     
#     
#     Genes_dendrogram <- 
#       group_enriched_expression %>%
#       select(from, to, expression) %>%
#       spread(key = from, value = expression) %>%
#       column_to_rownames("to") %>%
#       as.matrix() %>%
#       cor(method="spearman", use="pairwise.complete.obs")  %>%
#       {1 - .} %>%
#       as.dist() %>%
#       hclust("average")  %>%
#       as.dendrogram()
#     genes_dend_segments <- 
#       ggdendro::dendro_data(Genes_dendrogram)$segments
#     
#     
#     x_margin = 2.5
#     
#     g <-
#       group_enriched_expression %>%
#       #filter(from %in% .$from[1:100]) %>%
#       mutate(from = factor(from, levels = unique(from)[order.dendrogram(Genes_dendrogram)]),
#              to = factor(to, levels = unique(to)[order.dendrogram(GEN_dendrogram)])) %>%
#              {nfactors  <- length(levels(.$to));
#              n_sectors <- round((1 + 1/6)*nfactors)
#              y_sectors <- length(levels(.$from))
#              
#              
#              text_angle <- function(x) - (round((x + (1/6)*nfactors)) * 360/n_sectors - 360/(n_sectors*2))
#              ggplot(.) +
#                geom_tile(aes(to, from, fill = expression), alpha = 0.8) + 
#                geom_segment(data = dend_segments, aes(x = x, y = -0.5*y_sectors*yend, 
#                                                       xend = xend, yend = -0.5*y_sectors*y))+
#                
#                annotate(geom = "text", x = 1:nfactors, y = 1.15*y_sectors, label = levels(.$to),
#                         angle = text_angle(1:nfactors),
#                         size = 3)+
#                annotate("point", x = -x_margin, y = n_sectors, color = NA)+
#                #coord_fixed() +
#                coord_polar()+
#                #simple_theme + 
#                theme_minimal()+
#                theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank()) + 
#                scale_fill_gradientn(colors = c("yellow", "orangered", "#800026"))} 
#     
#     if(y_dendrogram){
#       genes_dend_segments_scaled <-
#       {max_y = max(c(genes_dend_segments$y,
#                      genes_dend_segments$yend))
#       min_y = min(c(genes_dend_segments$y,
#                     genes_dend_segments$yend))
#       genes_dend_segments %>%
#         mutate(y = ((y - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin),
#                yend = ((yend - min_y)/(max_y - min_y))*(0 - (-x_margin)) + (-x_margin))
#       }
#       g <- 
#         g + 
#         geom_segment(data = genes_dend_segments_scaled, aes(x = y, y = x,
#                                                             xend = yend, yend = xend))
#       
#     }
#     
#     
#     ggsave(plot = g, paste(outpath, paste0(prefix, '_all_elevated_gene_expression_heatmap_circle.pdf'),sep='/'), width=15, height=10)
#   }

make_dendrogram <- function(var1, var2, value, method = "spearman"){
  if(method == "spearman") {
    data.frame(var1, var2, value) %>%
      spread(key = var1, value = value) %>%
      column_to_rownames("var2") %>% 
      as.matrix() %>%
      cor(method="spearman", use="pairwise.complete.obs")  %>%
      {1 - .} %>%
      as.dist() %>%
      hclust("average")  %>%
      as.dendrogram() 
  } else if(method == "euclidean") {
    data.frame(var1, var2, value) %>%
      spread(key = var1, value = value) %>%
      column_to_rownames("var2") %>% 
      as.matrix() %>% 
      t() %>% 
      na.omit() %>%
      dist() %>%
      hclust("average")  %>%
      as.dendrogram() 
  }
  
}

get_dendrogram_segments <- function(dendrogram) {
  dendrogram %>%
    ggdendro::dendro_data() %$%
    left_join(segments, labels, by = c("xend" = "x", "yend" = "y"))
}

range_scale_manual <- function(x, xmax, xmin, span, dodge = 0) { 
  ((x - xmin)/(xmax - xmin))*
    (span[2] - span[1]) + span[1] + 
    dodge
}
range_scale <- function(x, span) { 
  xmin = min(x)
  xmax = max(x)
  ((x - xmin)/(xmax - xmin))*
    (span[2] - span[1]) + span[1]
}
ggdendroheat <- function(x, y, value, fill_factor = NA, show.legend = T, xdendrogram = T, ydendrogram = T,
                         x_margin = 0.2, y_margin = 0.2, xaxis_order = NULL, yaxis_order = NULL, 
                         xdendrogram_color = NULL, ydendrogram_color = NULL,
                         range_scale_x = F, range_scale_y = F, x_clustering_method = "spearman", y_clustering_method = "spearman",
                         under_lim_tile_color = NA, under_lim = 1){
  
  x_dendrogram <- 
    make_dendrogram(x, y, value, method = x_clustering_method)
  x_dendrogram_segments <- get_dendrogram_segments(x_dendrogram)
  
  y_dendrogram <- 
    make_dendrogram(y, x, value, method = y_clustering_method)
  y_dendrogram_segments <- get_dendrogram_segments(y_dendrogram)
  
  g <- 
    tibble(x, y, value, 
           under_limit = value < under_lim,
           fill_factor = fill_factor) 
  
  if(range_scale_x){
    g <- g %>%
      group_by(x) %>%
      mutate(value = range_scale(value, c(0,1))) %>%
      ungroup()
  }
  
  if(range_scale_y){
    g <- g %>%
      group_by(y) %>%
      mutate(value = range_scale(value, c(0,1))) %>%
      ungroup()
  }
  
  if(is.null(xaxis_order)){ 
    g <- g %>%
      mutate(x = factor(x, levels = labels(x_dendrogram))) 
  } else {
    if(xdendrogram) warning("displaying dendrogram and giving explicit axis order is not recommended! Axis and dendrogram will show conflicting order.")
    g <- g %>%
      mutate(x = factor(x, levels = xaxis_order)) 
  }
  
  if(is.null(yaxis_order)){ 
    g <- g %>%
      mutate(y = factor(y, levels = labels(y_dendrogram))) 
  } else {
    if(ydendrogram) warning("displaying dendrogram and giving explicit axis order is not recommended! Axis and dendrogram will show conflicting order.")
    g <- g %>%
      mutate(y = factor(y, levels = yaxis_order)) 
  }
  
    
  if(is.na(fill_factor[1])) {
    g <- g %>%
      ggplot() +
      geom_tile(aes(x, y, fill = value), show.legend = show.legend) + 
      geom_tile(data = filter(g, under_limit), aes(x, y), fill = under_lim_tile_color, show.legend = F)
  } else {
    g <- g %>%
      ggplot() +
      geom_tile(aes(x, y, fill = fill_factor), show.legend = show.legend) + 
      geom_tile(data = filter(g, under_limit), aes(x, y), fill = under_lim_tile_color, show.legend = F)
  }
  
  
  xp <- 0
  yp <- 0
  
  
  x_factors <- length(labels(x_dendrogram))
  y_factors <- length(labels(y_dendrogram))
  
  if(xdendrogram) {
    
    if(!is.null(xdendrogram_color)) {
      g <- g + 
        geom_segment(data = x_dendrogram_segments %>%
                       mutate(label = xdendrogram_color[match(x_dendrogram_segments$label, names(xdendrogram_color))]), 
                     aes(x = x, 
                         xend = xend,
                         y = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(y_factors, y_factors/(1 - y_margin)), 
                                                dodge = 0.5), 
                         yend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(y_factors, y_factors/(1 - y_margin)), 
                                                   dodge = 0.5),
                         color = label))
    } else{
      g <- g + 
        geom_segment(data = x_dendrogram_segments, 
                     aes(x = x, 
                         xend = xend,
                         y = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(y_factors, y_factors/(1 - y_margin)), 
                                                dodge = 0.5), 
                         yend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(y_factors, y_factors/(1 - y_margin)), 
                                                   dodge = 0.5)))
    }
    
    
    xp <- c(x_factors + x_margin)
  } 
  
  if(ydendrogram) {
    if(!is.null(ydendrogram_color)) {
      g <- g + 
        geom_segment(data = y_dendrogram_segments %>%
                       mutate(label = ydendrogram_color[match(y_dendrogram_segments$label, names(ydendrogram_color))]), 
                     aes(x = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(x_factors, x_factors/(1 - x_margin)), 
                                                dodge = 0.5), 
                         xend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(x_factors, x_factors/(1 - x_margin)), 
                                                   dodge = 0.5),
                         y = x, 
                         yend = xend,
                         color = label))
    } else{
      g <- g + 
        geom_segment(data = y_dendrogram_segments, 
                     aes(x = range_scale_manual(yend, max(y, yend), min(y, yend), 
                                                span = c(x_factors, x_factors/(1 - x_margin)), 
                                                dodge = 0.5), 
                         xend = range_scale_manual(y, max(y, yend), min(y, yend), 
                                                   span = c(x_factors, x_factors/(1 - x_margin)), 
                                                   dodge = 0.5),
                         y = x, 
                         yend = xend))
      
    }
    
    
    
    yp <- c(y_factors + y_margin)
  }
  
  g + 
    annotate("point", x = xp, y = yp, color = NA) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8), 
          panel.grid = element_blank())
  
  
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
    make_dendrogram(content_method, ensg_id, ex)
  
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




make_expression_heatmaps <- function(atlas.max.tb, atlas.cat, maxEx_column, tissue_column, 
                                     ensemblanno.table, proteinclass.table = NULL, proteinclass.table_ensg_id_column = "", 
                                     proteinclass.table_class_column = "", outpath, prefix, range_scale_x = F) {
  genes <- 
    atlas.max.tb %>%
    rename("content" = tissue_column, 
           "expression" = maxEx_column) %>%
    mutate(expression = log10(expression + 1)) %>% 
  
    left_join(select(atlas.cat, ensg_id, elevated.category, `enriched tissues`), by = "ensg_id") %>%
    left_join(ensemblanno.table, by = "ensg_id") %>% 
    group_by(gene_name, content) %>%
    mutate(gene_name_occurence = 1:length(unique(ensg_id))) %>%
    ungroup() %>%
    mutate(nonunique_gene_name = gene_name %in% {filter(., gene_name_occurence != 1) %$% 
                                                  unique(gene_name)},
           gene_name_2 = ifelse(nonunique_gene_name, 
                                paste(gene_name, gene_name_occurence),
                                gene_name),
           
           ###
           
           from = gene_name_2,
           to = content)
  
    
  if(!is.null(proteinclass.table)) {
    genes <- 
      genes %>%
      left_join(select(proteinclass.table, 
                       ensg_id = proteinclass.table_ensg_id_column, 
                       protein_class = proteinclass.table_class_column), by = "ensg_id")
  }
  
  
  genes_median_expression <- 
    genes %>%
    {mapply(unique(.$content), 
            FUN = function(cell1) mapply(unique(.$content), 
                                         FUN = function(cell2) {
                                           genes <- 
                                             filter(., content == cell1)$ensg_id
                                           atlas.max.tb %>%
                                             filter(eval(parse(text = tissue_column)) == cell2 & 
                                                      ensg_id %in% genes) %$%
                                             median(eval(parse(text = maxEx_column)), na.rm = T)
                                         }))} %>%
    as_tibble(., rownames = "from") %>%
    gather(key = "to", value = "median expression", -from)
  #-----------------------------------------------------------------
  # Tissue enriched genes
  plot.data <- 
    genes %>%
    filter(elevated.category == "tissue enriched") 
  
  ## 1
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1))+
    scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
    theme(axis.text.x = element_blank(), axis.title.y = element_blank()) + 
    ggtitle(paste0("Tissue enriched genes (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__tissue_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  
  ## 2
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1),
                 xdendrogram = F, 
                 xaxis_order = unique(from[order(`enriched tissues`)]))+
    scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
    theme(axis.text.x = element_blank(), axis.title.y = element_blank()) + 
    ggtitle(paste0("Tissue enriched genes (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__tissue_enriched_expression_dendroheatmap2.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  ## Group enriched genes
  plot.data <- 
    genes %>%
    filter(elevated.category == "group enriched")
  ## 1
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1))+
    scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 1), axis.title.y = element_blank()) + 
    ggtitle(paste0("Group enriched genes (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__group_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  ## All elevated genes
  plot.data <- 
    genes %>%
    filter(elevated.category %in% c("tissue enriched", "group enriched", "tissue enhanced"))
  ## 1
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1))+
    scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
    theme(axis.text.x = element_blank(), axis.title.y = element_blank()) + 
    ggtitle(paste0("All elevated genes (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_elevated_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  
  ## categories
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1), 
                 fill_factor = elevated.category)+
    scale_fill_manual(values = elevated.cat.cols) + 
    theme(axis.text.x = element_blank(), axis.title.y = element_blank()) + 
    ggtitle(paste0("All elevated genes (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_elevated_category_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  ## potential CD marker genes
  if(!is.null(proteinclass.table)) {
    plot.data <- 
      genes %>%
      filter(elevated.category %in% c("tissue enriched", "group enriched", "tissue enhanced") &
               protein_class == "membrane")
    ## All elevated
    plot.data %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), axis.title.y = element_blank()) + 
      ggtitle(paste0("Potential CD markers: All elevated genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_pot_CD_all_elevated_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
    
    ## Group enriched
    plot.data %>%
      filter(elevated.category == "group enriched") %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.title.y = element_blank()) + 
      ggtitle(paste0("Potential CD markers: Group enriched genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_pot_CD_group_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
    
    ## Tissue enriched
    plot.data %>%
      filter(elevated.category == "tissue enriched") %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.title.y = element_blank()) + 
      ggtitle(paste0("Potential CD markers: Tissue enriched genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_pot_CD_tissue_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  }
  #-----------------------------------------------------------------
  ## Known CD marker genes
  if(!is.null(proteinclass.table)) {
    plot.data <- 
      genes %>%
      filter(protein_class == "cd_marker")
    ## All elevated
    plot.data %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), axis.title.y = element_blank()) + 
      ggtitle(paste0("Known CD markers: All elevated genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_known_CD_all_elevated_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
    
    ## Group enriched
    plot.data %>%
      filter(elevated.category == "group enriched") %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.title.y = element_blank()) + 
      ggtitle(paste0("Known CD markers: Group enriched genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_known_CD_group_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
    
    ## Tissue enriched
    plot.data %>%
      filter(elevated.category == "tissue enriched") %$%
      {ggdendroheat(from, to, expression, 
                   x_clustering_method = "euclidean", 
                   range_scale_x = range_scale_x, 
                   under_lim_tile_color = "white", 
                   under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), axis.title.y = element_blank()) + 
      ggtitle(paste0("Known CD markers: Tissue enriched genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
    ggsave(paste(outpath, paste0(prefix, '_known_CD_tissue_enriched_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  }
  
  ## categories
  plot.data %$%
    {ggdendroheat(from, to, expression, 
                 x_clustering_method = "euclidean", 
                 range_scale_x = range_scale_x, 
                 under_lim_tile_color = "white", 
                 under_lim = log10(1 + 1), 
                 fill_factor = elevated.category)+
    scale_fill_manual(values = elevated.cat.cols) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4), axis.title.y = element_blank()) + 
    ggtitle(paste0("Known CD markers: Categories (n = ", length(unique(from)), ")")) + 
    xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '_known_CD_category_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  
  
}

make_immunodeficiency_expression_heatmaps <- function(atlas.max.tb, atlas.cat, immunodeficiency.table, maxEx_column, tissue_column, 
                                     ensemblanno.table, proteinclass.table = NULL, proteinclass.table_ensg_id_column = "", 
                                     proteinclass.table_class_column = "", outpath, prefix, range_scale_x = F) {
  genes <- 
    atlas.max.tb %>%
    filter(ensg_id %in% immunodeficiency.table$ensg_id) %>%
    left_join(immunodeficiency.table %>%
                select(1, 2, ensg_id) %>%
                mutate(Major_ID = gsub("\\.", "", str_extract( `Major groups of PIDs`, ".{1,2}\\."))) %>% 
                group_by(ensg_id) %>% 
                summarise(`Major groups of PIDs` = paste(unique(`Major groups of PIDs`), collapse = ", "),
                          `Subgroups of PIDs` = paste(unique(`Subgroups of PIDs`), collapse = ", "),
                          Major_ID = paste(unique(Major_ID), collapse = ", ")), by = "ensg_id") %>%
    rename("content" = tissue_column, 
           "expression" = maxEx_column) %>%
    mutate(expression = log10(expression + 1)) %>% 
    
    left_join(select(atlas.cat, ensg_id, elevated.category, `enriched tissues`), by = "ensg_id") %>%
    left_join(ensemblanno.table, by = "ensg_id") %>% 
    group_by(gene_name, content) %>%
    mutate(gene_name_occurence = 1:length(unique(ensg_id))) %>%
    ungroup() %>%
    mutate(nonunique_gene_name = gene_name %in% {filter(., gene_name_occurence != 1) %$% 
        unique(gene_name)},
        gene_name_2 = ifelse(nonunique_gene_name, 
                             paste(gene_name, gene_name_occurence),
                             gene_name),
        
        ###
        
        from = gene_name_2,
        to = content)
  
  
  if(!is.null(proteinclass.table)) {
    genes <- 
      genes %>%
      left_join(select(proteinclass.table, 
                       ensg_id = proteinclass.table_ensg_id_column, 
                       protein_class = proteinclass.table_class_column), by = "ensg_id")
  }
  
  
  
  
  #-----------------------------------------------------------------
  ## All elevated genes
  plot.data <- 
    genes %>%
    filter(elevated.category %in% c("tissue enriched", "group enriched", "tissue enhanced"))
  ## 1
  plot.data %$%
  {ggdendroheat(from, to, expression, 
                x_clustering_method = "euclidean", 
                range_scale_x = range_scale_x, 
                under_lim_tile_color = "white", 
                under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(size = 5), axis.title.y = element_blank()) + 
      ggtitle(paste0("All elevated genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_elevated_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  
  ## categories
  plot.data %$%
  {ggdendroheat(from, to, expression, 
                x_clustering_method = "euclidean", 
                range_scale_x = range_scale_x, 
                under_lim_tile_color = "white", 
                under_lim = log10(1 + 1), 
                fill_factor = elevated.category)+
      scale_fill_manual(values = elevated.cat.cols) + 
      theme(axis.text.x = element_text(size = 5), axis.title.y = element_blank()) + 
      ggtitle(paste0("All elevated genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_elevated_category_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  #-----------------------------------------------------------------
  ## All genes
  plot.data <- 
    genes 
  ## 1
  plot.data %$%
  {ggdendroheat(from, to, expression, 
                x_clustering_method = "euclidean", 
                range_scale_x = range_scale_x, 
                under_lim_tile_color = "white", 
                under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = c("white", "yellow", "orangered", "#800026")) + 
      theme(axis.text.x = element_text(size = 5), axis.title.y = element_blank()) + 
      ggtitle(paste0("All PID genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_genes_expression_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  
  plot.data %$%
  {ggdendroheat(from, to, expression, 
                x_clustering_method = "euclidean", 
                range_scale_x = range_scale_x, 
                under_lim_tile_color = "black", 
                under_lim = log10(1 + 1))+
      scale_fill_gradientn(colors = viridis::magma(4)) + 
      theme(axis.text.x = element_text(size = 5), axis.title.y = element_blank()) + 
      ggtitle(paste0("All PID genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_genes_expression_dendroheatmap_magma.pdf'),sep='/'), width=15, height=10)
  
  ## categories
  plot.data %$%
  {ggdendroheat(from, to, expression, 
                x_clustering_method = "euclidean", 
                range_scale_x = range_scale_x, 
                under_lim_tile_color = "white", 
                under_lim = log10(1 + 1), 
                fill_factor = elevated.category)+
      scale_fill_manual(values = elevated.cat.cols) + 
      theme(axis.text.x = element_text(size = 5), axis.title.y = element_blank()) + 
      ggtitle(paste0("All elevated genes (n = ", length(unique(from)), ")")) + 
      xlab("genes")}
  ggsave(paste(outpath, paste0(prefix, '__all_genes_category_dendroheatmap.pdf'),sep='/'), width=15, height=10)
  #-----------------------------------------------------------------
  
}

make_score_expression_scatter <- function(atlas.max.tb, atlas.cat, maxEx_column, tissue_column, ensemblanno.table, plot.order = NULL, outpath, prefix) {
  plot.data <- 
    atlas.max.tb %>%
    ungroup() %>%
    rename(tissue_column = tissue_column,
           maxEx_column = maxEx_column) %>%
    left_join(atlas.cat, by = "ensg_id") %>%
    left_join(select(ensemblanno.table, ensg_id, gene_name) , by = "ensg_id") %>%
    
    # Filter so that elevated genes are only kept for the tissues they are elevated in
    filter(mapply(paste0("(^|, )", tissue_column, "(, |$)"), 
                  `enriched tissues`, FUN = function(x,y) grepl(x, y))) %>%
    mutate(score = as.numeric(`tissue/group specific score`)) %>%
    filter(!is.na(score)) 
  
  if(!is.null(plot.order)){
    plot.data <- 
      plot.data %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.order))
  }
  
  plot.data %>%
    ggplot(aes(maxEx_column, score, color = elevated.category)) +
    geom_point() +
    facet_wrap(tissue_column ~ .) + 
    xlab("expression") +
    ylab("score") + 
    scale_x_log10() + 
    scale_y_log10() + 
    scale_color_manual(values = elevated.cat.cols) +
    simple_theme
  
  ggsave(paste(outpath, paste0(prefix, '_score_expression_scatter.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    group_by(tissue_column) %>%
    mutate(label_size = range_scale(range_scale(score, 1:2) * 
                                      range_scale(maxEx_column, 1:2), 
                                    span = c(1,4))) %>%
    ggplot(aes(maxEx_column, score, color = elevated.category, label = gene_name, size = label_size)) +
    geom_text() +
    facet_wrap(tissue_column ~ .) + 
    xlab("expression") +
    ylab("score") + 
    scale_x_log10() + 
    scale_y_log10() + 
    scale_color_manual(values = elevated.cat.cols) +
    simple_theme
  
  ggsave(paste(outpath, paste0(prefix, '_score_expression_labels.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    ggplot(aes(maxEx_column, fill = elevated.category)) +
    geom_density(alpha = 0.5) +
    facet_wrap(tissue_column ~ .) + 
    xlab("expression") +
    scale_x_log10() + 
    scale_fill_manual(values = elevated.cat.cols) +
    simple_theme
  ggsave(paste(outpath, paste0(prefix, '_expression_density.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    ggplot(aes(score, fill = elevated.category)) +
    geom_density(alpha = 0.5) +
    facet_wrap(tissue_column ~ .) + 
    xlab("score") +
    scale_x_log10() + 
    scale_fill_manual(values = elevated.cat.cols) +
    simple_theme
  ggsave(paste(outpath, paste0(prefix, '_score_density.pdf'),sep='/'), width=15, height=10)
  

}

elevated_NX_fraction_barplots <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, outpath, prefix) {
  atlas.max_temp <- 
    atlas.max %>% 
    rename(maxEx_column = maxEx_column,
           tissue_column = tissue_column)
  
  calculate_elevated_sum_NX <- function(tissue) {
    atlas.cat %>%
      filter(grepl(paste0("(^|, )", unique(tissue), "(, |$)"), `enriched tissues`)) %>%
      {filter(atlas.max_temp, ensg_id %in% .$ensg_id & tissue_column == unique(tissue))} %$%
      sum(maxEx_column)
  }
  
  plot.data <- 
    atlas.max_temp %>%
    group_by(tissue_column) %>% 
    summarise(total_NX = sum(maxEx_column),
              elevated_NX = calculate_elevated_sum_NX(tissue_column), 
              fraction = elevated_NX / total_NX) 
  
  plot.data %>%
    mutate(., tissue_column = factor(tissue_column, levels = tissue_column[order(total_NX)])) %>%
    ggplot(aes(tissue_column, total_NX))+
    geom_bar(stat = "identity", fill = "white", color = "darkgray") + 
    geom_bar(aes(tissue_column, elevated_NX), stat = "identity", fill = "red") + 
    geom_text(aes(tissue_column, elevated_NX, label = paste0(round(fraction*100, digits = 1), "%")), hjust = -0.1)+
    simple_theme + 
    coord_flip()
  
  ggsave(paste(outpath, paste0(prefix, '_total_NX_barplot1.pdf'),sep='/'), width=15, height=10)
  
  plot.data %>%
    mutate(., tissue_column = factor(tissue_column, levels = tissue_column[order(elevated_NX)])) %>%
    ggplot(aes(tissue_column, total_NX))+
    geom_bar(stat = "identity", fill = "white", color = "darkgray") + 
    geom_bar(aes(tissue_column, elevated_NX), stat = "identity", fill = "red") + 
    geom_text(aes(tissue_column, elevated_NX, label = paste0(round(fraction*100, digits = 1), "%")), hjust = -0.1)+
    simple_theme + 
    coord_flip()
  
  ggsave(paste(outpath, paste0(prefix, '_total_NX_barplot2.pdf'),sep='/'), width=15, height=10)

}
make_elevated_NX_fraction_barplots <- function(atlas.max, atlas.cat, maxEx_column, tissue_column, outpath, prefix) {
  elevated_NX_fraction_barplots(atlas.max, atlas.cat, maxEx_column, tissue_column, outpath, prefix)
  
  #### intrecellular proteins
  proteinlocalization.table %>% 
    filter(predicted_intracellular) %>%
    left_join(atlas.max, by = "ensg_id") %>%
    elevated_NX_fraction_barplots(atlas.cat, maxEx_column, tissue_column, outpath, prefix = paste0(prefix, "_intracellular"))
  
  #### secreted proteins
  proteinlocalization.table %>% 
    filter(predicted_secreted) %>%
    left_join(atlas.max, by = "ensg_id") %>%
    elevated_NX_fraction_barplots(atlas.cat, maxEx_column, tissue_column, outpath, prefix = paste0(prefix, "_secreted"))
  
  #### membrane proteins
  proteinlocalization.table %>% 
    filter(predicted_membrane) %>%
    left_join(atlas.max, by = "ensg_id") %>%
    elevated_NX_fraction_barplots(atlas.cat, maxEx_column, tissue_column, outpath, prefix = paste0(prefix, "_membrane"))
  
  #######
    
  atlas.max_temp <- 
    atlas.max %>% 
    rename(maxEx_column = maxEx_column,
           tissue_column = tissue_column)
  
  calculate_elevated_sum_NX <- function(tissue) {
    atlas.cat %>%
      filter(grepl(paste0("(^|, )", unique(tissue), "(, |$)"), `enriched tissues`)) %>%
      {filter(atlas.max_temp, ensg_id %in% .$ensg_id & tissue_column == unique(tissue))} %$%
      sum(maxEx_column)
  }
  
  plot.data <- 
    atlas.max_temp %>%
    left_join(proteinlocalization.table %>%
                mutate(location = case_when(!predicted_localization_class %in% c("membrane", "secreted", "intracellular") ~ "multiple locations",
                                            predicted_secreted ~ "secreted",
                                            predicted_membrane ~ "membrane",
                                            predicted_intracellular ~ "intracellular",
                                            T ~ "")),
              by = "ensg_id") %>%
    group_by(tissue_column, location) %>% 
    summarise(total_NX = sum(maxEx_column)) %>%
    left_join(group_by(., tissue_column) %>% summarise(sum_total_NX = sum(total_NX)), by = "tissue_column") %>%
    mutate(fraction = total_NX / sum_total_NX )
  
  write_delim(plot.data, paste(outpath, paste0(prefix, '_sum_abundance_location.txt'),sep='/'), delim = "\t")
  
  
  plot.data_order <- 
    plot.data %>%
    group_by(tissue_column) %>%
    summarise(total_NX = sum(total_NX)) %$% 
    tissue_column[order(total_NX)]
  
  plot.data %>%
    ungroup() %>%
    mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
           location = factor(location, levels = c("multiple locations", "intracellular", "membrane", "secreted"))) %>%
    ggplot(aes(tissue_column, total_NX, fill = location))+
    geom_bar(stat = "identity", color = "darkgray") + 
    simple_theme + 
    coord_flip() +
    scale_fill_manual(values = protein.localization.palette)
  
  ggsave(paste(outpath, paste0(prefix, '_total_NX_location.pdf'),sep='/'), width=15, height=10)
  
  # Classification elevated
  plot.data <- 
    atlas.max_temp %>%
    left_join(atlas.cat,
              by = "ensg_id") %>%
    group_by(tissue_column, elevated.category) %>% 
    summarise(total_NX = sum(maxEx_column)) %>%
    left_join(group_by(., tissue_column) %>% summarise(sum_total_NX = sum(total_NX)), by = "tissue_column") %>%
    mutate(fraction = total_NX / sum_total_NX )
  
  write_delim(plot.data, paste(outpath, paste0(prefix, '_sum_abundance_elevated_category.txt'),sep='/'), delim = "\t")
  
  plot.data_order <- 
    plot.data %>%
    group_by(tissue_column) %>%
    summarise(total_NX = sum(total_NX)) %$% 
    tissue_column[order(total_NX)]
  
  plot.data %>%
    ungroup() %>%
    mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
           elevated.category = factor(elevated.category, levels = c("not detected", "low tissue specificity", 
                                                                    "tissue enhanced", "group enriched", "tissue enriched"))) %>%
    ggplot(aes(tissue_column, total_NX, fill = elevated.category))+
    geom_bar(stat = "identity", color = "darkgray") + 
    simple_theme + 
    coord_flip() +
    scale_fill_manual(values = elevated.cat.cols)
  
  ggsave(paste(outpath, paste0(prefix, '_total_NX_elevation.pdf'),sep='/'), width=15, height=10)
  
  # Classification distribution
  plot.data <- 
    atlas.max_temp %>%
    left_join(atlas.cat,
              by = "ensg_id") %>%
    group_by(tissue_column, express.category.2) %>% 
    summarise(total_NX = sum(maxEx_column)) %>%
    left_join(group_by(., tissue_column) %>% summarise(sum_total_NX = sum(total_NX)), by = "tissue_column") %>%
    mutate(fraction = total_NX / sum_total_NX )
  
  write_delim(plot.data, paste(outpath, paste0(prefix, '_sum_abundance_distribution_category.txt'),sep='/'), delim = "\t")
  
  plot.data_order <- 
    plot.data %>%
    group_by(tissue_column) %>%
    summarise(total_NX = sum(total_NX)) %$% 
    tissue_column[order(total_NX)]
  
  plot.data %>%
    ungroup() %>%
    mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
           express.category.2 = factor(express.category.2, levels = c("not expressed", "expressed in all", "expressed in many",
                                                                      "expressed in some", "expressed in single"))) %>%
    ggplot(aes(tissue_column, total_NX, fill = express.category.2))+
    geom_bar(stat = "identity", color = "darkgray") + 
    simple_theme + 
    coord_flip() +
    scale_fill_manual(values = expressed.cat.cols)
  
  ggsave(paste(outpath, paste0(prefix, '_total_NX_distribution.pdf'),sep='/'), width=15, height=10)
  
}

make_sum_TPM_plot <- function(atlas, atlas.cat, tissue_column, method_column, outpath, prefix) {
  atlas_temp <- 
    atlas %>%
    rename(method_column = method_column, 
           tissue_column = tissue_column) %>%
    left_join(group_by(., method_column, tissue) %>%
                summarise(total_tpm = sum(expression, na.rm = T)),
              by = c("method_column", "tissue")) %>%
    mutate(expression = expression/(total_tpm / 1e6))
  
  for(method in unique(atlas_temp$method_column)) {
    plot.data <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(atlas.cat, by = "ensg_id") %>%
      group_by(tissue_column, elevated.category) %>%
      summarise(sum_TPM = sum(expression, na.rm = T))
    
    plot.data_order <- 
      plot.data %>%
      filter(elevated.category == "tissue enriched") %>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
      
    plot.data %>%
      ungroup() %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
             elevated.category = factor(elevated.category, levels = c("not detected", "low tissue specificity", 
                                                                      "tissue enhanced", "group enriched", "tissue enriched"))) %>%
      ggplot(aes(tissue_column, sum_TPM, fill = elevated.category))+
      geom_bar(stat = "identity", color = "darkgray") + 
      simple_theme + 
      coord_flip() +
      scale_fill_manual(values = elevated.cat.cols) + 
      ggtitle(method)
      
    ggsave(paste(outpath, paste0(prefix, "_", method, '_total_TPM_elevation.pdf'),sep='/'), width=15, height=10)
    
    plot.data <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(atlas.cat, by = "ensg_id") %>%
      group_by(tissue_column, express.category.2) %>%
      summarise(sum_TPM = sum(expression, na.rm = T))
    
    plot.data_order <- 
      plot.data %>%
      filter(express.category.2 == "expressed in single") %>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
    
    plot.data %>%
      ungroup() %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
             express.category.2 = factor(express.category.2, levels = c("not expressed", "expressed in all", "expressed in many",
                                                                        "expressed in some", "expressed in single"))) %>%
      ggplot(aes(tissue_column, sum_TPM, fill = express.category.2))+
      geom_bar(stat = "identity", color = "darkgray") + 
      simple_theme + 
      coord_flip() +
      scale_fill_manual(values = expressed.cat.cols)+ 
      ggtitle(method)
    
    ggsave(paste(outpath, paste0(prefix, "_", method, '_total_TPM_distribution.pdf'),sep='/'), width=15, height=10)
    
    
    # Localization 
    plot.data <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(proteinlocalization.table %>%
                  mutate(location = case_when(!predicted_localization_class %in% c("membrane", "secreted", "intracellular") ~ "multiple locations",
                                              predicted_secreted ~ "secreted",
                                              predicted_membrane ~ "membrane",
                                              predicted_intracellular ~ "intracellular",
                                              T ~ "")),
                by = "ensg_id") %>%
      group_by(tissue_column, location) %>%
      summarise(sum_TPM = sum(expression, na.rm = T))
    
    plot.data_order <- 
      plot.data %>%
      filter(location == "secreted")%>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
    
    plot.data %>%
      ungroup() %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.data_order),
             location = factor(location, levels = c("multiple locations", "intracellular", "membrane", "secreted"))) %>%
      ggplot(aes(tissue_column, sum_TPM, fill = location))+
      geom_bar(stat = "identity", color = "darkgray") + 
      simple_theme + 
      coord_flip() +
      scale_fill_manual(values = protein.localization.palette)+ 
      ggtitle(method)
    
    ggsave(paste(outpath, paste0(prefix, "_", method, '_total_TPM_location.pdf'),sep='/'), width=15, height=10)
    
    
    # Localization 2
    plot.data <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(proteinlocalization.table,
                by = "ensg_id") %>%
      mutate(location = predicted_localization_class) %>%
      group_by(tissue_column, location) %>%
      summarise(sum_TPM = sum(expression, na.rm = T))
    
    plot.data_order <- 
      plot.data %>%
      filter(location == "secreted")%>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
    
    plot.data %>%
      ungroup() %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.data_order)) %>%
      ggplot(aes(tissue_column, sum_TPM, fill = location))+
      geom_bar(stat = "identity", color = "darkgray") + 
      simple_theme + 
      coord_flip() +
      scale_fill_manual(values = protein.localization.palette)+ 
      ggtitle(method)
    
    ggsave(paste(outpath, paste0(prefix, "_", method, '_total_TPM_location_all.pdf'),sep='/'), width=15, height=10)
    
    
    # Localization hierarchy
    plot.data <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(proteinlocalization.table %>%
                  mutate(location = case_when(predicted_secreted ~ "secreted",
                                              predicted_membrane ~ "membrane",
                                              predicted_intracellular ~ "intracellular",
                                              T ~ "")),
                by = "ensg_id") %>%
      group_by(tissue_column, location) %>%
      summarise(sum_TPM = sum(expression, na.rm = T))
    
    plot.data_order <- 
      plot.data %>%
      filter(location == "secreted")%>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
    
    plot.data %>%
      ungroup() %>%
      mutate(tissue_column = factor(tissue_column, levels = plot.data_order)) %>%
      ggplot(aes(tissue_column, sum_TPM, fill = location))+
      geom_bar(stat = "identity", color = "darkgray") + 
      simple_theme + 
      coord_flip() +
      scale_fill_manual(values = protein.localization.palette)+ 
      ggtitle(method)
    
    ggsave(paste(outpath, paste0(prefix, "_", method, '_total_TPM_location_hierarchy.pdf'),sep='/'), width=15, height=10)
    
    
  }
    
}
####

make_sum_TPM_gene_bar_plot <- function(atlas, atlas.cat, tissue_column, method_column, outpath, prefix) {
  atlas_temp <- 
    atlas %>%
    rename(method_column = method_column, 
           tissue_column = tissue_column) %>%
    left_join(group_by(., method_column, tissue) %>%
                summarise(total_tpm = sum(expression, na.rm = T)),
              by = c("method_column", "tissue")) %>%
    mutate(expression = expression/(total_tpm / 1e6))
  
  for(method in unique(atlas_temp$method_column)) {
    plot.data_order <- 
      atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(proteinlocalization.table %>%
                  mutate(location = case_when(predicted_secreted ~ "secreted",
                                              predicted_membrane ~ "membrane",
                                              predicted_intracellular ~ "intracellular",
                                              T ~ "")),
                by = "ensg_id") %>%
      group_by(tissue_column, location) %>%
      summarise(sum_TPM = sum(expression, na.rm = T)) %>%
      filter(location == "secreted") %>%
      group_by(tissue_column) %>%
      summarise(sum_TPM = sum(sum_TPM)) %$% 
      tissue_column[order(sum_TPM)]
    
    atlas_temp %>%
      filter(method_column == method) %>% 
      left_join(proteinlocalization.table %>%
                  mutate(location = case_when(predicted_secreted ~ "secreted",
                                              predicted_membrane ~ "membrane",
                                              predicted_intracellular ~ "intracellular",
                                              T ~ "")),
                by = "ensg_id") %>%
                {.[order(.$location, .$expression),]} %>% 
      group_by(tissue_column) %>%
      mutate(gene_order = 1:length(unique(ensg_id))) %>%
      group_by(tissue_column, location) %>%
      mutate(cum_expression = cumsum(expression)) %>%
      #filter(location == "intracellular") %>%
      {ggplot(., aes(gene_order, cum_expression, fill = location, color = location))+
          #geom_line()+
          #geom_ribbon()+
          geom_area(alpha = 0.3, position = "identity")+
          facet_wrap(. ~ tissue_column)+
          simple_theme + 
          scale_fill_manual(values = protein.localization.palette)+ 
          scale_color_manual(values = protein.localization.palette)+ 
          ggtitle(method)}
    ggsave(paste(outpath, paste0(prefix, "_", method, '_TPM_location_cumsum_flagplot.pdf'),sep='/'), width=15, height=10)
    
    ####
    atlas_temp %>%
      filter(method_column == method) %>% 
      
      #filter(tissue_column == "pancreas") %>%
      
      left_join(proteinlocalization.table %>%
                  mutate(location = case_when(predicted_secreted ~ "secreted",
                                              predicted_membrane ~ "membrane",
                                              predicted_intracellular ~ "intracellular",
                                              T ~ "")),
                by = "ensg_id") %>%
      mutate(display = expression > 5000,
             display_label = expression > 10000,
             label = ifelse(display_label, display_name, NA)) %>%
             {rbind(.[.$display,] %>%
                      select(tissue_column, location, expression, label) %>%
                      mutate(stacked = "1"),
                    .[!.$display,] %>%
                      group_by(tissue_column, location) %>%
                      summarise(expression = sum(expression, na.rm = T),
                                label = paste0("n = ", length(ensg_id))) %>%
                      ungroup() %>%
                      mutate(stacked = "2"))} %>%
      mutate(label_size = ifelse(expression > 50000, 50000, expression)) %>%
                      {.[order(.$stacked, .$expression),]} %>%
      mutate(label = factor(label, levels = unique(label[order(location, stacked)])), 
             tissue_column = factor(tissue_column, levels = rev(plot.data_order))) %>%
             {ggplot(., aes(tissue_column, expression, fill = paste(location, stacked), color = paste(location, stacked)))+
                 geom_bar(stat = "identity", color = "white")+
                 geom_text(aes(label = label, size = label_size),
                           position = position_stack(vjust = 0.5), color = "white", fontface = "bold")+
                 simple_theme + 
                 
                 scale_size_continuous(range = c(0.25, 4))+
                 scale_fill_manual(values = c(setNames(protein.localization.palette,
                                                       # sapply(protein.localization.palette, 
                                                       #        FUN = function(x) colorRampPalette(c(x, "white"))(7)[2]),
                                                       paste(names(protein.localization.palette),
                                                             "2")), 
                                              setNames(sapply(protein.localization.palette, 
                                                              FUN = function(x) colorRampPalette(c(x, "white"))(6)[2]),
                                                       paste(names(protein.localization.palette),
                                                             "1"))))+ 
                 # scale_color_manual(values = c(setNames(protein.localization.palette,
                 #                                        paste(names(protein.localization.palette),
                 #                                              "1")), 
                 #                               setNames(sapply(protein.localization.palette, 
                 #                                               FUN = function(x) colorRampPalette(c(x, "white"))(3)[2]),
                 #                                        paste(names(protein.localization.palette),
                 #                                              "2"))))+ 
                 
                 ggtitle(method)}
    ggsave(paste(outpath, paste0(prefix, "_", method, '_TPM_bar_abundant_genes.pdf'),sep='/'), width=40, height=10)
  }
  
}


make_classification_pie_chart <- function(atlas.cat, outpath, prefix) {
  atlas.cat %>%
    group_by(elevated.category) %>%
    summarise(number = length(ensg_id)) %>%
    {.[match(c("tissue enriched",
               "group enriched",
               "tissue enhanced",
               "low tissue specificity",
               "not detected"),
             .$elevated.category), ]} %>%
    mutate(elevated.category = factor(elevated.category, levels = rev(c("tissue enriched",
                                                                        "group enriched",
                                                                        "tissue enhanced",
                                                                        "low tissue specificity",
                                                                        "not detected"))),
           label_y = cumsum(number) - number / 2) %>%
    ggplot(aes("", number, fill = elevated.category)) +
    geom_bar(stat = "identity", show.legend = T, color = "white", size = 1)+
    #geom_text(aes(x = 1.8, y = label_y, label = elevated.category))+
    geom_text(aes(x = 1.3, 
                  y = label_y, 
                  label = paste0(number, "\n", 
                                 "(", round(100 * number / sum(number), digits = 1), "%)")), 
              color = "black",
              size = 4)+
    coord_polar("y")+
    scale_fill_manual(values = elevated.cat.cols) + 
    theme_void() + 
    theme(legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.5, "cm"),
          legend.title = element_text(face="italic"))
  ggsave(paste(outpath, paste0(prefix, '_elevated_pie1.pdf'),sep='/'), width=8, height=8)
  
  atlas.cat %>%
    group_by(elevated.category) %>%
    summarise(number = length(ensg_id)) %>%
    {.[match(c("tissue enriched",
               "group enriched",
               "tissue enhanced",
               "low tissue specificity",
               "not detected"),
             .$elevated.category), ]} %>%
    mutate(elevated.category = factor(elevated.category, levels = rev(c("tissue enriched",
                                                                        "group enriched",
                                                                        "tissue enhanced",
                                                                        "low tissue specificity",
                                                                        "not detected"))),
           label_y = cumsum(number) - number / 2) %>%
    ggplot(aes("", number, fill = elevated.category)) +
    geom_bar(stat = "identity", show.legend = F, color = "white", size = 1)+
    geom_text(aes(x = 1.8, 
                  y = label_y, 
                  label = elevated.category),
              size = 4)+
    geom_text(aes(x = 1.3, 
                  y = label_y, 
                  label = paste0(number, "\n", 
                                 "(", round(100 * number / sum(number), digits = 1), "%)")), 
              color = "black",
              size = 4)+
    coord_polar("y")+
    scale_fill_manual(values = elevated.cat.cols) + 
    theme_void() 
  ggsave(paste(outpath, paste0(prefix, '_elevated_pie2.pdf'),sep='/'), width=8, height=8)
  
  
  ### Distribution
  
  atlas.cat %>%
    group_by(express.category.2) %>%
    summarise(number = length(ensg_id)) %>%
    {.[match(c("expressed in single",
               "expressed in some",
               "expressed in many",
               "expressed in all",
               "not expressed"),
             .$express.category.2), ]} %>%
    filter(!is.na(number)) %>%
    mutate(express.category.2 = factor(express.category.2, levels = rev(c("expressed in single",
                                                                        "expressed in some",
                                                                        "expressed in many",
                                                                        "expressed in all",
                                                                        "not expressed"))),
           label_y = cumsum(number) - number / 2) %>%
    ggplot(aes("", number, fill = express.category.2)) +
    geom_bar(stat = "identity", show.legend = T, color = "white", size = 1)+
    #geom_text(aes(x = 1.8, y = label_y, label = elevated.category))+
    geom_text(aes(x = 1.3, 
                  y = label_y, 
                  label = paste0(number, "\n", 
                                 "(", round(100 * number / sum(number), digits = 1), "%)")), 
              color = "black",
              size = 4)+
    coord_polar("y")+
    scale_fill_manual(values = expressed.cat.cols) + 
    theme_void() + 
    theme(legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.size= unit(0.5, "cm"),
          legend.title = element_text(face="italic"))
  ggsave(paste(outpath, paste0(prefix, '_distribution_pie1.pdf'),sep='/'), width=8, height=8)
  
  atlas.cat %>%
    group_by(express.category.2) %>%
    summarise(number = length(ensg_id)) %>%
    {.[match(c("expressed in single",
               "expressed in some",
               "expressed in many",
               "expressed in all",
               "not expressed"),
             .$express.category.2), ]} %>%
    filter(!is.na(number)) %>%
    mutate(express.category.2 = factor(express.category.2, levels = rev(c("expressed in single",
                                                                          "expressed in some",
                                                                          "expressed in many",
                                                                          "expressed in all",
                                                                          "not expressed"))),
           label_y = cumsum(number) - number / 2) %>%
    ggplot(aes("", number, fill = express.category.2)) +
    geom_bar(stat = "identity", show.legend = F, color = "white", size = 1)+
    geom_text(aes(x = 1.8, 
                  y = label_y, 
                  label = express.category.2),
              size = 4)+
    geom_text(aes(x = 1.3, 
                  y = label_y, 
                  label = paste0(number, "\n", 
                                 "(", round(100 * number / sum(number), digits = 1), "%)")), 
              color = "black",
              size = 4)+
    coord_polar("y")+
    scale_fill_manual(values = expressed.cat.cols) + 
    theme_void() 
  ggsave(paste(outpath, paste0(prefix, '_distribution_pie2.pdf'),sep='/'), width=8, height=8)
}




  



make_double_donut_chord <- function(cat1, cat2, regional_class_dictionary, grid.col, global_name, prefix, outpath) {
  
  elevated_name <- paste(global_name, "elevated")
  
  global_local_elevation <- 
    left_join(cat1 %>%
                select(ensg_id, `enriched tissues`, elevated.category), 
              cat2 %>%
                select(ensg_id, `enriched tissues`, elevated.category),
              by = "ensg_id") %>%
    filter(elevated.category.x != "not detected") %>%
    mutate(global_elevation = case_when(elevated.category.x == "low tissue specificity" ~ "low tissue specificity",
                                        grepl(global_name, `enriched tissues.x`) ~ elevated_name,
                                        elevated.category.x %in% c("tissue enhanced",
                                                                   "tissue enriched",
                                                                   "group enriched") ~ "elevated in other tissue"),
           local_elevation = case_when(elevated.category.y == "low tissue specificity" ~ regional_class_dictionary[1],
                                       elevated.category.y == "not detected" ~ regional_class_dictionary[2],
                                       elevated.category.y %in% c("tissue enhanced",
                                                                  "tissue enriched",
                                                                  "group enriched") ~ regional_class_dictionary[3])) 
  
  global_local_elevation_summary <- 
    global_local_elevation %>%
    group_by(global_elevation, local_elevation) %>%
    summarise(number = length(ensg_id))
  
  pdf(file = paste(outpath, paste0(prefix, '_double_donut_chord.pdf'),sep='/'), width=10, height=10, useDingbats = F)
  global_local_elevation_summary %$%
    chord_classification(from = global_elevation, 
                         to = local_elevation, 
                         sizes = number, 
                         groups = c(1, 1, 1, 2, 2, 2), 
                         grid.col = grid.col,
                         plot.order = c(elevated_name, "elevated in other tissue", "low tissue specificity", 
                                        regional_class_dictionary[2], regional_class_dictionary[1], regional_class_dictionary[3]), 
                         size_labels = T)
  dev.off()
  
}



make_FACS_boxplot <- function(sample_atlas, FACS_markers, plot.order, Ex_column, content_hierarchy, content_column, plot_colors, x_title = "expression") {
  
  # FACS_genes <- 
  #   tibble(marker = c('CCR4', 'CD3', 'CD3', 'CD3', 
  #                     'CD4', 'CD8', 'CD8', 'CD8', 
  #                     'CD11c', 'CD14', 'CD15', 'CD15', 
  #                     'CD16', 'CD16', 'CD19', 'CD20', 
  #                     'CD25', 'CD27', 'CD38', 'CD45', 
  #                     'CD45RA', 'CD56', 'CD123', 'CD127', 
  #                     'CD161', 'CD193', 'CD203', 'HLA-DR', 
  #                     'HLA-DR', 'HLA-DR', 'PTCRA'),
  #          gene_name = c('CCR4', 'CD3G', 'CD3D', 'CD3E', 
  #                        'CD4', 'CD8A', 'CD8B', 'CD8B2', 
  #                        'ITGAX', 'CD14', 'FUT4', 'FUT9', 
  #                        'FCGR3A', 'FCGR3B', 'CD19', 'MS4A1',
  #                        'IL2RA', 'CD27', 'CD38', 'PTPRC', 
  #                        'PTPRC', 'NCAM1', 'IL3RA', 'IL7R', 
  #                        'KLRB1', 'CCR3', 'ENPP3', 'HLA-DRA', 
  #                        'HLA-DRB1', 'HLA-DRB5', 'PTCRA'),
  #          comment = c('','','','',
  #                      '','','','',
  #                      '','','synthesizes CD15','synthesizes CD15',
  #                      '','','','',
  #                      '','','','',
  #                      '','','','',
  #                      '','','','',
  #                      '','','')) %>% 
  #   left_join(ensemblanno.table %>% select(1, 2, 6), by = "gene_name")
  # 
  #  
  # FACS %>% 
  #   filter(!is.na(celltype)) %>% 
  #   select(-2) %>% 
  #   gather(key = "marker", value = "expression", -1) %>% 
  #   left_join(FACS_genes, by = "marker") %>%
  #   readr::write_delim("./ref/20190123FACS cell CD markers long.txt", delim = "\t")
  
  
  
  # if(!all(unique(sample_atlas$tissue) %in% plot.order)) {
  #   missing_levels = unique(sample_atlas$tissue)[!unique(sample_atlas$tissue) %in% plot.order]
  #   warning(paste0("Missing levels in plot.order: '", paste(missing_levels, collapse = "','"), 
  #                  "'\nMissing levels added last in plot.order"))
  #   plot.order_ <- c(plot.order, missing_levels)
  # } else plot.order_ <- plot.order
  
  plot_data <- 
    sample_atlas %>% 
    rename(content_column = content_column, 
           Ex_column = Ex_column) %>%
    
    left_join(content_hierarchy, by = c("content_column" = "content")) %>%
    mutate(content_column = factor(content_column, levels = plot.order)) %>%
    filter(!is.na(content_column)) #%>%
    #left_join(FACS_markers, by = "ensg_id")
  
  
  markers <- 
    FACS_markers %>% 
    select(marker, gene_name, ensg_id, gene_description, comment) %>%
    unique()
  
  for(i in 1:nrow(markers)){
    marker <- markers$marker[i] 
    gene_name <- markers$gene_name[i]  
    ensg <- markers$ensg_id[i]
    gene_description <- markers$gene_description[i]  
    commen <- markers$comment[i]  
    
    if(is.na(ensg)) next
    
    palet <- colorRampPalette(c("red3", "white", "darkgreen"))(3)
    
    label_colors <- 
      FACS_markers %>%
      filter(ensg_id == ensg) %>%
      mutate(color = case_when(is.na(expression) ~ "black",
                               expression == "high" ~ palet[3],
                               expression == "+" ~ palet[3],
                               expression == "low" ~ palet[1],
                               expression == "-" ~ palet[1],
                               expression =="either" ~ "black")) %>%
      {.[match(plot.order, .$celltype),]} %$% 
      setNames(color, celltype)
    
    titl <- ifelse(is.na(commen), 
                   paste(unique(c(marker, gene_name)), collapse = " / "),
                   paste0(paste(unique(c(marker, gene_name)), collapse = " / "), " (", commen, ")"))
    
    
    plot_data %>%
      filter(ensg_id == ensg) %>%
      {ggplot(., aes(content_column, Ex_column, fill = content_l1))+
          
          geom_boxplot(show.legend = F) + 
          geom_rect(data = data.frame(), 
                    aes(xmin = 1:length(label_colors)-0.5, 
                        xmax = 1:length(label_colors)+0.5, 
                        ymin = -Inf, ymax = Inf), 
                    fill = ifelse(label_colors == "black", NA, label_colors), 
                    alpha = 0.3, inherit.aes = F) + 
          geom_boxplot(show.legend = F) +
          simple_theme +
          ggtitle(label = titl, subtitle = gene_description) + 
          coord_flip() + 
          scale_fill_manual(name = "", 
                            values = plot_colors) + 
          theme(#axis.text.y = element_text(color = label_colors, face = ifelse(label_colors != "black", "bold", "plain")), 
                axis.title.y = element_blank())}
    
    
    
    ggsave(paste(outpath, paste0(prefix, "_FACS_", gene_name, '_boxplot.pdf'),sep='/'), width = 10, height = 10)
    
  }
  
  
  
}

 
