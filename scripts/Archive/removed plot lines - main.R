# =========== *All altas ===========

make_sum_TPM_plot(all.atlas, all.atlas.category, tissue_column = "content_name", method_column = "method", outpath = result_folder, prefix = "atlas")

make_classification_pie_chart(atlas.cat = all.atlas.category,
                              outpath = result_folder,
                              prefix = "all_atlas")

all.atlas.max.wide <- generate_wide(all.atlas.max, ensg_column='ensg_id',
                                    group_column='consensus_content_name',
                                    max_column="limma_gene_dstmm.zero.impute.expression_maxEx")


# Tissue distribution
tissues_to_plot <- c('tongue', 'thalamus', 'skeletal muscle', 'thyroid gland', 'appendix', 'esophagus',
                     'postcentral gyrus', 'liver', 'Myeloid DCs', 'olfactory region', 'Memory CD8 T-cells',
                     'pancreas', 'hypothalamus', 'Memory B-cells', 'spleen', 'duodenum', 'rectum', 'cerebellum',
                     'prostate', 'Naive B-cells', 'vagina', 'endometrium 1', 'heart muscle', 'adrenal gland',
                     'skin 1', 'salivary gland', 'thymus', 'small intestine', 'frontal lobe', 'kidney',
                     'gallbladder', 'putamen', 'cervix, uterine', 'Eosinophils', 'seminal vesicle', 'pons',
                     'placenta', 'ductus deferens', 'amygdala')

make_tissue_distributions_plot(atlas.tb = all.atlas,
                               Ex_column = "limma_gene_dstmm.zero.impute.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = tissues_to_plot,
                               outpath = result_folder,
                               prefix = "All atlas NX")

make_tissue_distributions_plot(atlas.tb = all.atlas,
                               Ex_column = "expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = tissues_to_plot,
                               outpath = result_folder,
                               prefix = "X")

make_tissue_distributions_plot(atlas.tb = all.atlas,
                               Ex_column = "dstmm.zero.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = tissues_to_plot,
                               outpath = result_folder,
                               prefix = "All atlas TMM")

make_tissue_distributions_plot(atlas.tb = all.atlas,
                               Ex_column = "gene_dstmm.zero.impute.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = tissues_to_plot,
                               outpath = result_folder,
                               prefix = "All atlas TMM pareto")

## Spearman method cluster
make_spearman_method_dendrogram(all.atlas.tb = all.atlas,
                                Ex_column = "limma_gene_dstmm.zero.impute.expression",
                                content_column = "content_name",
                                named_color_replacement = dataset.colors,
                                outpath = result_folder,
                                prefix = "All_atlas_norm_method_color")

make_spearman_method_dendrogram(all.atlas.tb = all.atlas,
                                Ex_column = "limma_gene_dstmm.zero.impute.expression",
                                content_column = "content_name",
                                named_color_replacement = tissue.colors,
                                outpath = result_folder,
                                prefix = "All_atlas_norm_tissue_color")

make_spearman_method_dendrogram(all.atlas.tb = all.atlas,
                                Ex_column = "expression",
                                content_column = "content_name",
                                named_color_replacement = dataset.colors,
                                outpath = result_folder,
                                prefix = "All_atlas_exp_method_color")

make_spearman_method_dendrogram(all.atlas.tb = all.atlas,
                                Ex_column = "expression",
                                content_column = "content_name",
                                named_color_replacement = tissue.colors,
                                outpath = result_folder,
                                prefix = "All_atlas_exp_tissue_color")

## tissue distribution of normalized values
make_tissue_distribution_plot(tb.atlas = all.atlas,
                              expr_column = "limma_gene_dstmm.zero.impute.expression",
                              outpath = result_folder,
                              prefix = 'all_tissues')

## PCA and clustering plots
all.atlas.max.pca.values <- pca.cal(all.atlas.max.wide)
scores <- all.atlas.max.pca.values[[1]]
loadings <-
  all.atlas.max.pca.values[[2]] %>%
  as.tibble(rownames = "ensg_id") %>%
  mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(all.atlas.max.pca.values[[1]]), rownames(all.atlas.max.pca.values[[1]])),
               groups.color = tissue.colors,
               outpath = result_folder,
               prefix = 'all_tissues')

make_clustering_plot(tb.wide = all.atlas.max.wide,
                     colors = tissue.colors,
                     outpath = result_folder,
                     prefix = 'all_tissues')

## tissue elevated plot
all.atlas.elevated.table <- calc_elevated.table(tb.wide = all.atlas.max.wide,
                                                atlas.categories = all.atlas.category)
all.atlas.elevated.summary.table <- calc_elevated.summary.table(all.atlas.elevated.table)
make_elevated_bar_plot(elevated.summary.table = all.atlas.elevated.summary.table,
                       outpath = result_folder,
                       prefix = 'all_tissues')

## specificity distribution
make_specificity_distribution_plot(atlas.cat = all.atlas.category,
                                   type = "Tissue",
                                   outpath = result_folder,
                                   prefix = 'all_tissues')

## chord plot
make_classification_chord_plot(atlas.cat = all.atlas.category,
                               outpath = result_folder,
                               prefix = 'all_tissues')



## swarm plot
make_swarm_expression_plot(atlas.max = all.atlas.max,
                           atlas.cat = all.atlas.category,
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues')

# group enriched chord diagram
all_atlas_hierarchy <-
  contenthierarchy.table.tissue %>%
  select(1:2) %>%
  rename(content = 1, content_l1 = 2)
make_chord_group_enriched(all.atlas.elevated.table,
                          grid.col = tissue.colors,
                          #tissue_hierarcy = all_atlas_hierarchy,
                          tissue_hierarcy = rbind(all_atlas_hierarchy, mutate(all_atlas_hierarchy, content = paste(content, 1))),
                          palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
                          outpath = result_folder,
                          prefix = "all_atlas")

make_heatmap_group_enriched(all.atlas.elevated.table,
                            outpath = result_folder,
                            prefix = "all_atlas")

make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
                                                     all.atlas.max.tb = all.atlas.max,
                                                     maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                                     tissue_column = "consensus_content_name",
                                                     outpath = result_folder,
                                                     prefix = "all_atlas")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = all.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         prefix = "all atlas")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = all.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         range_scale_x = T,
                         prefix = "all atlas range scaled")



# make_heatmap_all_elevated_expression_levels_circle(elevated.table = all.atlas.elevated.table,
#                                                    all.atlas.max.tb = all.atlas.max,
#                                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                    tissue_column = "consensus_content_name",
#                                                    outpath = result_folder,
#                                                    prefix = "all_atlas")

# make_heatmap_group_and_enhanced_expression_levels_circle(elevated.table = all.atlas.elevated.table,
#                                                          all.atlas.max.tb = all.atlas.max,
#                                                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                          tissue_column = "consensus_content_name",
#                                                          outpath = result_folder,
#                                                          prefix = "all_atlas")

# make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
#                                                      all.atlas.max.tb = all.atlas.max,
#                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                      tissue_column = "consensus_content_name",
#                                                      outpath = result_folder,
#                                                      prefix = "all_atlas",
#                                                      y_dendrogram = F)
#
# make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
#                                                      all.atlas.max.tb = all.atlas.max,
#                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                      tissue_column = "consensus_content_name",
#                                                      outpath = result_folder,
#                                                      prefix = "all_atlas_dendro",
#                                                      y_dendrogram = T)

# Number of expressed genes
make_number_detected_genes_barplot(all.atlas.max.tb = all.atlas.max,
                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                   tissue_column = "consensus_content_name",
                                   outpath = result_folder,
                                   prefix = "all_atlas")

# Total elevated expression fraction
make_elevated_NX_fraction_barplots(atlas.max = all.atlas.max,
                                   atlas.cat = all.atlas.category,
                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                   tissue_column = "consensus_content_name",
                                   outpath = result_folder,
                                   prefix = "all_atlas")

# =========== *Brain altas* ===========

brain.atlas.max.wide <- generate_wide(brain.atlas.max, ensg_column='ensg_id', group_column='subgroup',
                                      max_column="limma_gene_dstmm.zero.impute.expression_maxEx")

brain.atlas.max.wide_all_regions <- generate_wide(brain.atlas.max_all_regions, ensg_column='ensg_id',
                                                  group_column='content_name',
                                                  max_column="limma_gene_dstmm.zero.impute.expression_maxEx")

make_classification_pie_chart(atlas.cat = brain.atlas.category,
                              outpath = result_folder,
                              prefix = "brain_atlas")


## tissue distribution of normalized values
make_tissue_distribution_plot(tb.atlas = brain.atlas,
                              expr_column = "limma_gene_dstmm.zero.impute.expression",
                              outpath = result_folder,
                              prefix = 'brain_regions')

## PCA and clustering plots
brain.atlas.max.pca.values <- pca.cal(brain.atlas.max.wide)
scores <- brain.atlas.max.pca.values[[1]]
loadings <-
  brain.atlas.max.pca.values[[2]] %>%
  as.tibble(rownames = "ensg_id") %>%
  mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
tissue.colors <- with(brainregions.table, setNames(subgroup.color, subgroup))

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(brain.atlas.max.pca.values[[1]]), rownames(brain.atlas.max.pca.values[[1]])),
               groups.color = tissue.colors,
               outpath = result_folder,
               prefix = 'brain_regions')

make_clustering_plot(tb.wide = brain.atlas.max.wide,
                     colors = tissue.colors,
                     outpath = result_folder,
                     prefix = 'brain_regions')

cell.colors <- with(brainregions.table, setNames(subgroup.color, tissue.type))
make_clustering_plot(tb.wide = brain.atlas.max.wide_all_regions,
                     colors = cell.colors,
                     outpath = result_folder,
                     prefix = 'brain_all_cells')

## tissue elevated plot
brain.atlas.elevated.table <- calc_elevated.table(tb.wide = brain.atlas.max.wide,
                                                  atlas.categories = brain.atlas.category)
brain.atlas.elevated.summary.table <- calc_elevated.summary.table(brain.atlas.elevated.table)
make_elevated_bar_plot(elevated.summary.table = brain.atlas.elevated.summary.table,
                       outpath = result_folder,
                       prefix = 'brain_regions')

## specificity distribution
make_specificity_distribution_plot(atlas.cat = brain.atlas.category,
                                   type = "Tissue",
                                   outpath = result_folder,
                                   prefix = 'brain_regions')

## chord plot
make_classification_chord_plot(atlas.cat = brain.atlas.category,
                               outpath = result_folder,
                               prefix = 'brain_tissues')

## swarm plot
make_swarm_expression_plot(atlas.max = brain.atlas.max,
                           atlas.cat = brain.atlas.category,
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                           tissue_column = "subgroup",
                           outpath = result_folder,
                           prefix = 'brain_regions')

# group enriched chord diagram
#brain_atlas_hierarchy <- readr::read_delim("ref/brain_atlas_hierarchy.txt", delim = "\t")
contenthierarchy.table.brain <- contenthierarchy.table %>% filter(type=='brain')
tissue.colors.brain <- with(contenthierarchy.table.brain, setNames(c(color, color, color), c(tissue_name, organ_name, paste(organ_name, 1))))

brain.atlas.elevated.table_all_regions <-
  calc_elevated.table(tb.wide = brain.atlas.max.wide_all_regions,
                      atlas.categories = brain.atlas.category_all_regions)

brain.atlas.elevated.table <-
  calc_elevated.table(tb.wide = brain.atlas.max.wide,
                      atlas.categories = brain.atlas.category)

brain_atlas_hierarchy <-
  contenthierarchy.table.brain %>%
  select(content=organ_name)

make_chord_group_enriched(brain.atlas.elevated.table,
                          grid.col = tissue.colors.brain,
                          tissue_hierarcy = rbind(brain_atlas_hierarchy, mutate(brain_atlas_hierarchy, content = paste(content, 1))),
                          palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
                          outpath = result_folder,
                          prefix = "brain_atlas")

make_heatmap_group_enriched(brain.atlas.elevated.table_all_regions,
                            outpath = result_folder,
                            prefix = "brain_atlas_all_regions")

make_heatmap_median_expression_levels(elevated.table = brain.atlas.elevated.table,
                                      all.atlas.max.tb = brain.atlas.max,
                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                      tissue_column = "subgroup",
                                      enrichment = c(3),
                                      outpath = result_folder,
                                      prefix = "brain_atlas_group_enriched")

make_heatmap_median_expression_levels(elevated.table = brain.atlas.elevated.table,
                                      all.atlas.max.tb = brain.atlas.max,
                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                      tissue_column = "subgroup",
                                      enrichment = c(2, 3, 4),
                                      outpath = result_folder,
                                      prefix = "brain_atlas_all_elevated")

make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
                               all.atlas.max.tb = brain.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "subgroup",
                               enrichment = c(3),
                               outpath = result_folder,
                               prefix = "brain_atlas_group_enriched")

make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
                               all.atlas.max.tb = brain.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "subgroup",
                               enrichment = c(2,3,4),
                               outpath = result_folder,
                               prefix = "brain_atlas_all_elevated")

make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
                               all.atlas.max.tb = brain.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "subgroup",
                               enrichment = c(2),
                               outpath = result_folder,
                               prefix = "brain_atlas_tissue_enriched")

# =========== *Blood altas* ===========

blood.atlas.max.wide <- generate_wide(blood.atlas.max, ensg_column='ensg_id', group_column='content_name',
                                      max_column="limma_gene_dstmm.zero.impute.expression_maxEx")

make_classification_pie_chart(atlas.cat = blood.atlas.category,
                              outpath = result_folder,
                              prefix = "blood_atlas")

blood.atlas.category %>%
  left_join(proteinclass.table, by = c("ensg_id" = "rna.genes")) %>%
  filter(category %in% 2:4) %$%
  table(proteinclass.vec.single)


blood.atlas.max %>%
  filter(ensg_id %in% unique(ensg_id)[1:100]) %>%
  make_gene_expression_barplot(maxEx_columns = c("Raw" = "expression_maxEx", "TMM" = "dstmm.zero.expression_maxEx", "TMM + Pareto" = "gene_dstmm.zero.impute.expression_maxEx"),
                               content_column = "content_name",
                               content_color = with(blood_atlas_colors, setNames(color, content)))

# Tissue distribution
make_tissue_distributions_plot(atlas.tb = blood.atlas,
                               Ex_column = "limma_gene_dstmm.zero.impute.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = "all",
                               outpath = result_folder,
                               prefix = "Blood atlas NX")

make_tissue_distributions_plot(atlas.tb = blood.atlas,
                               Ex_column = "expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = "all",
                               outpath = result_folder,
                               prefix = "Blood atlas X")

make_tissue_distributions_plot(atlas.tb = blood.atlas,
                               Ex_column = "dstmm.zero.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = "all",
                               outpath = result_folder,
                               prefix = "Blood atlas TMM")

make_tissue_distributions_plot(atlas.tb = blood.atlas,
                               Ex_column = "gene_dstmm.zero.impute.expression",
                               content_column = "content_name",
                               und.lim = 1,
                               do.tissues = "all",
                               outpath = result_folder,
                               prefix = "Blood atlas TMM pareto")

## tissue distribution of normalized values
make_tissue_distribution_plot(tb.atlas = blood.atlas,
                              expr_column = "limma_gene_dstmm.zero.impute.expression",
                              outpath = result_folder,
                              prefix = 'blood_cells')

## PCA and clustering plots
blood.atlas.max.pca.values <- pca.cal(blood.atlas.max.wide)
scores <- blood.atlas.max.pca.values[[1]]
loadings <-
  blood.atlas.max.pca.values[[2]] %>%
  as.tibble(rownames = "ensg_id") %>%
  mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
tissue.colors <- with(contenthierarchy.table, setNames(c(color, color, color), c(tissue_name, organ_name, paste(tissue_name, 1))))

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(blood.atlas.max.pca.values[[1]]), rownames(blood.atlas.max.pca.values[[1]])),
               groups.color = tissue.colors,
               outpath = result_folder,
               prefix = 'blood_celltypes')

make_clustering_plot(tb.wide = blood.atlas.max.wide,
                     colors = tissue.colors,
                     outpath = result_folder,
                     prefix = 'blood_celltypes')

## tissue elevated plot
blood.atlas.elevated.table <- calc_elevated.table(tb.wide = blood.atlas.max.wide,
                                                  atlas.categories = blood.atlas.category,
                                                  cat.colum = "category")
blood.atlas.elevated.summary.table <- calc_elevated.summary.table(blood.atlas.elevated.table)
make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table,
                       outpath = result_folder,
                       translate_categories = c("Celltype" = "Tissue"),
                       prefix = 'blood_celltypes')

make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table,
                       outpath = result_folder,
                       translate_categories = c("Celltype" = "Tissue"),
                       prefix = 'blood_celltypes')

make_score_expression_scatter(atlas.max.tb = blood.atlas.max,
                              atlas.cat = blood.atlas.category,
                              maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                              tissue_column = "content_name",
                              ensemblanno.table = ensemblanno.table,
                              plot.order = blood_atlas_hierarchy %>%
                                filter(!content %in% c("blood", "Total PBMCs")) %$%
                                content[order(content_l1)],
                              outpath = result_folder,
                              prefix = "blood_celltypes")

## specificity distribution
make_specificity_distribution_plot(atlas.cat = blood.atlas.category,
                                   type = "Tissue",
                                   outpath = result_folder,
                                   prefix = 'blood_celltypes')

## chord plot
make_classification_chord_plot(atlas.cat = blood.atlas.category,
                               outpath = result_folder,
                               prefix = 'blood_tissues')

## swarm plot
make_swarm_expression_plot(atlas.max = blood.atlas.max,
                           atlas.cat = blood.atlas.category,
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                           tissue_column = "content_name",
                           plot.order = blood_atlas_hierarchy %>%
                             filter(!content %in% c("blood", "Total PBMCs")) %$%
                             content[order(content_l1)],
                           outpath = result_folder,
                           prefix = 'blood_celltypes')


make_swarm_expression_circle_plot(atlas.max = blood.atlas.max,
                                  atlas.cat = blood.atlas.category,
                                  maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                  tissue_column = "content_name",
                                  plot.order = blood_atlas_hierarchy %>%
                                    filter(!content %in% c("blood", "Total PBMCs")) %$%
                                    content[order(content_l1)],
                                  outpath = result_folder,
                                  prefix = 'blood_celltypes')
# group enriched chord diagram
blood_atlas_hierarchy <- readr::read_delim("ref/blood_atlas_hierarchy.txt", delim = "\t")
blood_atlas_colors <- readr::read_delim("ref/blood_atlas_colors.txt", delim = "\t")

make_chord_group_enriched(blood.atlas.elevated.table,
                          grid.col = with(blood_atlas_colors, setNames(color, content)),
                          tissue_hierarcy = blood_atlas_hierarchy,
                          palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
                          outpath = result_folder,reverse = T,
                          prefix = "blood_atlas")

make_heatmap_group_enriched(elevated.table = blood.atlas.elevated.table,
                            outpath = result_folder,
                            prefix = "blood_atlas")

make_expression_heatmaps(atlas.max.tb = blood.atlas.max,
                         atlas.cat = blood.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         prefix = "blood atlas")

make_expression_heatmaps(atlas.max.tb = blood.atlas.max,
                         atlas.cat = blood.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         range_scale_x = T,
                         prefix = "blood atlas range scaled")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = blood.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         prefix = "blood atlas cat on all atlas")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = blood.atlas.category,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         range_scale_x = T,
                         prefix = "blood atlas cat on all atlas range scaled")


make_immunodeficiency_expression_heatmaps(atlas.max.tb = blood.atlas.max,
                                          atlas.cat = blood.atlas.category,
                                          immunodeficiency.table = immunodeficiency.table,
                                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                          tissue_column = "content_name",
                                          ensemblanno.table = ensemblanno.table,
                                          proteinclass.table = proteinclass.table,
                                          proteinclass.table_ensg_id_column = "rna.genes",
                                          proteinclass.table_class_column = "proteinclass.vec.single",
                                          outpath = result_folder,
                                          range_scale_x = F,
                                          prefix = "blood atlas immuno deficiency")

make_immunodeficiency_expression_heatmaps(atlas.max.tb = blood.atlas.max,
                                          atlas.cat = blood.atlas.category,
                                          immunodeficiency.table = immunodeficiency.table,
                                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                          tissue_column = "content_name",
                                          ensemblanno.table = ensemblanno.table,
                                          proteinclass.table = proteinclass.table,
                                          proteinclass.table_ensg_id_column = "rna.genes",
                                          proteinclass.table_class_column = "proteinclass.vec.single",
                                          outpath = result_folder,
                                          range_scale_x = T,
                                          prefix = "blood atlas immuno deficiency range scaled")

# make_heatmap_group_enriched_expression_levels_circle(elevated.table = blood.atlas.elevated.table,
#                                                      all.atlas.max.tb = blood.atlas.max,
#                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                      tissue_column = "content_name",
#                                                      outpath = result_folder,
#                                                      prefix = "blood_atlas")

make_heatmap_median_expression_levels(elevated.table = blood.atlas.elevated.table,
                                      all.atlas.max.tb = blood.atlas.max,
                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                      tissue_column = "content_name",
                                      enrichment = c(3),
                                      outpath = result_folder,
                                      prefix = "blood_atlas_group_enriched")

make_heatmap_median_expression_levels(elevated.table = blood.atlas.elevated.table,
                                      all.atlas.max.tb = blood.atlas.max,
                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                      tissue_column = "content_name",
                                      enrichment = c(2, 3, 4),
                                      outpath = result_folder,
                                      prefix = "blood_atlas_all_elevated")

make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
                               all.atlas.max.tb = blood.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "content_name",
                               enrichment = c(3),
                               outpath = result_folder,
                               prefix = "blood_atlas_group_enriched")

make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
                               all.atlas.max.tb = blood.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "content_name",
                               enrichment = c(2,3,4),
                               outpath = result_folder,
                               prefix = "blood_atlas_all_elevated")

make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
                               all.atlas.max.tb = blood.atlas.max,
                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                               tissue_column = "content_name",
                               enrichment = c(2),
                               outpath = result_folder,
                               prefix = "blood_atlas_tissue_enriched")

# Categories between blood and all atlas
make_class_comparison_chord(cat1 = blood.atlas.category,
                            cat2 = all.atlas.category,
                            outpath = result_folder, prefix = "blood")

# Number of expressed genes
make_number_detected_genes_barplot(all.atlas.max.tb = blood.atlas.max,
                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                   tissue_column = "content_name",
                                   outpath = result_folder,
                                   prefix = "blood_atlas")

# Comparison of elevated genes to tissue atlas
make_elevated_organ_total_chord(cat1 = blood.atlas.category,
                                cat2 = all.atlas.category,
                                grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))),
                                elevated_cats = c(2,3,4),
                                direction = 1,
                                cat1_name = "celltypes",
                                cat2_name = "tissues",
                                outpath = result_folder,
                                prefix = "Blood to tissue elevated")


make_elevated_organ_total_chord(cat1 = blood.atlas.category,
                                cat2 = all.atlas.category,
                                grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))),
                                elevated_cats = c(2,3,4),
                                direction = 2,
                                cat1_name = "celltypes",
                                cat2_name = "tissues",
                                outpath = result_folder,
                                prefix = "Tissue to Blood elevated")

make_elevated_organ_total_chord(cat1 = blood.atlas.category,
                                cat2 = all.atlas.category,
                                grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))),
                                elevated_cats = c(2),
                                direction = 1,
                                cat1_name = "celltypes",
                                cat2_name = "tissues",
                                outpath = result_folder,
                                prefix = "Blood to tissue tissue enriched")

make_elevated_organ_total_chord(cat1 = blood.atlas.category,
                                cat2 = all.atlas.category,
                                grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))),
                                elevated_cats = c(2),
                                direction = 2,
                                cat1_name = "celltypes",
                                cat2_name = "tissues",
                                outpath = result_folder,
                                prefix = "Tissue to Blood tissue enriched")

# Total elevated expression fraction
make_elevated_NX_fraction_barplots(atlas.max = blood.atlas.max,
                                   atlas.cat = blood.atlas.category,
                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                   tissue_column = "content_name",
                                   outpath = result_folder,
                                   prefix = "blood_atlas")

# =========== *Blood altas (6 cells)* ===========

blood.atlas.max.wide.6 <- generate_wide(blood.atlas.max.6, ensg_column='ensg_id', group_column='content_name',
                                        max_column="limma_gene_dstmm.zero.impute.expression_maxEx")

make_classification_pie_chart(atlas.cat = blood.atlas.category.6,
                              outpath = result_folder,
                              prefix = "blood_atlas_6")


## PCA and clustering plots
blood.atlas.max.pca.values.6 <- pca.cal(blood.atlas.max.wide.6)
scores <- blood.atlas.max.pca.values.6[[1]]
loadings <-
  blood.atlas.max.pca.values.6[[2]] %>%
  as.tibble(rownames = "ensg_id") %>%
  mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(blood.atlas.max.pca.values.6[[1]]), rownames(blood.atlas.max.pca.values.6[[1]])),
               groups.color = with(blood_atlas_colors, setNames(color, content)),
               outpath = result_folder,
               prefix = 'blood_celltypes_6')

make_clustering_plot(tb.wide = blood.atlas.max.wide.6,
                     colors = with(blood_atlas_colors, setNames(color, content)),
                     outpath = result_folder,
                     prefix = 'blood_celltypes_6')


## tissue elevated plot
blood.atlas.elevated.table.6 <- calc_elevated.table(tb.wide = blood.atlas.max.wide.6,
                                                    atlas.categories = blood.atlas.category.6)
blood.atlas.elevated.summary.table.6 <- calc_elevated.summary.table(blood.atlas.elevated.table.6)
make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table.6,
                       outpath = result_folder,
                       prefix = 'blood_celltypes_6')

## specificity distribution
make_specificity_distribution_plot(atlas.cat = blood.atlas.category.6,
                                   type = "Tissue",
                                   outpath = result_folder,
                                   prefix = 'blood_celltypes_6')


make_swarm_expression_plot(atlas.max = blood.atlas.max.6,
                           atlas.cat = blood.atlas.category.6,
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                           tissue_column = "content_name",
                           plot.order = blood_atlas_hierarchy %>%
                             filter(!content %in% c("blood", "Total PBMCs")) %$%
                             unique(content_l1[order(content_l2)]),
                           outpath = result_folder,
                           prefix = 'blood_cell_lineages_6')


## chord plot
# make_classification_chord_plot(atlas.cat = blood.atlas.category.6,
#                                outpath = result_folder,
#                                prefix = 'blood_tissues_6')

# group enriched chord diagram


make_chord_group_enriched(blood.atlas.elevated.table.6,
                          grid.col = with(blood_atlas_colors, setNames(color, content)),
                          tissue_hierarcy = blood_atlas_hierarchy %>%
                            select(content_l1, content_l2, content_l3) %>%
                            rename(content = content_l1,
                                   content_l1 = content_l2,
                                   content_l2 = content_l3),
                          palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
                          outpath = result_folder,reverse = T,
                          prefix = "blood_atlas_6")

make_heatmap_group_enriched(elevated.table = blood.atlas.elevated.table.6,
                            outpath = result_folder,
                            prefix = "blood_atlas_6")


make_expression_heatmaps(atlas.max.tb = blood.atlas.max.6,
                         atlas.cat = blood.atlas.category.6,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         prefix = "blood atlas 6")

make_expression_heatmaps(atlas.max.tb = blood.atlas.max.6,
                         atlas.cat = blood.atlas.category.6,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         range_scale_x = T,
                         prefix = "blood atlas 6 range scaled")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = blood.atlas.category.6,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         prefix = "blood atlas 6 cat on all atlas")

make_expression_heatmaps(atlas.max.tb = all.atlas.max,
                         atlas.cat = blood.atlas.category.6,
                         maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                         tissue_column = "consensus_content_name",
                         ensemblanno.table = ensemblanno.table,
                         proteinclass.table = proteinclass.table,
                         proteinclass.table_ensg_id_column = "rna.genes",
                         proteinclass.table_class_column = "proteinclass.vec.single",
                         outpath = result_folder,
                         range_scale_x = T,
                         prefix = "blood atlas 6 cat on all atlas range scaled")

# Total elevated expression fraction
make_elevated_NX_fraction_barplots(atlas.max = blood.atlas.max.6,
                                   atlas.cat = blood.atlas.category.6,
                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                   tissue_column = "content_name",
                                   outpath = result_folder,
                                   prefix = "blood_atlas_6")