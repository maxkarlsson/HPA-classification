#---
#title: "Blood Atlas Pilot"
#author: "Max Karlsson"
#created date: "2018 August 28"
#modified by:	"Per Oksvold"
#modified date: "2018 Nov 13"
#modified by: "Wen Zhong"
#modified date: "2018 Nov 28"
#---

library('impute', quietly = TRUE)
library('tidyverse', quietly = TRUE)
library('limma', quietly = TRUE)
library('NOISeq', quietly = TRUE)
library('magrittr', quietly = TRUE)

setwd('/Users/wen.zhong/Work/BloodAtlas/Bloodatlas_work/')
setwd('/Users/max.karlsson/Documents/Scilifelab/Projects/HPA-classification/')
source('./scripts/function.R')
source('./scripts/utility.R')
source('./scripts/plots.R')
source('./scripts/colors.R')
source('./scripts/classification.R')
source('./scripts/themes.R')

#Create result directory if non-existent:
result_folder <- paste0("./results/", Sys.Date())
dir.create(result_folder, showWarnings = FALSE)

## annotation
geneinfo_path <- './data/anno/gene_info.tsv'  ## coding genes list
genetable_path <- './data/anno/gene_table.txt'  ## lims geneid to ensemblid
ensemblanno_path <- './data/anno/ensembl_anno.tsv' ## ensemblid to genesymbol
contenttable_path <- './data/anno/content_table.txt' ## lims tissueid to tissue name
consensustissue_path <- './data/anno/consensus_tissue.tsv' ## lims tissue to organ
brainregions_path <- './data/anno/brain_regions.txt' ## brain_regions_anno
proteinclass_path <- './data/anno/gene.classes.txt' ## protein classification

## datasets
hpa_path <- './data/lims/rna_hpa.tsv'
gtex_path <- './data/lims/rna_gtex.tsv'
fantom_path <- './data/lims/rna_fantom.tsv'
blood_path <- './data/lims/rna_blood.tsv'

###################################################################
## Step 1. data wrangling
###################################################################

## annotation
ensg.info <- 
  geneinfo_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character()
		))

gene.table <- 
  genetable_path %>%
  readr::read_delim(delim = "\t",
                    col_types = cols(eg_id = col_character()))

ensemblanno.table <-
  ensemblanno_path %>%
  readr::read_delim(delim = "\t")

proteinclass.table <-
  proteinclass_path %>%
  readr::read_delim(delim = "\t")

content.table <-
  contenttable_path %>%
  readr::read_delim(delim = "\t",
                    col_types = cols(content_id = col_character()))

consensustissue.table <-
  consensustissue_path %>%
  readr::read_delim(delim = "\t",
                    col_types = cols(content_id = col_character(),
                                     consensus_content_id = col_character()))

brainregions.table <-
  brainregions_path %>%
  readr::read_delim(delim = "\t")

## input datasets
hpa.atlas <-
  hpa_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(lims_id = 1,
                tissue = 2,
                expression = 3)

fantom.atlas <-
  fantom_path %>%
  readr::read_delim(delim = "\t", 
		col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(lims_id = 1,
                tissue = 2,
                expression = 3) %>%
  right_join(crossing(lims_id = ensg.info$ensg_id, tissue = unique(.$tissue)), 
             by = c("lims_id", "tissue"))

gtex.atlas <-
  gtex_path %>%
  readr::read_delim(delim = "\t",
			col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(lims_id = 1,
                tissue = 2,
                expression = 3) %>%
  right_join(crossing(lims_id = ensg.info$ensg_id, tissue = unique(.$tissue)), 
             by = c("lims_id", "tissue"))

blood.atlas <-
  blood_path %>%
  readr::read_delim(delim = "\t", 
			col_types = cols(
			ensg_id = col_character(),
			tissue = col_character(),
			expression = col_double(),
			sample_type_id = col_integer()
		)) %>%
  dplyr::rename(lims_id = 1,
                tissue = 2,
                expression = 3)

## combine atlas datasets
#write("-- Combine atlas datasets", stdout())
all.atlas.raw <- 
  rbind(dplyr::select(hpa.atlas, lims_id, tissue, expression, sample_type_id), 
        dplyr::select(gtex.atlas, lims_id, tissue, expression, sample_type_id), 
        dplyr::select(fantom.atlas, lims_id, tissue, expression, sample_type_id),
        dplyr::select(blood.atlas, lims_id, tissue, expression, sample_type_id) )%>%
  mutate(method = factor(c(rep("HPA", nrow(hpa.atlas)), 
                           rep("GTEx", nrow(gtex.atlas)), 
                           rep("FANTOM", nrow(fantom.atlas)),
                           rep("Blood", nrow(blood.atlas))),
                         levels = c("Blood", "HPA", "GTEx", "FANTOM"))) %>%
  mutate(expression = round(expression, 4), 
         tissue.method = paste(tissue, method, sep = ".")) 

#readr::write_delim(all.atlas.raw, path = paste('./', "all.atlas.raw.tsv", sep = "/"), delim = "\t")


all.atlas.raw <- 
  all.atlas.raw %>%
  left_join(gene.table, by = c("lims_id" = "eg_id")) %>%
  left_join(content.table, by = c("tissue" = "content_id")) %>%
  left_join(consensustissue.table, by = c("tissue" = "content_id", "content_name" = "content_name"))

###################################################################
## Step 2. normalization
###################################################################
all.atlas <- 
  all.atlas.raw %>%
  mutate(
    # Impute missing values
    imputed.expression = impute_expression(expression, tissue.method, lims_id),
    imputed = case_when(is.na(expression) ~ TRUE,
                        TRUE ~ FALSE),
    # TMM scaling of data with imputation
    dstmm.expression = tmm_method_normalization(imputed.expression, method, tissue.method, lims_id),
    gene_dstmm.expression = pareto_scale_method_gene(dstmm.expression, method, lims_id),

    limma_gene_dstmm.expression = limma_method_correction(gene_dstmm.expression, method,
                                                          tissue.method, lims_id,
                                                          filtered.methods = "Blood"),
    limma_gene_dstmm.expression_temp = limma_gene_dstmm.expression)  %>%
  # Scale so that under limit is 1
  mutate_at(.funs = funs(. / under_limit(., expression, method)), 
            .vars = grep(".expression$", colnames(.), value = T)) 

# Remove PBMCs
all.atlas <-
 all.atlas %>%
 filter(!content_name=='total PBMC')

# readr::write_delim(all.atlas, path = paste(result_folder, paste0('all.atlas.txt'),sep='/'), delim = "\t")

###################################################################
## Step 3. Consensus
###################################################################

all.atlas.max <-
  all.atlas %>%
  # Remove genes that are imputed
  filter(!imputed) %>%
  group_by(consensus_content_name, ensg_id) %>% 
  mutate(method = as.character(method)) %>%
  dplyr::summarise_at(.funs = funs(maxEx = max(., na.rm = T),
                                   method = get_method(method, ., max(., na.rm = T))),
                      .vars = grep("expression$", colnames(.), value = T)) 

# readr::write_delim(all.atlas.max, path = paste(result_folder, paste0('all.atlas.max.txt'),sep='/'), delim = "\t")

###################################################################
## Step 4. Category
###################################################################

all.atlas.category <- get.categories.with.num.expressed(all.atlas.max,
                                                        max_column = "limma_gene_dstmm.expression_maxEx",
                                                        cat_column = "consensus_content_name",
                                                        enrich.fold = 5,
                                                        under.lim = 1,
                                                        group.num = 6)
write.table(all.atlas.category,
            file=paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'),
            row.names = F,
            sep='\t',
            quote=F)

# print number of different categories of genes
table(all.atlas.category$category.text)


###################################################################
## Blood atlas classification
###################################################################
blood.atlas <- 
  all.atlas %>%
  filter(method == "Blood")# %>% 
  # Remove Total for classification
  #filter(consensus_content_name != "total PBMC") 

blood.atlas.max <- 
  blood.atlas %>%
  group_by(content_name, ensg_id) %>% 
  filter(!is.na(limma_gene_dstmm.expression)) %>%
  dplyr::summarise(limma_gene_dstmm.expression_maxEx = max(limma_gene_dstmm.expression)) 


blood.atlas.category <- get.categories.with.num.expressed(blood.atlas.max,
                                                          max_column = "limma_gene_dstmm.expression_maxEx", 
                                                          cat_column = "content_name",
                                                          enrich.fold = 5, 
                                                          under.lim = 1, 
                                                          group.num = 11)
write.table(blood.atlas.category,
            file=paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'),
            row.names = F,
            sep='\t',
            quote=F)

# print number of different categories of genes
table(blood.atlas.category$category.text)

###################################################################
## Brain atlas classification
###################################################################
brain.atlas <- 
  all.atlas %>%
  filter(!method=='HPA') %>%
  filter(consensus_content_name == "brain")

brain.atlas <- 
  brain.atlas %>%
  left_join(brainregions.table, by=c("content_name" = "tissue.type" ,"consensus_content_name" = "group"))

## overlap genes in all tissues
brain.genelist <- 
  brain.atlas %>%
  group_by(ensg_id) %>%
  summarise(n = length(content_name)) %>%
  filter(n == 53) %>%
  select(ensg_id) %>%
  as.matrix()

brain.atlas.filter <- 
  brain.atlas[brain.atlas$ensg_id %in% brain.genelist,]


brain.atlas.max <- 
  brain.atlas.filter %>%
  group_by(subgroup, ensg_id) %>% 
  filter(!is.na(limma_gene_dstmm.expression)) %>%
  dplyr::summarise(limma_gene_dstmm.expression_maxEx = max(limma_gene_dstmm.expression)) 
# 
# write.table(brain.atlas.max,
#             file=paste(result_folder, paste0('consensus_brain_regions.txt'),sep='/'),
#             row.names = F,
#             sep='\t',
#             quote=F)

brain.atlas.category <- get.categories.with.num.expressed(brain.atlas.max,
                                                          max_column = "limma_gene_dstmm.expression_maxEx", 
                                                          cat_column = "subgroup",
                                                          enrich.fold = 5, 
                                                          group.num = 6)

# write.table(brain.atlas.category,
#             file=paste(result_folder, paste0('gene_categories_brain_regions.txt'),sep='/'),
#             row.names = F,
#             sep='\t',
#             quote=F)

# print number of different categories of genes
table(brain.atlas.category$category.text)


###################################################################
## Visulization
###################################################################
all.atlas.max.wide <- generate_wide(all.atlas.max, ensg_column='ensg_id', group_column='consensus_content_name', max_column="limma_gene_dstmm.expression_maxEx")

## tissue distribution of normalized values
make_tissue_distribution_plot(tb.atlas = all.atlas, 
                              expr_column = "limma_gene_dstmm.expression",
                              outpath = result_folder,
                              prefix = 'all_tissues')

## PCA and clustering plots
all.atlas.max.pca.values <- pca.cal(all.atlas.max.wide)
scores <- all.atlas.max.pca.values[[1]]
loadings <- as.data.frame(all.atlas.max.pca.values[[2]])
loadings$labels <- ensemblanno.table$Gene[match(rownames(loadings), ensemblanno.table$ENSG)]
tissue.colors <- with(all.tissue.table, setNames(Grouping.color, Grouping))

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(all.atlas.max.pca.values[[1]]), rownames(all.atlas.max.pca.values[[1]])),
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
                           maxEx_column = "limma_gene_dstmm.expression_maxEx", 
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues')


###### plots for brain
brain.atlas.max.wide <- generate_wide(brain.atlas.max, ensg_column='ensg_id', group_column='subgroup', max_column="limma_gene_dstmm.expression_maxEx")

## tissue distribution of normalized values
make_tissue_distribution_plot(tb.atlas = brain.atlas, 
                              expr_column = "limma_gene_dstmm.expression",
                              outpath = result_folder,
                              prefix = 'brain_regions')

## PCA and clustering plots
brain.atlas.max.pca.values <- pca.cal(brain.atlas.max.wide)
scores <- brain.atlas.max.pca.values[[1]]
loadings <- as.data.frame(brain.atlas.max.pca.values[[2]])
loadings$labels <- ensemblanno.table$gene_name[match(rownames(loadings), ensemblanno.table$ensg_id)]
tissue.colors <- with(brainregions.table, setNames(subgroup.color, tissue.type))

make_PCA_plots(scores = scores,
               loadings = loadings,
               groups = setNames(rownames(brain.atlas.max.pca.values[[1]]), rownames(brain.atlas.max.pca.values[[1]])),
               outpath = result_folder,
               prefix = 'brain_regions')

make_clustering_plot(tb.wide = brain.atlas.max.wide, 
                     colors = tissue.colors, 
                     outpath = result_folder,
                     prefix = 'brain_regions')

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
                           maxEx_column = "limma_gene_dstmm.expression_maxEx", 
                           tissue_column = "subgroup",
                           outpath = result_folder,
                           prefix = 'brain_regions')
