#---
#title: "HPA classification - raw sample processing"
#author: "Wen Zhong"
#created date: "2019 January 1st"
#---

#
# ----------- set up ----------- 
#

library('tidyverse', quietly = TRUE)
library('Biobase')


setwd('/Users/wen.zhong/Work/localGit/HPA-classification/')
setwd('/Users/max.karlsson/Documents/Scilifelab/Projects/HPA-classification/')
setwd('C:/Data/Fagerberg lab/Projects/HPA-classification')
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
tissuehierarchy_path <- './ref/colors_92.tsv' # tissue colors and hierarchy

## datasets
mouse_raw_path <- './data/lims/rna_mousebrain_mouse_sample_92.tsv'
pig_raw_path <- './data/lims/rna_pigbrain_pig_sample_92.tsv'

#
# ----------- Step 1. data wrangling ----------- 
#

## annotation

## input datasets
mouse.raw <-
  mouse_raw_path %>%
  readr::read_delim(delim = "\t", 
                    col_types = cols(
                      ensmusg_id = col_character(),
                      tissue = col_character(),
                      tpm = col_double(),
                      sample = col_character()
                    )) %>%
  dplyr::rename(ensg_id = 1,
                tissue = 2,
                sample = 3,
                expression = 4) %>%
  mutate(method = factor(c(rep("mouse")),
                         levels = c("mouse"))) %>%
  mutate(tissue.method = paste(sample, method, sep = "."))

#
# ----------- Step 2. normalization ----------- 
#

mouse.nx <- 
  mouse.raw %>%
  mutate(
    # TMM scaling of data
    dstmm.expression = tmm_method_normalization(expression, method, tissue.method, ensg_id),
    # Gene pareto
    gene_dstmm.expression = pareto_scale_method_gene(dstmm.expression, method, ensg_id))


## convert raw datasets into ExpressionSet
mouse.exprs <- as.data.frame(generate_wide(mouse.nx, ensg_column='ensg_id', group_column='sample', max_column='gene_dstmm.expression'))
mouse.exprs <- mouse.exprs[,sort(colnames(mouse.exprs))]
filter <- apply(mouse.exprs, 1, function(x) mean(x)>=1)
mouse.exprs <- mouse.exprs[filter,]
mouse.pdata <- 
  mouse.nx %>%
  select(sample=sample,tissue=tissue) %>%
  distinct(sample, tissue)
mouse.pdata <- as.data.frame(mouse.pdata)
rownames(mouse.pdata) <- mouse.pdata$sample
mouse.pdata <- mouse.pdata[sort(rownames(mouse.pdata)),]

mouse.eset <- ExpressionSet(as.matrix(mouse.exprs), AnnotatedDataFrame(mouse.pdata))

#
# ----------- Step 3. visualization ----------- 
#

make_umap_plot(eset = mouse.eset,
               outpath = result_folder,
               prefix = 'mouse_samples_gene_dstmm')
