#---
#title: "HPA classification"
#author: "Max Karlsson and Wen Zhong"
#created date: "2018 August 28"
#---

#
# ----------- set up ----------- 
#

library('impute', quietly = TRUE)
library('tidyverse', quietly = TRUE)
library('limma', quietly = TRUE)
library('NOISeq', quietly = TRUE)
library('magrittr', quietly = TRUE)

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
proteinlocalization_path <- './data/anno/new_proteinclass_all_19670.txt' ## protein localization classification
tissuehierarchy_path <- './ref/colors_92.tsv' # tissue colors and hierarchy
celllines_path <- './data/lims/Normalized/consensus_celline_hpa_92.tsv'
immunodeficiency_genes <- './ref/PID_genes_Petter_190202.txt'

blood_atlas_hierarchy <- readr::read_delim("ref/blood_atlas_hierarchy.txt", delim = "\t")
blood_atlas_colors <- readr::read_delim("ref/blood_atlas_colors.txt", delim = "\t")



## datasets
hpa_path <- './data/lims/rna_hpa.tsv'
gtex_path <- './data/lims/rna_gtex.tsv'
fantom_path <- './data/lims/rna_fantom.tsv'
blood_path <- './data/lims/rna_blood.tsv'
mouse_path <- './data/lims/rna_mousebrain_mouse_92.tsv'
pig_path <- './data/lims/rna_pigbrain_pig_92.tsv'

### normalized
hpa_norm_path <- "data/lims/Normalized/consensus_hpa_92.tsv"
gtex_norm_path <- "data/lims/Normalized/consensus_gtex_92.tsv"
fantom_norm_path <- "data/lims/Normalized/consensus_fantom_92.tsv"
blood_norm_path <- "data/lims/Normalized/consensus_bloodcells_hpa_92.tsv"


### Consensus
consensus_path <- "data/lims/Normalized/consensus_max_of_all_groups_92.tsv"
consensus_blood_path <- "data/lims/Normalized/consensus_bloodcells_hpa_92.tsv"
### category
category_path <- "data/lims/Normalized/consensus_all_category_92.tsv"

#
# ----------- Step 1. data wrangling ----------- 
#
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

proteinlocalization.table <-
  proteinlocalization_path %>%
  readr::read_delim(delim = "\t") %>%
  mutate(predicted_localization_class = case_when(class == "Predicted intracellular proteins,Predicted membrane proteins,Predicted secreted proteins" ~ "intracellular, membrane, secreted isoforms",
                                                  class == "Predicted membrane proteins,Predicted secreted proteins" ~ "membrane and secreted isoforms",
                                                  class == "Predicted intracellular proteins,Predicted secreted proteins" ~ "intracellular and secreted isoforms",
                                                  class == "Predicted intracellular proteins,Predicted membrane proteins" ~ "intracellular and membrane isoforms",
                                                  class == "Predicted membrane proteins" ~ "membrane",
                                                  class == "Predicted intracellular proteins" ~ "intracellular",
                                                  class == "Predicted secreted proteins" ~ "secreted"),
         predicted_intracellular = case_when(predicted_localization_class %in% c("intracellular, membrane, secreted isoforms",
                                                                                 "intracellular and secreted isoforms",
                                                                                 "intracellular and membrane isoforms",
                                                                                 "intracellular") ~ T,
                                             T ~ F),
         predicted_membrane = case_when(predicted_localization_class %in% c("intracellular, membrane, secreted isoforms",
                                                                            "membrane and secreted isoforms",
                                                                            "intracellular and membrane isoforms",
                                                                            "membrane") ~ T,
                                        T ~ F),
         predicted_secreted = case_when(predicted_localization_class %in% c("intracellular, membrane, secreted isoforms",
                                                                            "membrane and secreted isoforms",
                                                                            "intracellular and secreted isoforms",
                                                                            "secreted") ~ T,
                                        T ~ F))

         
immunodeficiency.table <- 
  immunodeficiency_genes %>%
  readr::read_delim(delim = "\t") %>%
  mutate(Gene = trimws(Gene),
         Gene = case_when(Gene == "IL7RA" ~ "IL7R",
                          Gene == "CD3Z" ~ "CD247",
                          Gene == "MRE11A" ~ "MRE11",
                          T ~ Gene)) %>%
  left_join(ensemblanno.table, by = c("Gene" = "gene_name")) %>%
  
  # Filter genes that are not protein coding:
  filter(!Gene %in% c("RMRP", "TERC")) %>%
  
  # Filter IgG and T-cell receptor related genes (not in Atlas):
  filter(!Gene %in% c("TRAC", "IGHM", "IGKC"))



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

contenthierarchy.table <- 
  tissuehierarchy_path %>%
  readr::read_delim(delim = "\t")

contenthierarchy.table.tissue <- contenthierarchy.table %>% filter(type=='tissue')
tissue.colors <- with(contenthierarchy.table.tissue, setNames(c(color, color, color), 
                                                              c(tissue_name, organ_name, paste(tissue_name, 1))))

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


# 
cell.lines.atlas <-
  celllines_path %>%
  readr::read_delim(delim = "\t") 
  


#
# ----------- Step 2. normalization ----------- 
#
  
blood.atlas <- 
  read_delim(blood_norm_path, delim = "\t") %>%
  mutate(method = "Blood") %>% 
  # Remove Total for classification
  filter(content_name != "total PBMC") 

blood.atlas.total.PBMC <- 
  read_delim(blood_norm_path, delim = "\t") %>%
  mutate(method = "Blood") 

blood.atlas.6 <- 
  blood.atlas %>%
  left_join(blood_atlas_hierarchy, by = c("content_name" = "content")) %>%
  mutate(content_name = content_l1)  

all.atlas <- 
  bind_rows(read_delim(hpa_norm_path, delim = "\t") %>%
              mutate(method = "HPA"),
            read_delim(gtex_norm_path, delim = "\t") %>%
              mutate(method = "GTEx"),
            read_delim(fantom_norm_path, delim = "\t") %>%
              mutate(method = "FANTOM"),
            blood.atlas)






#
# ----------- Step 3. Consensus ----------- 
#

all.atlas.max <- 
  read_delim(consensus_path, delim = "\t") %>%
  rename(consensus_content_name = 2)

# Blood Atlas 
blood.atlas.max <- 
  read_delim(consensus_blood_path, delim = "\t") %>%
  rename(consensus_content_name = 2, 
         max_norm_exp = 3)
  


#
# ----------- Step 4. Category ----------- 
#

  
all.atlas.category <- get.categories.with.num.expressed(all.atlas.max,
                                                        max_column = "max_norm_exp",
                                                        cat_column = "consensus_content_name",
                                                        enrich.fold = 4,
                                                        under.lim = 1,
                                                        group.num = 6)
readr::write_delim(all.atlas.category, path = paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'), delim = "\t")

blood.atlas.category <- get.categories.with.num.expressed(blood.atlas.max,
                                                          max_column = "max_norm_exp", 
                                                          cat_column = "consensus_content_name",
                                                          enrich.fold = 4, 
                                                          under.lim = 1, 
                                                          group.num = 11)
readr::write_delim(blood.atlas.category, path = paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'), delim = "\t")



# print number of different categories of genes
table(all.atlas.category$category.text)
all.atlas.category %>%
  group_by(elevated.category, express.category.2) %>% 
  summarise(n = length(elevated.category)) %>% 
  spread(key = express.category.2, value = n) %>%
  readr::write_delim(paste(result_folder, paste0('gene_categories_summarized.txt'),sep='/'), delim = "\t")


for(celltype in unique(blood.atlas.max$consensus_content_name)) {
  blood.atlas.category %>%
    filter(mapply(paste0("(^|, )", celltype, "(, |$)"), 
                  `enriched tissues`, FUN = function(x,y) grepl(x, y))) %>%
    write_csv(paste(result_folder, paste0(celltype, "_elevated_genes.csv"),sep='/'))
}

#
# ----------- Step 5. Cytoscape files ----------- 
#

blood.atlas.category.cytoscape.nodes <- 
  blood.atlas.category %>% 
  group_by(elevated.category, `enriched tissues`) %>% 
  summarise(n = length(ensg_id)) %>%
  ungroup() %>%
  mutate(node_id = 1:nrow(.))

first <- T
for(content_name in unique(blood.atlas$content_name)) {
  temp <- 
    blood.atlas.category.cytoscape.nodes %>%
    filter(grepl(paste0("(^|, )", content_name, "(, |$)"), `enriched tissues`)) %>% 
    mutate(content_name = content_name, 
           edge = 1)
  
  if(first){
    blood.atlas.category.cytoscape.nodes.full <- temp
    
    first <- F
  } else {
    blood.atlas.category.cytoscape.nodes.full <- 
      rbind(blood.atlas.category.cytoscape.nodes.full, temp)
  }
  
  
}

blood.atlas.category.cytoscape.nodes.full %>%
  left_join(blood_atlas_hierarchy %>% 
              filter(!content %in% c("blood", "Total PBMCs")) %>% 
              {.[with(., order(content_l3, content_l2, content_l1, content)),]} %>% 
              mutate(circle_order_1 = 1:nrow(.), 
                     circle_order_2 = 1:nrow(.),
                     content_name = content), 
            by = "content_name") %>%
  filter(n > 2 & elevated.category %in% c("group enriched", "tissue enriched")) %>%
  write_delim(paste(result_folder, paste0('gene_categories_blood_cells_summarised.txt'),sep='/'), delim = "\t")

#
# ----------- Brain atlas classification ----------- 
#
# 
# brain.atlas <- 
#   all.atlas %>%
#   filter(!method=='HPA') %>%
#   filter(consensus_content_name == "brain")
# 
# brain.atlas <- 
#   brain.atlas %>%
#   left_join(brainregions.table, by=c("content_name" = "tissue.type" ,"consensus_content_name" = "group"))
# 
# ## overlap genes in all tissues
# brain.genelist <- 
#   brain.atlas %>%
#   group_by(ensg_id) %>%
#   summarise(n = length(content_name)) %>%
#   filter(n == 53) %>%
#   select(ensg_id) %>%
#   as.matrix()
# 
# brain.atlas.filter <- 
#   brain.atlas[brain.atlas$ensg_id %in% brain.genelist,]
# 
# 
# brain.atlas.max <- 
#   brain.atlas.filter %>%
#   group_by(subgroup, ensg_id) %>% 
#   filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
#   dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 
# 
# brain.atlas.max_all_regions <- 
#   brain.atlas.filter %>% 
#   group_by(content_name, ensg_id) %>% 
#   filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
#   dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 
# 
# write.table(brain.atlas.max,
#             file=paste(result_folder, paste0('consensus_brain_regions.txt'),sep='/'),
#             row.names = F,
#             sep='\t',
#             quote=F)
# 
# brain.atlas.category <- get.categories.with.num.expressed(brain.atlas.max,
#                                                           max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                                                           cat_column = "subgroup",
#                                                           enrich.fold = 4, 
#                                                           group.num = 6)
# 
# brain.atlas.category_all_regions <-
#   get.categories.with.num.expressed(brain.atlas.max_all_regions,
#                                     max_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                     cat_column = "content_name",
#                                     enrich.fold = 4,
#                                     group.num = 6)
# write.table(brain.atlas.category,
#             file=paste(result_folder, paste0('gene_categories_brain_regions.txt'),sep='/'),
#             row.names = F,
#             sep='\t',
#             quote=F)
# 
# # print number of different categories of genes
# table(brain.atlas.category$category.text)


#
# ----------- Visulization ----------- 
#



# =========== *All altas =========== 
plots <- c("",
           "all",
           # "spearman dendrogram",
           # "tissue distribution",
           # "PCA",
           # "cluster",
           # "elevated bar",
           # "specificity distribution",
           # "class chord",
           # "group chord",
           # "heatmaps",
           # "swarm expression",
           # "number detected bar",
           # "NX fraction bar",
           # "classification pie",
           # "TPM NX example genes bar",
           # "score plots",
           "sum TPM",
           # "class comparison chord",
           # "class comparison chord",
           # "spearman dendrogram",
           "blood class tissue expression")

make_plots(atlas = all.atlas, 
           atlas.max = all.atlas.max, 
           atlas.cat = all.atlas.category, 
           Ex_column = "limma_gene_dstmm.zero.impute.expression", 
           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",  
           content_column = "consensus_content_name", 
           content_hierarchy = contenthierarchy.table.tissue %>% 
             select(1:2) %>% 
             rename(content = tissue_name, content_l1 = organ_name),  
           content_colors = tissue.colors, 
           plots = plots, 
           plot.atlas = "tissue", 
           plot.order = unique(all.atlas.max$consensus_content_name),
           subatlas_unit = "tissue",
           outpath = result_folder, 
           prefix = "tissue")

#####
atlas = all.atlas
atlas.max = all.atlas.max
atlas.cat = all.atlas.category
Ex_column = "limma_gene_dstmm.zero.impute.expression"
maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx"
content_column = "content_name"
consensus_content_column = "consensus_content_name"
content_hierarchy = NULL
content_colors = tissue.colors
plots = plots
plot.atlas = "tissue"
plot.order = unique(all.atlas.max$consensus_content_name)
subatlas_unit = "tissue"
outpath = result_folder
prefix = "tissue"
#####


# =========== *Brain altas* =========== 



# =========== *Blood altas* =========== 
plots <- c("",
           # "spearman dendrogram",
           # "tissue distribution",
           # "PCA",
           # "cluster",
           # "elevated bar",
           # "specificity distribution",
           # "class chord",
           # "group chord",
           # "heatmaps",
           # "swarm expression",
           # "number detected bar",
           # "NX fraction bar",
           # "classification pie",
           # "TPM NX example genes bar",
           # "score plots",
           "sum TPM",
           # "class comparison chord",
           # "class comparison chord",
           # "spearman dendrogram",
           "blood class tissue expression", 
           "double donut chord")

make_plots(atlas = blood.atlas, 
           atlas.max = blood.atlas.max, 
           atlas.cat = blood.atlas.category, 
           Ex_column = "limma_gene_dstmm.zero.impute.expression", 
           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",  
           content_column = "content_name", 
           content_hierarchy = blood_atlas_hierarchy, 
           content_colors = with(blood_atlas_colors, setNames(color, content)), 
           plots = plots, 
           plot.atlas = "blood", 
           plot.order = blood_atlas_hierarchy %>%
             filter(!content %in% c("blood", "Total PBMCs")) %$% 
             content[order(content_l1)],
           subatlas_unit = "celltype",
           outpath = result_folder, 
           prefix = "blood_cells")


###
atlas = blood.atlas 
atlas.max = blood.atlas.max 
atlas.cat = blood.atlas.category 
Ex_column = "limma_gene_dstmm.zero.impute.expression" 
maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx"  
content_column = "content_name" 
content_hierarchy = blood_atlas_hierarchy 
content_colors = with(blood_atlas_colors, setNames(color, content))
plots = "all" 
plot.atlas = "blood" 
plot.order = blood_atlas_hierarchy %>%
  filter(!content %in% c("blood", "Total PBMCs")) %$% 
  content[order(content_l1)]
subatlas_unit = "celltype"
outpath = result_folder 
prefix = "blood_cells"
###

# =========== *Blood altas (6 cells)* =========== 



# 

