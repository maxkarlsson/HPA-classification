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
celllines_path <- './data/lims/consensus_celline_hpa_92.tsv'
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

if(!file.exists(paste(result_folder, paste0('all.atlas.txt'),sep='/'))) {
  all.atlas <- 
    all.atlas.raw %>%
    mutate(
      # Impute missing values
      imputed = case_when(is.na(expression) ~ TRUE,
                          TRUE ~ FALSE),
      # TMM scaling of data with imputation (set to 0)
      
      imputed.zero.expression = ifelse(imputed, 0, expression),
      dstmm.zero.expression = tmm_method_normalization(imputed.zero.expression, method, tissue.method, lims_id),
      
      # Gene pareto. "Imputed" values are set to NA to not count.
      gene_dstmm.zero.impute.expression = pareto_scale_method_gene(ifelse(imputed, NA, dstmm.zero.expression), 
                                                                   method, lims_id),
      
      # Limma
      limma_gene_dstmm.zero.impute.expression = limma_method_correction(gene_dstmm.zero.impute.expression, method,
                                                                        tissue.method, lims_id,
                                                                        filtered.methods = "Blood"))  %>%
    # Scale so that under limit is 1, scale by 3
    mutate_at(.funs = funs(. / (under_limit(., expression, method, scale_by = 3))), 
              .vars = grep(".expression$", colnames(.), value = T)) 
  
  
  readr::write_delim(all.atlas, path = paste(result_folder, paste0('all.atlas.txt'),sep='/'), delim = "\t")
} else {
  all.atlas <- readr::read_delim(paste(result_folder, paste0('all.atlas.txt'),sep='/'), delim = "\t")
  }




#
# ----------- Step 3. Consensus ----------- 
#

if(!file.exists(paste(result_folder, paste0('all.atlas.max.txt'),sep='/'))) {
  all.atlas.max <-
    all.atlas %>%
    # Remove genes that are imputed
    filter(!imputed) %>%
    group_by(consensus_content_name, ensg_id) %>% 
    mutate(method = as.character(method)) %>%
    dplyr::summarise_at(.funs = funs(maxEx = max(., na.rm = T),
                                     method = get_method(method, ., max(., na.rm = T))),
                        .vars = grep("expression$", colnames(.), value = T)) 
  readr::write_delim(all.atlas.max, path = paste(result_folder, paste0('all.atlas.max.txt'),sep='/'), delim = "\t")
} else {
  all.atlas.max <- readr::read_delim(paste(result_folder, paste0('all.atlas.max.txt'),sep='/'), delim = "\t")
}


 

#
# ----------- Step 4. Category ----------- 
#

if(!file.exists(paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'))) {
  all.atlas.category <- get.categories.with.num.expressed(all.atlas.max,
                                                          max_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                                          cat_column = "consensus_content_name",
                                                          enrich.fold = 4,
                                                          under.lim = 1,
                                                          group.num = 6)
  readr::write_delim(all.atlas.category, path = paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'), delim = "\t")
} else {
  all.atlas.category <- readr::read_delim(paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'), delim = "\t")
  } 


# print number of different categories of genes
table(all.atlas.category$category.text)
all.atlas.category %>%
  group_by(elevated.category, express.category.2) %>% 
  summarise(n = length(elevated.category)) %>% 
  spread(key = express.category.2, value = n) %>%
  readr::write_delim(paste(result_folder, paste0('gene_categories_summarized.txt'),sep='/'), delim = "\t")

#
# ----------- Blood atlas classification (18 cells) ----------- 
#

blood.atlas <- 
  all.atlas %>%
  filter(method == "Blood") %>% 
  # Remove Total for classification
  filter(content_name != "total PBMC") 

if(!file.exists(paste(result_folder, paste0('blood.atlas.max.txt'),sep='/'))) {
  # blood.atlas.max <- 
  #   blood.atlas %>%
  #   group_by(content_name, ensg_id) %>% 
  #   filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
  #   dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 
  blood.atlas.max <-
    blood.atlas %>%
    # Remove genes that are imputed
    filter(!imputed) %>%
    group_by(content_name, ensg_id) %>% 
    mutate(method = as.character(method)) %>%
    dplyr::summarise_at(.funs = funs(maxEx = max(., na.rm = T),
                                     method = get_method(method, ., max(., na.rm = T))),
                        .vars = grep("expression$", colnames(.), value = T)) 
  readr::write_delim(blood.atlas.max, path = paste(result_folder, paste0('blood.atlas.max.txt'),sep='/'), delim = "\t")
} else {
  blood.atlas.max <- readr::read_delim(paste(result_folder, paste0('blood.atlas.max.txt'),sep='/'), delim = "\t")
} 


if(!file.exists(paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'))) {
  blood.atlas.category <- get.categories.with.num.expressed(blood.atlas.max,
                                                            max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                                                            cat_column = "content_name",
                                                            enrich.fold = 4, 
                                                            under.lim = 1, 
                                                            group.num = 11)
  
  readr::write_delim(blood.atlas.category, path = paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'), delim = "\t")
  
  for(celltype in unique(blood.atlas.max$content_name)) {
    blood.atlas.category %>%
      filter(mapply(paste0("(^|, )", celltype, "(, |$)"), 
                    `enriched tissues`, FUN = function(x,y) grepl(x, y))) %>%
      write_csv(paste(result_folder, paste0(celltype, "_elevated_genes.csv"),sep='/'))
  }
} else {
  blood.atlas.category <- readr::read_delim(paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'), delim = "\t")
} 



# print number of different categories of genes
table(blood.atlas.category$category.text)

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
# ----------- Blood atlas classification (6 cells) ----------- 
#

blood.atlas.6 <- 
  all.atlas %>%
  filter(method == "Blood") %>% 
  left_join(blood_atlas_hierarchy, by = c("content_name" = "content")) %>%
  mutate(content_name = content_l1) %>%
  # Remove Total for classification
  filter(content_name != "total PBMC") 

if(!file.exists(paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'))) {
  blood.atlas.max.6 <- 
    blood.atlas.6 %>%
    group_by(content_name, ensg_id) %>% 
    filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
    dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 
  
  readr::write_delim(blood.atlas.max.6, path = paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'), delim = "\t")
} else {
  blood.atlas.max.6 <- readr::read_delim(paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'), delim = "\t")
} 


if(!file.exists(paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'))) {
  blood.atlas.category.6 <- get.categories.with.num.expressed(blood.atlas.max.6,
                                                              max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                                                              cat_column = "content_name",
                                                              enrich.fold = 4, 
                                                              under.lim = 1, 
                                                              group.num = 5)
  
  readr::write_delim(blood.atlas.category.6, path = paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'), delim = "\t")
} else {
  blood.atlas.category.6 <- readr::read_delim(paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'), delim = "\t")
} 



# print number of different categories of genes
table(blood.atlas.category.6$category.text)

#
# ----------- Brain atlas classification ----------- 
#

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
  filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
  dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 

brain.atlas.max_all_regions <- 
  brain.atlas.filter %>% 
  group_by(content_name, ensg_id) %>% 
  filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
  dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 

write.table(brain.atlas.max,
            file=paste(result_folder, paste0('consensus_brain_regions.txt'),sep='/'),
            row.names = F,
            sep='\t',
            quote=F)

brain.atlas.category <- get.categories.with.num.expressed(brain.atlas.max,
                                                          max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                                                          cat_column = "subgroup",
                                                          enrich.fold = 4, 
                                                          group.num = 6)

brain.atlas.category_all_regions <-
  get.categories.with.num.expressed(brain.atlas.max_all_regions,
                                    max_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                    cat_column = "content_name",
                                    enrich.fold = 4,
                                    group.num = 6)
write.table(brain.atlas.category,
            file=paste(result_folder, paste0('gene_categories_brain_regions.txt'),sep='/'),
            row.names = F,
            sep='\t',
            quote=F)

# print number of different categories of genes
table(brain.atlas.category$category.text)


#
# ----------- Visulization ----------- 
#


plot.data <- 
  all.atlas.max %>%
  select(ensg_id, 
         consensus_content_name, 
         norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
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
                 consensus_content_name = content_name, 
                 norm_exp = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
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


# =========== *All altas =========== 



# =========== *Brain altas* =========== 



# =========== *Blood altas* =========== 
make_plots(atlas = blood.atlas, 
           atlas.max = blood.atlas.max, 
           atlas.cat = blood.atlas.category, 
           Ex_column = "limma_gene_dstmm.zero.impute.expression", 
           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",  
           content_column = "content_name", 
           content_hierarchy = blood_atlas_hierarchy, 
           content_colors = with(blood_atlas_colors, setNames(color, content)), 
           plots = "all", 
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

# # =========== *All altas =========== 
# 
# make_sum_TPM_plot(all.atlas, all.atlas.category, tissue_column = "content_name", method_column = "method", outpath = result_folder, prefix = "atlas")
# 
# make_classification_pie_chart(atlas.cat = all.atlas.category, 
#                               outpath = result_folder, 
#                               prefix = "all_atlas")
# 
# all.atlas.max.wide <- generate_wide(all.atlas.max, ensg_column='ensg_id', 
#                                     group_column='consensus_content_name', 
#                                     max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
# 
# 
# # Tissue distribution
# tissues_to_plot <- c('tongue', 'thalamus', 'skeletal muscle', 'thyroid gland', 'appendix', 'esophagus', 
#                      'postcentral gyrus', 'liver', 'Myeloid DCs', 'olfactory region', 'Memory CD8 T-cells', 
#                      'pancreas', 'hypothalamus', 'Memory B-cells', 'spleen', 'duodenum', 'rectum', 'cerebellum', 
#                      'prostate', 'Naive B-cells', 'vagina', 'endometrium 1', 'heart muscle', 'adrenal gland', 
#                      'skin 1', 'salivary gland', 'thymus', 'small intestine', 'frontal lobe', 'kidney', 
#                      'gallbladder', 'putamen', 'cervix, uterine', 'Eosinophils', 'seminal vesicle', 'pons', 
#                      'placenta', 'ductus deferens', 'amygdala')
# 
# make_tissue_distributions_plot(atlas.tb = all.atlas, 
#                                Ex_column = "limma_gene_dstmm.zero.impute.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = tissues_to_plot, 
#                                outpath = result_folder, 
#                                prefix = "All atlas NX")
# 
# make_tissue_distributions_plot(atlas.tb = all.atlas, 
#                                Ex_column = "expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = tissues_to_plot, 
#                                outpath = result_folder, 
#                                prefix = "X")
# 
# make_tissue_distributions_plot(atlas.tb = all.atlas, 
#                                Ex_column = "dstmm.zero.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = tissues_to_plot, 
#                                outpath = result_folder, 
#                                prefix = "All atlas TMM")
# 
# make_tissue_distributions_plot(atlas.tb = all.atlas, 
#                                Ex_column = "gene_dstmm.zero.impute.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = tissues_to_plot, 
#                                outpath = result_folder, 
#                                prefix = "All atlas TMM pareto")
# 
# ## Spearman method cluster
# make_spearman_method_dendrogram(all.atlas.tb = all.atlas, 
#                                 Ex_column = "limma_gene_dstmm.zero.impute.expression", 
#                                 content_column = "content_name", 
#                                 named_color_replacement = dataset.colors, 
#                                 outpath = result_folder, 
#                                 prefix = "All_atlas_norm_method_color")
# 
# make_spearman_method_dendrogram(all.atlas.tb = all.atlas, 
#                                 Ex_column = "limma_gene_dstmm.zero.impute.expression", 
#                                 content_column = "content_name", 
#                                 named_color_replacement = tissue.colors, 
#                                 outpath = result_folder, 
#                                 prefix = "All_atlas_norm_tissue_color")
# 
# make_spearman_method_dendrogram(all.atlas.tb = all.atlas, 
#                                 Ex_column = "expression", 
#                                 content_column = "content_name", 
#                                 named_color_replacement = dataset.colors, 
#                                 outpath = result_folder, 
#                                 prefix = "All_atlas_exp_method_color")
# 
# make_spearman_method_dendrogram(all.atlas.tb = all.atlas, 
#                                 Ex_column = "expression", 
#                                 content_column = "content_name", 
#                                 named_color_replacement = tissue.colors, 
#                                 outpath = result_folder, 
#                                 prefix = "All_atlas_exp_tissue_color")
# 
# ## tissue distribution of normalized values
# make_tissue_distribution_plot(tb.atlas = all.atlas, 
#                               expr_column = "limma_gene_dstmm.zero.impute.expression",
#                               outpath = result_folder,
#                               prefix = 'all_tissues')
# 
# ## PCA and clustering plots
# all.atlas.max.pca.values <- pca.cal(all.atlas.max.wide)
# scores <- all.atlas.max.pca.values[[1]]
# loadings <- 
#   all.atlas.max.pca.values[[2]] %>%
#   as.tibble(rownames = "ensg_id") %>%
#   mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
# 
# make_PCA_plots(scores = scores,
#                loadings = loadings,
#                groups = setNames(rownames(all.atlas.max.pca.values[[1]]), rownames(all.atlas.max.pca.values[[1]])),
#                groups.color = tissue.colors,
#                outpath = result_folder,
#                prefix = 'all_tissues')
# 
# make_clustering_plot(tb.wide = all.atlas.max.wide, 
#                      colors = tissue.colors, 
#                      outpath = result_folder,
#                      prefix = 'all_tissues')
# 
# ## tissue elevated plot
# all.atlas.elevated.table <- calc_elevated.table(tb.wide = all.atlas.max.wide, 
#                                                 atlas.categories = all.atlas.category)
# all.atlas.elevated.summary.table <- calc_elevated.summary.table(all.atlas.elevated.table)
# make_elevated_bar_plot(elevated.summary.table = all.atlas.elevated.summary.table, 
#                        outpath = result_folder,
#                        prefix = 'all_tissues')
# 
# ## specificity distribution
# make_specificity_distribution_plot(atlas.cat = all.atlas.category, 
#                                    type = "Tissue",
#                                    outpath = result_folder,
#                                    prefix = 'all_tissues')
# 
# ## chord plot
# make_classification_chord_plot(atlas.cat = all.atlas.category,
#                                outpath = result_folder,
#                                prefix = 'all_tissues')
# 
# 
# 
# ## swarm plot
# make_swarm_expression_plot(atlas.max = all.atlas.max, 
#                            atlas.cat = all.atlas.category, 
#                            maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                            tissue_column = "consensus_content_name",
#                            outpath = result_folder,
#                            prefix = 'all_tissues')
# 
# # group enriched chord diagram
# all_atlas_hierarchy <- 
#   contenthierarchy.table.tissue %>%
#   select(1:2) %>%
#   rename(content = 1, content_l1 = 2)
# make_chord_group_enriched(all.atlas.elevated.table, 
#                           grid.col = tissue.colors, 
#                           #tissue_hierarcy = all_atlas_hierarchy,
#                           tissue_hierarcy = rbind(all_atlas_hierarchy, mutate(all_atlas_hierarchy, content = paste(content, 1))),
#                           palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
#                           outpath = result_folder, 
#                           prefix = "all_atlas")
# 
# make_heatmap_group_enriched(all.atlas.elevated.table, 
#                             outpath = result_folder,
#                             prefix = "all_atlas")
# 
# make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
#                                                      all.atlas.max.tb = all.atlas.max, 
#                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                                      tissue_column = "consensus_content_name",
#                                                      outpath = result_folder,
#                                                      prefix = "all_atlas") 
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = all.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          prefix = "all atlas")
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = all.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          range_scale_x = T,
#                          prefix = "all atlas range scaled")
# 
# 
# 
# # make_heatmap_all_elevated_expression_levels_circle(elevated.table = all.atlas.elevated.table,
# #                                                    all.atlas.max.tb = all.atlas.max, 
# #                                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
# #                                                    tissue_column = "consensus_content_name",
# #                                                    outpath = result_folder,
# #                                                    prefix = "all_atlas") 
# 
# # make_heatmap_group_and_enhanced_expression_levels_circle(elevated.table = all.atlas.elevated.table,
# #                                                          all.atlas.max.tb = all.atlas.max, 
# #                                                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
# #                                                          tissue_column = "consensus_content_name",
# #                                                          outpath = result_folder,
# #                                                          prefix = "all_atlas") 
# 
# # make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
# #                                                      all.atlas.max.tb = all.atlas.max, 
# #                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
# #                                                      tissue_column = "consensus_content_name",
# #                                                      outpath = result_folder,
# #                                                      prefix = "all_atlas",
# #                                                      y_dendrogram = F) 
# # 
# # make_heatmap_group_enriched_expression_levels_circle(elevated.table = all.atlas.elevated.table,
# #                                                      all.atlas.max.tb = all.atlas.max, 
# #                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
# #                                                      tissue_column = "consensus_content_name",
# #                                                      outpath = result_folder,
# #                                                      prefix = "all_atlas_dendro",
# #                                                      y_dendrogram = T) 
# 
# # Number of expressed genes
# make_number_detected_genes_barplot(all.atlas.max.tb = all.atlas.max, 
#                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                    tissue_column = "consensus_content_name",
#                                    outpath = result_folder,
#                                    prefix = "all_atlas")
# 
# # Total elevated expression fraction
# make_elevated_NX_fraction_barplots(atlas.max = all.atlas.max, 
#                                    atlas.cat = all.atlas.category, 
#                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                    tissue_column = "consensus_content_name",
#                                    outpath = result_folder, 
#                                    prefix = "all_atlas")
# 
# # =========== *Brain altas* =========== 
# 
# brain.atlas.max.wide <- generate_wide(brain.atlas.max, ensg_column='ensg_id', group_column='subgroup', 
#                                       max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
# 
# brain.atlas.max.wide_all_regions <- generate_wide(brain.atlas.max_all_regions, ensg_column='ensg_id',
#                                                   group_column='content_name',
#                                                   max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
# 
# make_classification_pie_chart(atlas.cat = brain.atlas.category, 
#                               outpath = result_folder, 
#                               prefix = "brain_atlas")
# 
# 
# ## tissue distribution of normalized values
# make_tissue_distribution_plot(tb.atlas = brain.atlas, 
#                               expr_column = "limma_gene_dstmm.zero.impute.expression",
#                               outpath = result_folder,
#                               prefix = 'brain_regions')
# 
# ## PCA and clustering plots
# brain.atlas.max.pca.values <- pca.cal(brain.atlas.max.wide)
# scores <- brain.atlas.max.pca.values[[1]]
# loadings <- 
#   brain.atlas.max.pca.values[[2]] %>%
#   as.tibble(rownames = "ensg_id") %>%
#   mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
# tissue.colors <- with(brainregions.table, setNames(subgroup.color, subgroup))
# 
# make_PCA_plots(scores = scores,
#                loadings = loadings,
#                groups = setNames(rownames(brain.atlas.max.pca.values[[1]]), rownames(brain.atlas.max.pca.values[[1]])),
#                groups.color = tissue.colors,
#                outpath = result_folder,
#                prefix = 'brain_regions')
# 
# make_clustering_plot(tb.wide = brain.atlas.max.wide, 
#                      colors = tissue.colors, 
#                      outpath = result_folder,
#                      prefix = 'brain_regions')
# 
# cell.colors <- with(brainregions.table, setNames(subgroup.color, tissue.type))
# make_clustering_plot(tb.wide = brain.atlas.max.wide_all_regions, 
#                      colors = cell.colors, 
#                      outpath = result_folder,
#                      prefix = 'brain_all_cells')
# 
# ## tissue elevated plot
# brain.atlas.elevated.table <- calc_elevated.table(tb.wide = brain.atlas.max.wide, 
#                                                   atlas.categories = brain.atlas.category)
# brain.atlas.elevated.summary.table <- calc_elevated.summary.table(brain.atlas.elevated.table)
# make_elevated_bar_plot(elevated.summary.table = brain.atlas.elevated.summary.table, 
#                        outpath = result_folder,
#                        prefix = 'brain_regions')
# 
# ## specificity distribution
# make_specificity_distribution_plot(atlas.cat = brain.atlas.category, 
#                                    type = "Tissue",
#                                    outpath = result_folder,
#                                    prefix = 'brain_regions')
# 
# ## chord plot
# make_classification_chord_plot(atlas.cat = brain.atlas.category,
#                                outpath = result_folder,
#                                prefix = 'brain_tissues')
# 
# ## swarm plot
# make_swarm_expression_plot(atlas.max = brain.atlas.max, 
#                            atlas.cat = brain.atlas.category, 
#                            maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                            tissue_column = "subgroup",
#                            outpath = result_folder,
#                            prefix = 'brain_regions')
# 
# # group enriched chord diagram
# #brain_atlas_hierarchy <- readr::read_delim("ref/brain_atlas_hierarchy.txt", delim = "\t")
# contenthierarchy.table.brain <- contenthierarchy.table %>% filter(type=='brain')
# tissue.colors.brain <- with(contenthierarchy.table.brain, setNames(c(color, color, color), c(tissue_name, organ_name, paste(organ_name, 1))))
# 
# brain.atlas.elevated.table_all_regions <-
#   calc_elevated.table(tb.wide = brain.atlas.max.wide_all_regions,
#                       atlas.categories = brain.atlas.category_all_regions)
# 
# brain.atlas.elevated.table <-
#   calc_elevated.table(tb.wide = brain.atlas.max.wide,
#                       atlas.categories = brain.atlas.category)
# 
# brain_atlas_hierarchy <- 
#   contenthierarchy.table.brain %>%
#   select(content=organ_name)
# 
# make_chord_group_enriched(brain.atlas.elevated.table, 
#                           grid.col = tissue.colors.brain, 
#                           tissue_hierarcy = rbind(brain_atlas_hierarchy, mutate(brain_atlas_hierarchy, content = paste(content, 1))),
#                           palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
#                           outpath = result_folder, 
#                           prefix = "brain_atlas")
# 
# make_heatmap_group_enriched(brain.atlas.elevated.table_all_regions, 
#                             outpath = result_folder,
#                             prefix = "brain_atlas_all_regions")
#    
# make_heatmap_median_expression_levels(elevated.table = brain.atlas.elevated.table,
#                                       all.atlas.max.tb = brain.atlas.max, 
#                                       maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                       tissue_column = "subgroup",
#                                       enrichment = c(3),
#                                       outpath = result_folder,
#                                       prefix = "brain_atlas_group_enriched")
# 
# make_heatmap_median_expression_levels(elevated.table = brain.atlas.elevated.table,
#                                       all.atlas.max.tb = brain.atlas.max, 
#                                       maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                       tissue_column = "subgroup",
#                                       enrichment = c(2, 3, 4),
#                                       outpath = result_folder,
#                                       prefix = "brain_atlas_all_elevated")
# 
# make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
#                                all.atlas.max.tb = brain.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "subgroup",
#                                enrichment = c(3),
#                                outpath = result_folder,
#                                prefix = "brain_atlas_group_enriched")
# 
# make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
#                                all.atlas.max.tb = brain.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "subgroup",
#                                enrichment = c(2,3,4),
#                                outpath = result_folder,
#                                prefix = "brain_atlas_all_elevated")
# 
# make_heatmap_expression_levels(elevated.table = brain.atlas.elevated.table,
#                                all.atlas.max.tb = brain.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "subgroup",
#                                enrichment = c(2),
#                                outpath = result_folder,
#                                prefix = "brain_atlas_tissue_enriched")
# 
# # =========== *Blood altas* =========== 
# 
# blood.atlas.max.wide <- generate_wide(blood.atlas.max, ensg_column='ensg_id', group_column='content_name', 
#                                       max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
# 
# make_classification_pie_chart(atlas.cat = blood.atlas.category, 
#                               outpath = result_folder, 
#                               prefix = "blood_atlas")
# 
# blood.atlas.category %>%
#   left_join(proteinclass.table, by = c("ensg_id" = "rna.genes")) %>%
#   filter(category %in% 2:4) %$%
#   table(proteinclass.vec.single)
# 
# 
# blood.atlas.max %>%
#   filter(ensg_id %in% unique(ensg_id)[1:100]) %>%
#   make_gene_expression_barplot(maxEx_columns = c("Raw" = "expression_maxEx", "TMM" = "dstmm.zero.expression_maxEx", "TMM + Pareto" = "gene_dstmm.zero.impute.expression_maxEx"),
#                                content_column = "content_name", 
#                                content_color = with(blood_atlas_colors, setNames(color, content)))
# 
# # Tissue distribution
# make_tissue_distributions_plot(atlas.tb = blood.atlas, 
#                                Ex_column = "limma_gene_dstmm.zero.impute.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = "all", 
#                                outpath = result_folder, 
#                                prefix = "Blood atlas NX")
# 
# make_tissue_distributions_plot(atlas.tb = blood.atlas, 
#                                Ex_column = "expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = "all", 
#                                outpath = result_folder, 
#                                prefix = "Blood atlas X")
# 
# make_tissue_distributions_plot(atlas.tb = blood.atlas, 
#                                Ex_column = "dstmm.zero.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = "all", 
#                                outpath = result_folder, 
#                                prefix = "Blood atlas TMM")
# 
# make_tissue_distributions_plot(atlas.tb = blood.atlas, 
#                                Ex_column = "gene_dstmm.zero.impute.expression", 
#                                content_column = "content_name",
#                                und.lim = 1, 
#                                do.tissues = "all", 
#                                outpath = result_folder, 
#                                prefix = "Blood atlas TMM pareto")
# 
# ## tissue distribution of normalized values
# make_tissue_distribution_plot(tb.atlas = blood.atlas, 
#                               expr_column = "limma_gene_dstmm.zero.impute.expression",
#                               outpath = result_folder,
#                               prefix = 'blood_cells')
# 
# ## PCA and clustering plots
# blood.atlas.max.pca.values <- pca.cal(blood.atlas.max.wide)
# scores <- blood.atlas.max.pca.values[[1]]
# loadings <- 
#   blood.atlas.max.pca.values[[2]] %>%
#   as.tibble(rownames = "ensg_id") %>%
#   mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
# tissue.colors <- with(contenthierarchy.table, setNames(c(color, color, color), c(tissue_name, organ_name, paste(tissue_name, 1))))
# 
# make_PCA_plots(scores = scores,
#                loadings = loadings,
#                groups = setNames(rownames(blood.atlas.max.pca.values[[1]]), rownames(blood.atlas.max.pca.values[[1]])),
#                groups.color = tissue.colors,
#                outpath = result_folder,
#                prefix = 'blood_celltypes')
# 
# make_clustering_plot(tb.wide = blood.atlas.max.wide, 
#                      colors = tissue.colors, 
#                      outpath = result_folder,
#                      prefix = 'blood_celltypes')
# 
# ## tissue elevated plot
# blood.atlas.elevated.table <- calc_elevated.table(tb.wide = blood.atlas.max.wide, 
#                                                   atlas.categories = blood.atlas.category, 
#                                                   cat.colum = "category")
# blood.atlas.elevated.summary.table <- calc_elevated.summary.table(blood.atlas.elevated.table)
# make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table, 
#                        outpath = result_folder,
#                        translate_categories = c("Celltype" = "Tissue"),
#                        prefix = 'blood_celltypes')
# 
# make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table, 
#                        outpath = result_folder,
#                        translate_categories = c("Celltype" = "Tissue"),
#                        prefix = 'blood_celltypes')
# 
# make_score_expression_scatter(atlas.max.tb = blood.atlas.max, 
#                               atlas.cat = blood.atlas.category, 
#                               maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                               tissue_column = "content_name", 
#                               ensemblanno.table = ensemblanno.table,
#                               plot.order = blood_atlas_hierarchy %>%
#                                 filter(!content %in% c("blood", "Total PBMCs")) %$% 
#                                 content[order(content_l1)],  
#                               outpath = result_folder, 
#                               prefix = "blood_celltypes")
# 
# ## specificity distribution
# make_specificity_distribution_plot(atlas.cat = blood.atlas.category, 
#                                    type = "Tissue",
#                                    outpath = result_folder,
#                                    prefix = 'blood_celltypes')
# 
# ## chord plot
# make_classification_chord_plot(atlas.cat = blood.atlas.category,
#                                outpath = result_folder,
#                                prefix = 'blood_tissues')
# 
# ## swarm plot
# make_swarm_expression_plot(atlas.max = blood.atlas.max, 
#                            atlas.cat = blood.atlas.category, 
#                            maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                            tissue_column = "content_name",
#                            plot.order = blood_atlas_hierarchy %>%
#                              filter(!content %in% c("blood", "Total PBMCs")) %$% 
#                              content[order(content_l1)],
#                            outpath = result_folder,
#                            prefix = 'blood_celltypes')
# 
#   
# make_swarm_expression_circle_plot(atlas.max = blood.atlas.max, 
#                                   atlas.cat = blood.atlas.category, 
#                                   maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                                   tissue_column = "content_name",
#                                   plot.order = blood_atlas_hierarchy %>%
#                                     filter(!content %in% c("blood", "Total PBMCs")) %$% 
#                                     content[order(content_l1)],
#                                   outpath = result_folder,
#                                   prefix = 'blood_celltypes')
# # group enriched chord diagram
# blood_atlas_hierarchy <- readr::read_delim("ref/blood_atlas_hierarchy.txt", delim = "\t")
# blood_atlas_colors <- readr::read_delim("ref/blood_atlas_colors.txt", delim = "\t")
# 
# make_chord_group_enriched(blood.atlas.elevated.table, 
#                           grid.col = with(blood_atlas_colors, setNames(color, content)), 
#                           tissue_hierarcy = blood_atlas_hierarchy,
#                           palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
#                           outpath = result_folder,reverse = T, 
#                           prefix = "blood_atlas")
# 
# make_heatmap_group_enriched(elevated.table = blood.atlas.elevated.table, 
#                             outpath = result_folder,
#                             prefix = "blood_atlas")
# 
# make_expression_heatmaps(atlas.max.tb = blood.atlas.max, 
#                          atlas.cat = blood.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          prefix = "blood atlas")
# 
# make_expression_heatmaps(atlas.max.tb = blood.atlas.max, 
#                          atlas.cat = blood.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          range_scale_x = T, 
#                          prefix = "blood atlas range scaled")
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = blood.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          prefix = "blood atlas cat on all atlas")
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = blood.atlas.category, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          range_scale_x = T,
#                          prefix = "blood atlas cat on all atlas range scaled")
# 
# 
# make_immunodeficiency_expression_heatmaps(atlas.max.tb = blood.atlas.max, 
#                                           atlas.cat = blood.atlas.category, 
#                                           immunodeficiency.table = immunodeficiency.table,
#                                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                                           tissue_column = "content_name", 
#                                           ensemblanno.table = ensemblanno.table,
#                                           proteinclass.table = proteinclass.table, 
#                                           proteinclass.table_ensg_id_column = "rna.genes", 
#                                           proteinclass.table_class_column = "proteinclass.vec.single", 
#                                           outpath = result_folder, 
#                                           range_scale_x = F,
#                                           prefix = "blood atlas immuno deficiency")
# 
# make_immunodeficiency_expression_heatmaps(atlas.max.tb = blood.atlas.max, 
#                                           atlas.cat = blood.atlas.category, 
#                                           immunodeficiency.table = immunodeficiency.table,
#                                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                                           tissue_column = "content_name", 
#                                           ensemblanno.table = ensemblanno.table,
#                                           proteinclass.table = proteinclass.table, 
#                                           proteinclass.table_ensg_id_column = "rna.genes", 
#                                           proteinclass.table_class_column = "proteinclass.vec.single", 
#                                           outpath = result_folder, 
#                                           range_scale_x = T,
#                                           prefix = "blood atlas immuno deficiency range scaled")
# 
# # make_heatmap_group_enriched_expression_levels_circle(elevated.table = blood.atlas.elevated.table,
# #                                                      all.atlas.max.tb = blood.atlas.max, 
# #                                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
# #                                                      tissue_column = "content_name",
# #                                                      outpath = result_folder,
# #                                                      prefix = "blood_atlas") 
# 
# make_heatmap_median_expression_levels(elevated.table = blood.atlas.elevated.table,
#                                       all.atlas.max.tb = blood.atlas.max, 
#                                       maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                       tissue_column = "content_name",
#                                       enrichment = c(3),
#                                       outpath = result_folder,
#                                       prefix = "blood_atlas_group_enriched")
# 
# make_heatmap_median_expression_levels(elevated.table = blood.atlas.elevated.table,
#                                       all.atlas.max.tb = blood.atlas.max, 
#                                       maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                       tissue_column = "content_name",
#                                       enrichment = c(2, 3, 4),
#                                       outpath = result_folder,
#                                       prefix = "blood_atlas_all_elevated")
# 
# make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
#                                all.atlas.max.tb = blood.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "content_name",
#                                enrichment = c(3),
#                                outpath = result_folder,
#                                prefix = "blood_atlas_group_enriched")
# 
# make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
#                                all.atlas.max.tb = blood.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "content_name",
#                                enrichment = c(2,3,4),
#                                outpath = result_folder,
#                                prefix = "blood_atlas_all_elevated")
# 
# make_heatmap_expression_levels(elevated.table = blood.atlas.elevated.table,
#                                all.atlas.max.tb = blood.atlas.max, 
#                                maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                tissue_column = "content_name",
#                                enrichment = c(2),
#                                outpath = result_folder,
#                                prefix = "blood_atlas_tissue_enriched")
# 
# # Categories between blood and all atlas
# make_class_comparison_chord(cat1 = blood.atlas.category, 
#                             cat2 = all.atlas.category,
#                             outpath = result_folder, prefix = "blood")
# 
# # Number of expressed genes
# make_number_detected_genes_barplot(all.atlas.max.tb = blood.atlas.max, 
#                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                    tissue_column = "content_name",
#                                    outpath = result_folder,
#                                    prefix = "blood_atlas")
# 
# # Comparison of elevated genes to tissue atlas
# make_elevated_organ_total_chord(cat1 = blood.atlas.category, 
#                                 cat2 = all.atlas.category, 
#                                 grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))), 
#                                 elevated_cats = c(2,3,4), 
#                                 direction = 1, 
#                                 cat1_name = "celltypes", 
#                                 cat2_name = "tissues",
#                                 outpath = result_folder, 
#                                 prefix = "Blood to tissue elevated")
# 
# 
# make_elevated_organ_total_chord(cat1 = blood.atlas.category, 
#                                 cat2 = all.atlas.category, 
#                                 grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))), 
#                                 elevated_cats = c(2,3,4), 
#                                 direction = 2, 
#                                 cat1_name = "celltypes", 
#                                 cat2_name = "tissues",
#                                 outpath = result_folder, 
#                                 prefix = "Tissue to Blood elevated")
# 
# make_elevated_organ_total_chord(cat1 = blood.atlas.category, 
#                                 cat2 = all.atlas.category, 
#                                 grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))), 
#                                 elevated_cats = c(2), 
#                                 direction = 1, 
#                                 cat1_name = "celltypes", 
#                                 cat2_name = "tissues",
#                                 outpath = result_folder, 
#                                 prefix = "Blood to tissue tissue enriched")
# 
# make_elevated_organ_total_chord(cat1 = blood.atlas.category, 
#                                 cat2 = all.atlas.category, 
#                                 grid.col = c(tissue.colors, with(blood_atlas_colors, setNames(color, content))), 
#                                 elevated_cats = c(2), 
#                                 direction = 2, 
#                                 cat1_name = "celltypes", 
#                                 cat2_name = "tissues",
#                                 outpath = result_folder, 
#                                 prefix = "Tissue to Blood tissue enriched")
# 
# # Total elevated expression fraction
# make_elevated_NX_fraction_barplots(atlas.max = blood.atlas.max, 
#                                    atlas.cat = blood.atlas.category, 
#                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                    tissue_column = "content_name",
#                                    outpath = result_folder, 
#                                    prefix = "blood_atlas")
# 
# # =========== *Blood altas (6 cells)* =========== 
# 
# blood.atlas.max.wide.6 <- generate_wide(blood.atlas.max.6, ensg_column='ensg_id', group_column='content_name', 
#                                         max_column="limma_gene_dstmm.zero.impute.expression_maxEx")
# 
# make_classification_pie_chart(atlas.cat = blood.atlas.category.6, 
#                               outpath = result_folder, 
#                               prefix = "blood_atlas_6")
# 
# 
# ## PCA and clustering plots
# blood.atlas.max.pca.values.6 <- pca.cal(blood.atlas.max.wide.6)
# scores <- blood.atlas.max.pca.values.6[[1]]
# loadings <- 
#   blood.atlas.max.pca.values.6[[2]] %>%
#   as.tibble(rownames = "ensg_id") %>%
#   mutate(labels = ensemblanno.table$gene_name[match(ensg_id, ensemblanno.table$ensg_id)])
# 
# make_PCA_plots(scores = scores,
#                loadings = loadings,
#                groups = setNames(rownames(blood.atlas.max.pca.values.6[[1]]), rownames(blood.atlas.max.pca.values.6[[1]])),
#                groups.color = with(blood_atlas_colors, setNames(color, content)),
#                outpath = result_folder,
#                prefix = 'blood_celltypes_6')
# 
# make_clustering_plot(tb.wide = blood.atlas.max.wide.6, 
#                      colors = with(blood_atlas_colors, setNames(color, content)), 
#                      outpath = result_folder,
#                      prefix = 'blood_celltypes_6')
# 
# 
# ## tissue elevated plot
# blood.atlas.elevated.table.6 <- calc_elevated.table(tb.wide = blood.atlas.max.wide.6, 
#                                                     atlas.categories = blood.atlas.category.6)
# blood.atlas.elevated.summary.table.6 <- calc_elevated.summary.table(blood.atlas.elevated.table.6)
# make_elevated_bar_plot(elevated.summary.table = blood.atlas.elevated.summary.table.6, 
#                        outpath = result_folder,
#                        prefix = 'blood_celltypes_6')
# 
# ## specificity distribution
# make_specificity_distribution_plot(atlas.cat = blood.atlas.category.6, 
#                                    type = "Tissue",
#                                    outpath = result_folder,
#                                    prefix = 'blood_celltypes_6')
# 
# 
# make_swarm_expression_plot(atlas.max = blood.atlas.max.6, 
#                            atlas.cat = blood.atlas.category.6, 
#                            maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                            tissue_column = "content_name",
#                            plot.order = blood_atlas_hierarchy %>%
#                              filter(!content %in% c("blood", "Total PBMCs")) %$% 
#                              unique(content_l1[order(content_l2)]),
#                            outpath = result_folder,
#                            prefix = 'blood_cell_lineages_6')
# 
# 
# ## chord plot
# # make_classification_chord_plot(atlas.cat = blood.atlas.category.6,
# #                                outpath = result_folder,
# #                                prefix = 'blood_tissues_6')
# 
# # group enriched chord diagram
# 
# 
# make_chord_group_enriched(blood.atlas.elevated.table.6, 
#                           grid.col = with(blood_atlas_colors, setNames(color, content)), 
#                           tissue_hierarcy = blood_atlas_hierarchy %>%
#                             select(content_l1, content_l2, content_l3) %>%
#                             rename(content = content_l1,
#                                    content_l1 = content_l2,
#                                    content_l2 = content_l3),
#                           palet = colorRampPalette(colors = c("yellow", "orangered", "#800026")),
#                           outpath = result_folder,reverse = T, 
#                           prefix = "blood_atlas_6")
# 
# make_heatmap_group_enriched(elevated.table = blood.atlas.elevated.table.6, 
#                             outpath = result_folder,
#                             prefix = "blood_atlas_6")
# 
# 
# make_expression_heatmaps(atlas.max.tb = blood.atlas.max.6, 
#                          atlas.cat = blood.atlas.category.6, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          prefix = "blood atlas 6")
# 
# make_expression_heatmaps(atlas.max.tb = blood.atlas.max.6, 
#                          atlas.cat = blood.atlas.category.6, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          range_scale_x = T, 
#                          prefix = "blood atlas 6 range scaled")
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = blood.atlas.category.6, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          prefix = "blood atlas 6 cat on all atlas")
# 
# make_expression_heatmaps(atlas.max.tb = all.atlas.max, 
#                          atlas.cat = blood.atlas.category.6, 
#                          maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                          tissue_column = "consensus_content_name", 
#                          ensemblanno.table = ensemblanno.table,
#                          proteinclass.table = proteinclass.table, 
#                          proteinclass.table_ensg_id_column = "rna.genes", 
#                          proteinclass.table_class_column = "proteinclass.vec.single", 
#                          outpath = result_folder, 
#                          range_scale_x = T,
#                          prefix = "blood atlas 6 cat on all atlas range scaled")
# 
# # Total elevated expression fraction
# make_elevated_NX_fraction_barplots(atlas.max = blood.atlas.max.6, 
#                                    atlas.cat = blood.atlas.category.6, 
#                                    maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
#                                    tissue_column = "content_name",
#                                    outpath = result_folder, 
#                                    prefix = "blood_atlas_6")

