if(all(file.exists(c(hpa_norm_path, 
                     gtex_norm_path, 
                     fantom_norm_path, 
                     blood_norm_path)))) {
  
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
  
} else if(!file.exists(paste(result_folder, paste0('all.atlas.txt'),sep='/'))) {
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







# if(!file.exists(paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'))) {
#   blood.atlas.max.6 <- 
#     blood.atlas.6 %>%
#     group_by(content_name, ensg_id) %>% 
#     filter(!is.na(limma_gene_dstmm.zero.impute.expression)) %>%
#     dplyr::summarise(limma_gene_dstmm.zero.impute.expression_maxEx = max(limma_gene_dstmm.zero.impute.expression)) 
#   
#   readr::write_delim(blood.atlas.max.6, path = paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'), delim = "\t")
# } else {
#   blood.atlas.max.6 <- readr::read_delim(paste(result_folder, paste0('blood.atlas.max.6.txt'),sep='/'), delim = "\t")
# } 
# 
# 
# if(!file.exists(paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'))) {
#   blood.atlas.category.6 <- get.categories.with.num.expressed(blood.atlas.max.6,
#                                                               max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
#                                                               cat_column = "content_name",
#                                                               enrich.fold = 4, 
#                                                               under.lim = 1, 
#                                                               group.num = 5)
#   
#   readr::write_delim(blood.atlas.category.6, path = paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'), delim = "\t")
# } else {
#   blood.atlas.category.6 <- readr::read_delim(paste(result_folder, paste0('gene_categories_blood_cells.6.txt'),sep='/'), delim = "\t")
# } 
# 
# 
#


lims_files <- F

if(file.exists(consensus_path)) {
  
  lims_files <- T
  
  all.atlas.max <- 
    read_delim(consensus_path, delim = "\t") %>%
    rename(consensus_content_name = 2)
  
  # Blood Atlas 
  all.atlas.max <- 
    read_delim(consensus_blood_path, delim = "\t") %>%
    rename(consensus_content_name = 2)
  
  
} else if(!file.exists(paste(result_folder, paste0('all.atlas.max.txt'),sep='/'))) {
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
  all.atlas.max <- readr::read_delim(paste(result_folder, paste0('all.atlas.max.txt'),sep='/'), delim = "\t")
  blood.atlas.max <- readr::read_delim(paste(result_folder, paste0('blood.atlas.max.txt'),sep='/'), delim = "\t")
}










# if(file.exists(category_path)) {
#   
#   all.atlas.category_asasas <- 
#     read_delim(category_path, delim = "\t") 
#   
#   all.atlas.category %>% 
#     full_join(all.atlas.category_asasas, by = c("ensg_id")) %>% 
#     select("enriched tissues", "tissue/group specific score", "express.category.2", "elevated.category",
#            "specificity_category", "distribution_category", "ts_score", "enhanced_score", "enhanced_tissues") %>% 
#     select(1, 9)
#   
# } else 

if(lims_files) {
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
  
  
} else if(!file.exists(paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'))) {
  all.atlas.category <- get.categories.with.num.expressed(all.atlas.max,
                                                          max_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                                          cat_column = "consensus_content_name",
                                                          enrich.fold = 4,
                                                          under.lim = 1,
                                                          group.num = 6)
  readr::write_delim(all.atlas.category, path = paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'), delim = "\t")
  blood.atlas.category <- get.categories.with.num.expressed(blood.atlas.max,
                                                            max_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                                                            cat_column = "content_name",
                                                            enrich.fold = 4, 
                                                            under.lim = 1, 
                                                            group.num = 11)
  readr::write_delim(blood.atlas.category, path = paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'), delim = "\t")
  
} else {
  all.atlas.category <- readr::read_delim(paste(result_folder, paste0('gene_categories_all_tissues.txt'),sep='/'), delim = "\t")
  blood.atlas.category <- readr::read_delim(paste(result_folder, paste0('gene_categories_blood_cells.txt'),sep='/'), delim = "\t")
} 



# if(file.exists(category_path)) {
#   
#   all.atlas.category_asasas <- 
#     read_delim(category_path, delim = "\t") 
#   
#   all.atlas.category %>% 
#     full_join(all.atlas.category_asasas, by = c("ensg_id")) %>% 
#     select("enriched tissues", "tissue/group specific score", "express.category.2", "elevated.category",
#            "specificity_category", "distribution_category", "ts_score", "enhanced_score", "enhanced_tissues") %>% 
#     select(1, 9)
#   
# } else 

#######################################

all.atlas.raw %>% 
  group_by(method, content_name) %>%
  summarise(n = sum(expression))

all.atlas.raw %>%
  group_by(method, content_name) %>%
  mutate(expression = 1e6 * expression/sum(expression)) %>%
  summarise(n = sum(expression))

all.atlas.1 <- 
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
                                                                 method, lims_id))

all.atlas.2 <- 
  all.atlas.1 %>% 
  mutate(
    # Limma
    limma_gene_dstmm.zero.impute.expression = limma_method_correction(gene_dstmm.zero.impute.expression, method,
                                                                      tissue.method, lims_id,
                                                                      filtered.methods = "Blood"))  


all.atlas.2 %>% View
all.atlas.2 %>% 
  {.[sample(1:nrow(.), 100000), ]} %>%
  ggplot(aes(gene_dstmm.zero.impute.expression, limma_gene_dstmm.zero.impute.expression, color = method))+
  geom_point()

all.atlas.2 %>% 
  filter(method == "Blood") %$% 
  all(gene_dstmm.zero.impute.expression == limma_gene_dstmm.zero.impute.expression)
# LIMMA: Batch correction of method, method and gene pareto scaled values without blood
A <- tibble(expression, tissue.method, ensg_id)

B <- 
  filter(A, !method %in% filtered.methods) %>%
  spread(key = tissue.method, value = expression) %>%
  column_to_rownames("ensg_id") %>%
  as.matrix() %>%
  {.[apply(., 1, function(x) !any(is.na(x))),]} %>%
  {log(.+1, 10)} %>%
  limma::removeBatchEffect(batch = as.factor(str_extract(colnames(.), "(?<=\\.).*")))
  
C <- 
  B %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(ensg_id = names)} %>%
  gather(key = "tissue.method", value = "method.corrected.expression", -ensg_id) %>%
  mutate(method.corrected.expression = 10^method.corrected.expression - 1,
         method.corrected.expression = ifelse(is.na(method.corrected.expression),
                                              norm.dsscaled.expression, 
                                              method.corrected.expression),
         method.corrected.expression = ifelse(method.corrected.expression<0, 0,
                                              method.corrected.expression))

D <- 
  left_join(A, C, 
            by = c("ensg_id", "tissue.method")) 

D %>% filter(method %in% filtered.methods)

E <- 
  D %>% 
  mutate(method.corrected.expression = ifelse(method %in% filtered.methods, 
                                              expression, 
                                              method.corrected.expression)) 

E %>% filter(method %in% filtered.methods)
  