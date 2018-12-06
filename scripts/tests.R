
###################################################################
## Test to see bias in imputed values
###################################################################

ensg_id_imputation_status <- 
  all.atlas %>%
  select(ensg_id, imputed, method) %>%
  unique() %>%
  group_by(ensg_id) %>%
  summarise(imputed_for = paste(sort(method[which(imputed)]), collapse = " "), 
            imputed = any(imputed))

all.atlas.max_imputation_status <- 
  all.atlas.max %>%
  select("consensus_content_name", "ensg_id", "expression_maxEx", 
         "imputed.expression_maxEx", "limma_gene_dstmm.expression_maxEx", "limma_gene_dstmm.expression_method") %>%
  left_join(ensg_id_imputation_status, by = "ensg_id")

all.atlas.max_imputation_status %$%
  tibble(X = log10(expression_maxEx + 1), 
         NX = log10(limma_gene_dstmm.expression_maxEx + 1), 
         imputed = imputed, 
         imputed_for = imputed_for,
         method = limma_gene_dstmm.expression_method) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues")

ggsave(paste(result_folder, "BA NX vs X all tissues.png", sep = "/"), width = 8, height = 8)

all.atlas.max_imputation_status %>%
  group_by(ensg_id) %>%
  summarise(expression_maxEx = mean(expression_maxEx, na.rm = T),
            limma_gene_dstmm.expression_maxEx = mean(limma_gene_dstmm.expression_maxEx, na.rm = T),
            imputed = any(imputed),
            imputed_for = "") %$%
  tibble(X = log10(expression_maxEx + 1), 
         NX = log10(limma_gene_dstmm.expression_maxEx + 1), 
         imputed = imputed, 
         imputed_for = imputed_for) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues mean")

ggsave(paste(result_folder, "BA NX vs X all tissues mean.png", sep = "/"), width = 8, height = 8)

for(tissue in consensustissue.table$consensus_content_name) {
  all.atlas.max_imputation_status %>%
    filter(consensus_content_name == tissue) %$%
    tibble(X = log10(expression_maxEx + 1), 
           NX = log10(limma_gene_dstmm.expression_maxEx + 1), 
           imputed = imputed, 
           imputed_for = imputed_for,
           method = limma_gene_dstmm.expression_method) %>%
    mutate(x = (NX + X) / 2,
           y = (NX - X)) %$%
    make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = tissue) %>% 
    print()
  ggsave(paste(result_folder, paste0("BA NX vs X ", tissue, ".png"), sep = "/"), width = 8, height = 8)
}

# HPA
all.atlas %>%
  left_join(select(ensg_id_imputation_status, -imputed), by = "ensg_id") %>%
  filter(method == "HPA") %$%
  tibble(X = log10(imputed.expression + 1), 
         NX = log10(limma_gene_dstmm.expression + 1), 
         imputed = imputed, 
         imputed_for = imputed_for,
         method = method) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "HPA") %>% 
  print()
ggsave(paste(result_folder, "BA NX vs X HPA.png", sep = "/"), width = 8, height = 8)

# FANTOM
all.atlas %>%
  left_join(select(ensg_id_imputation_status, -imputed), by = "ensg_id") %>%
  filter(method == "FANTOM") %$%
  tibble(X = log10(imputed.expression + 1), 
         NX = log10(limma_gene_dstmm.expression + 1), 
         imputed = imputed, 
         imputed_for = imputed_for,
         method = method) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "FANTOM") %>% 
  print()
ggsave(paste(result_folder, "BA NX vs X FANTOM.png", sep = "/"), width = 8, height = 8)

# Blood
all.atlas %>%
  left_join(select(ensg_id_imputation_status, -imputed), by = "ensg_id") %>%
  filter(method == "Blood") %$%
  tibble(X = log10(imputed.expression + 1), 
         NX = log10(limma_gene_dstmm.expression + 1), 
         imputed = imputed, 
         imputed_for = imputed_for,
         method = method) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "Blood") %>% 
  print()
ggsave(paste(result_folder, "BA NX vs X Blood.png", sep = "/"), width = 8, height = 8)

# FANTOM
all.atlas %>%
  left_join(select(ensg_id_imputation_status, -imputed), by = "ensg_id") %>%
  filter(method == "FANTOM") %$%
  tibble(X = log10(imputed.expression + 1), 
         NX = log10(limma_gene_dstmm.expression + 1), 
         imputed = imputed, 
         imputed_for = imputed_for,
         method = method) %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "FANTOM") %>% 
  print()
ggsave(paste(result_folder, "BA NX vs X FANTOM.png", sep = "/"), width = 8, height = 8)
