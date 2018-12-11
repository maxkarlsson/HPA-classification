
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

table(ensg_id_imputation_status$imputed_for)

all.atlas.max_imputation_status <- 
  all.atlas.max %>%
  select("consensus_content_name", "ensg_id", 
         "expression_maxEx", "imputed.expression_maxEx", 
         "limma_gene_dstmm.expression_maxEx", "limma_gene_dstmm.expression_method",
         
         #"imputed.log.expression_maxEx", "limma_gene_dstmm.log.expression_maxEx",
         #"limma_gene_dstmm.zero.expression_maxEx", 
         "limma_gene_dstmm.zero.impute.expression_maxEx") %>%
  left_join(ensg_id_imputation_status, by = "ensg_id") %>%
  mutate(X = log10(expression_maxEx + 1), 
         NX = log10(limma_gene_dstmm.expression_maxEx + 1), 
         IX = log10(imputed.expression_maxEx + 1), 
         
         # Log
         #ILX = log10(imputed.log.expression_maxEx + 1),
         #NLX = log10(limma_gene_dstmm.log.expression_maxEx + 1), 
         
         # NA -> Zero
         #NZX = log10(limma_gene_dstmm.zero.expression_maxEx + 1), 
         
         # NA -> Zero -> TMM -> LIMMA
         NZTX = log10(limma_gene_dstmm.zero.impute.expression_maxEx + 1))

# imputed.zero.expression = ifelse(imputed, 0, expression),
# dstmm.zero.expression = tmm_method_normalization(imputed.zero.expression, method, tissue.method, lims_id),
# gene_dstmm.zero.expression = pareto_scale_method_gene(dstmm.zero.expression, method, lims_id),
# 
# limma_gene_dstmm.zero.expression = limma_method_correction(gene_dstmm.zero.expression, method,
#                                                            tissue.method, lims_id,
#                                                            filtered.methods = "Blood"),
# # 3. zero impute then impute again after tmm
# 
# gene_dstmm.zero.impute.expression = pareto_scale_method_gene(ifelse(imputed, NA, dstmm.zero.expression), method, lims_id),
# 
# limma_gene_dstmm.zero.impute.expression = limma_method_correction(gene_dstmm.zero.impute.expression, method,
#                                                                   tissue.method, lims_id,
#                                                                   filtered.methods = "Blood")



all.atlas %>% 
    left_join(select(ensg_id_imputation_status, -imputed), by = "ensg_id") %>%
    filter(imputed_for != "") %>%
  gather(key = "normalisation", value = "value", expression, imputed.expression, limma_gene_dstmm.expression, limma_gene_dstmm.zero.impute.expression) %>%
    ggplot(aes(value + 1, fill = method)) +
    geom_histogram(alpha= 0.5)+
    geom_vline(xintercept = 2) + 
    simple_theme + 
  facet_wrap( ~ normalisation, ncol = 4) +
    scale_x_log10()+
    scale_y_continuous()

ggsave(paste(result_folder, "distribution tissues.png", sep = "/"), width = 15, height = 8)

  




all.atlas.max_imputation_status  %>%
  mutate(x = (NZTX + NX) / 2,
         y = (NZTX - NX)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues")

all.atlas.max_imputation_status  %>%
  mutate(x = (NZTX + NX) / 2,
         y = (NZTX - NX)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues", Points = T)


all.atlas.max_imputation_status  %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues")

all.atlas.max_imputation_status  %>%
  mutate(x = (ILX + IX) / 2,
         y = (ILX - IX)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues", 
                         xl = "Mean of ILX and IX", yl = "ILX - IX")

all.atlas.max_imputation_status  %>%
  mutate(x = (NLX + X) / 2,
         y = (NLX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues", 
                         xl = "Mean of NLX and X", yl = "NLX - X")

all.atlas.max_imputation_status  %>%
  mutate(x = (NZX + X) / 2,
         y = (NZX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues", 
                         xl = "Mean of NZX and X", yl = "NZX - X")

all.atlas.max_imputation_status  %>%
  mutate(x = (NZTX + X) / 2,
         y = (NZTX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues", 
                         xl = "Mean of NZTX and X", yl = "NZTX - X")

all.atlas.max_imputation_status_mean <- 
  all.atlas.max_imputation_status %>%
  group_by(ensg_id) %>%
  summarise(X = mean(X, na.rm = T),
            NX = mean(NX, na.rm = T),
            IX = mean(IX, na.rm = T),
            NLX = mean(NLX, na.rm = T),
            ILX = mean(ILX, na.rm = T),
            NZX = mean(NZX, na.rm = T),
            NZTX = mean(NZTX, na.rm = T),
            imputed = any(imputed),
            imputed_for = "") 

grid.arrange(all.atlas.max_imputation_status_mean %>%
  mutate(x = (NX + X) / 2,
         y = (NX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues mean"),

all.atlas.max_imputation_status_mean %>%
  mutate(x = (NLX + X) / 2,
         y = (NLX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues mean"),
all.atlas.max_imputation_status_mean %>%
  mutate(x = (NZX + X) / 2,
         y = (NZX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues mean"),
all.atlas.max_imputation_status_mean %>%
  mutate(x = (NZTX + X) / 2,
         y = (NZTX - X)) %$%
  make_bland_altman_plot(x, y, imputed, fillname = "imputed", title = "All tissues mean"), ncol = 2)




# all.atlas %>% filter(imputed) %>% left_join(all.atlas.max_imputation_status, by = c("ensg_id", "consensus_content_name"))->A
# ggplot(A, aes(imputed.expression+1, imputed.expression_maxEx+1, color = method))+geom_point() + scale_x_log10() + scale_y_log10()








ggsave(paste(result_folder, "BA NX vs X all tissues.png", sep = "/"), width = 8, height = 8)



ggsave(paste(result_folder, "BA NX vs X all tissues mean.png", sep = "/"), width = 8, height = 8)

##
# Imputed vs expression
all.atlas.max_imputation_status %>% 
  filter(imputed) %$%
  {grid.arrange(mutate(.,
                       x = (NX + X) / 2,
                       y = (NX - X)) %$%
                  make_bland_altman_plot(x, y, imputed_for, fillname = "imputed", 
                                         title = "All tissues only imputed", Points = T, alpha = 0.05,
                                         xl = "Mean of NX and X", yl = "NX - X"),
                mutate(.,
                       x = (NZX + X) / 2,
                       y = (NZX - X)) %$%
                  make_bland_altman_plot(x, y, imputed_for, fillname = "imputed", 
                                         title = "All tissues only imputed", Points = T, alpha = 0.05,
                                         xl = "Mean of NZX and X", yl = "NZX - X"),
                mutate(.,
                       x = (NZTX + X) / 2,
                       y = (NZTX - X)) %$%
                  make_bland_altman_plot(x, y, imputed_for, fillname = "imputed", 
                                         title = "All tissues only imputed", Points = T, alpha = 0.05,
                                         xl = "Mean of NZTX and X", yl = "NZTX - X"))}



all.atlas.max_imputation_status %>% 
  filter(imputed) %>%
  gather(key = "Type", value = "Value", NX, X, IX, NZX, NZTX) %>%
  ggplot(aes(Value, fill = Type)) +
  geom_density(alpha = 0.2, color = "black") +
  simple_theme

all.atlas.max_imputation_status %>%
  gather(key = "Type", value = "Value", NX, X, IX, NZX, NZTX) %>%
  ggplot(aes(Value, fill = Type)) +
  geom_density(alpha = 0.5, color = "black") +
  simple_theme

# 
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


###################################################################
## Test for different categorization settings
###################################################################



library(pbapply)


enrich_fold_settings <- seq(3,7,0.5)
under_lim_settings <- seq(1,8,0.25)
group_num_settings <- c(6)#seq(5,11,1)

number_of_settings <- nrow(expand.grid(enrich_fold_settings, under_lim_settings, group_num_settings))

cat(paste(number_of_settings, "settings\nIt takes ~1 min per setting.\n"))


first <- T
pb <- timerProgressBar(min = 1, max = number_of_settings)
on.exit(close(pb))

i = 1
for(enrich_fold in enrich_fold_settings) {
  for(under_lim in under_lim_settings) {
    for(group_num in group_num_settings) {
      setTimerProgressBar(pb, i)
      i <- i + 1
      atlas_categories_temp <- 
        get.categories.with.num.expressed(all.atlas.max,
                                          max_column = "limma_gene_dstmm.zero.impute.expression_maxEx",
                                          cat_column = "consensus_content_name",
                                          enrich.fold = enrich_fold,
                                          under.lim = under_lim,
                                          group.num = group_num) %>%
        mutate(enrich_fold = enrich_fold,
               under_lim = under_lim,
               group_num = group_num)
      if(first){
        first <- F
        atlas_categories <- atlas_categories_temp
      } else {
        atlas_categories <- rbind(atlas_categories, atlas_categories_temp)
      }
    }
  }
}



atlas_categories %>%
  group_by(express.category.2, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(express.category.2)) %>%
  ggplot(aes(under_lim, number, label = number,color = express.category.2)) + 
  geom_point()+
  geom_line()+
  geom_text(vjust = -0.5)+
  scale_color_manual(values = expressed.cat.cols) + 
  simple_theme

ggsave(paste(result_folder, "settings expressed.png", sep = "/"), width = 8, height = 8)

atlas_categories %>%
  group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(elevated.category)) %>%
  ggplot(aes(under_lim, number, label = number,color = elevated.category)) + 
  geom_point()+
  geom_line()+
  geom_text(vjust = -0.5)+
  scale_color_manual(values = elevated.cat.cols) + 
  simple_theme

ggsave(paste(result_folder, "settings elevated.png", sep = "/"), width = 8, height = 8)

readr::write_delim(atlas_categories, paste(result_folder, "atlas_categories_test.txt", sep = "/"), delim = "\t")
