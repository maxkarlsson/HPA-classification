
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

# If file exists: 
atlas_categories <- readr::read_delim("./results/atlas_categories_test.txt", delim = "\t")


atlas_categories %>%
  group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(express.category.2)) %>%
  ggplot(aes(under_lim, number, label = number, color = elevated.category)) + 
  geom_point(aes(size = enrich_fold))+
  geom_line()+
  geom_text(vjust = -0.5)+
  scale_color_manual(values = elevated.cat.cols) + 
  simple_theme

library(rgl)

make_3d_gif_rgl <- function(x, y, z, Group, Color = "black", GroupColor = NA, GroupTextColor = "black", GroupTextcex = 2, GroupTextJust = c(0,0,0),
                            bgcol = NA, 
                            sphere.r = 0.2, sphere = F,
                            legend.pos = "bottom", legend = F,  
                            legend.cex = 4, legend.ncol = 1, legend.bty = "n",
                            legend.text.col = "black",
                            axes.col = "black", 
                            save = F, savefile = "", 
                            ellipse = T, 
                            ellipse.subdivisions = 5, 
                            groupLabels = F, 
                            duration = 10, spin.axis = c(0, -1, 1), 
                            size = 1000, zoom = 1,
                            plot.axes = T, 
                            ellipse.alpha = 0.2, point.alpha = 1) {
  require(rgl)
  
  # Assert that input are correct formats
  stopifnot(!anyNA(x), is.numeric(x))
  stopifnot(!anyNA(y), is.numeric(y))
  stopifnot(!anyNA(z), is.numeric(z))
  
  # If GroupColor is NA, use Color
  if(is.na(GroupColor[[1]]) & length(GroupColor)) GroupColor <- Color
  
  if(save) {
    stopifnot(nchar(savefile) != 0, is.character(savefile))
    temp_savedir <- paste(savefile, "TEMP", sep = "_")
    dir.create(temp_savedir, showWarnings = FALSE)
  }
  
  # Calculate group means
  data_means <-
    tibble(x, y, z, Group) %>% 
    group_by(Group) %>%
    dplyr::summarise(x.mean = mean(x) + GroupTextJust[1],
                     y.mean = mean(y) + GroupTextJust[2],
                     z.mean = mean(z) + GroupTextJust[3])
  
  # Make clear plot space
  rgl::clear3d()
  
  par3d(windowRect = c(0, 0, size, size)) # make the window large
  par3d(zoom = zoom)
  
  if(sphere) {
    rgl.spheres(x, y, z, r = sphere.r, color = Color) 
  } else {
    plot3d(x, y, z, 
           col = Color, 
           size = 15, 
           axes = plot.axes, add = T, 
           alpha = point.alpha)
  }
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  
  if(plot.axes) {
    rgl.lines(lim(x), c(0, 0), c(0, 0), color = axes.col)
    rgl.lines(c(0, 0), lim(y), c(0, 0), color = axes.col)
    rgl.lines(c(0, 0), c(0, 0), lim(z), color = axes.col)
  }
  
  
  
  if(legend){
    bgplot3d({
      par(bg = bgcol)
      plot.new()
      legend(legend.pos, legend = unique(Group), pch = 16, col = unique(GroupColor), 
             cex=legend.cex, inset=c(0.02), ncol = legend.ncol, bty = legend.bty,
             text.col = legend.text.col)})
  } else {
    rgl.bg(color = bgcol)
  }
  
  if(ellipse){
    for (g in unique(Group)){
      tibble(x, y, z) %>%
        filter(Group==g) %>%
        {shade3d(ellipse3d(x = cov(.), centre = colMeans(.), level = 0.95, smooth = T, subdivide = ellipse.subdivisions), 
                 col = GroupColor[Group == g][1], alpha = 0.2)}
    }
  }
  
  
  if(groupLabels) {
    data_means %$%
      text3d(x.mean, y.mean, z.mean, texts = Group, col = GroupTextColor, cex = GroupTextcex)
  }
  
  spinner <- spin3d(axis = spin.axis, rpm = 60/duration)
  
  if(save) {
    movie3d(spinner, 
            duration = duration,
            dir = temp_savedir, 
            fps = 20)
  } else {
    play3d(spinner, 
           duration = duration)
  }
  
  if(save) {
    file.rename(from = paste(temp_savedir, "movie.gif", sep = "/"), 
                to = savefile)
    
    unlink(temp_savedir, recursive = T)
  }
  
}


atlas_categories %>%
  group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(express.category.2)) %>%
  ggplot(aes(under_lim, number, label = number, color = elevated.category)) + 
  geom_point(aes(size = enrich_fold))+
  geom_line()+
  geom_text(vjust = -0.5)+
  scale_color_manual(values = elevated.cat.cols) + 
  simple_theme

atlas_categories %>%
  group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(express.category.2)) %$% 
  make_3d_gif_rgl(x = enrich_fold/max(enrich_fold), y = under_lim/max(under_lim), z = number/max(number), 
                  Group = elevated.category, 
                  Color = elevated.cat.cols[match(elevated.category, names(elevated.cat.cols))], ellipse = F)

atlas_categories %>%
  group_by(express.category.2, enrich_fold, under_lim, group_num) %>%
  summarise(number = length(express.category.2)) %$% 
  make_3d_gif_rgl(x = enrich_fold/max(enrich_fold), y = under_lim/max(under_lim), z = number/max(number), 
                  Group = express.category.2, 
                  Color = expressed.cat.cols[match(express.category.2, names(expressed.cat.cols))], ellipse = F)


for(enrich_fold_setting in unique(atlas_categories$enrich_fold)) {
  atlas_categories %>%
    filter(enrich_fold == enrich_fold_setting) %>%
    group_by(express.category.2, enrich_fold, under_lim, group_num) %>%
    summarise(number = length(express.category.2)) %>%
    ggplot(aes(under_lim, number, label = number,color = express.category.2)) + 
    geom_point()+
    geom_line()+
    geom_text(vjust = -0.5)+
    scale_color_manual(values = expressed.cat.cols) + 
    simple_theme
  
  ggsave(paste(result_folder, paste0("settings expressed fold ", format(enrich_fold_setting, nsmall = 2), ".png"), sep = "/"), width = 16, height = 8)
  
  atlas_categories %>%
    filter(enrich_fold == enrich_fold_setting) %>%
    group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
    summarise(number = length(elevated.category)) %>%
    ggplot(aes(under_lim, number, label = number,color = elevated.category)) + 
    geom_point()+
    geom_line()+
    geom_text(vjust = -0.5)+
    scale_color_manual(values = elevated.cat.cols) + 
    simple_theme
  
  ggsave(paste(result_folder, paste0("settings elevated fold ", format(enrich_fold_setting, nsmall = 2), ".png"), sep = "/"), width = 16, height = 8)
}

for(under_lim_setting in unique(atlas_categories$under_lim)) {
  atlas_categories %>%
    filter(under_lim == under_lim_setting) %>%
    group_by(express.category.2, enrich_fold, under_lim, group_num) %>%
    summarise(number = length(express.category.2)) %>%
    ggplot(aes(enrich_fold, number, label = number,color = express.category.2)) + 
    geom_point()+
    geom_line()+
    geom_text(vjust = -0.5)+
    scale_color_manual(values = expressed.cat.cols) + 
    simple_theme
  
  ggsave(paste(result_folder, paste0("settings expressed under lim ", format(under_lim_setting, nsmall = 2), ".png"), sep = "/"), width = 16, height = 8)
  
  atlas_categories %>%
    filter(under_lim == under_lim_setting) %>%
    group_by(elevated.category, enrich_fold, under_lim, group_num) %>%
    summarise(number = length(elevated.category)) %>%
    ggplot(aes(enrich_fold, number, label = number,color = elevated.category)) + 
    geom_point()+
    geom_line()+
    geom_text(vjust = -0.5)+
    scale_color_manual(values = elevated.cat.cols) + 
    simple_theme
  
  ggsave(paste(result_folder, paste0("settings elevated under lim ", format(under_lim_setting, nsmall = 2), ".png"), sep = "/"), width = 16, height = 8)
}




#-------- Demonstrate the differences with different settings ------

atlas_categories <- readr::read_delim("./results/atlas_categories_test.txt", delim = "\t")

all.atlas.category.4 <- 
  atlas_categories %>%
  filter(., under_lim == 3 & enrich_fold == 4)

all.atlas.category.5 <- 
  atlas_categories %>%
  filter(., under_lim == 3 & enrich_fold == 5)


make_swarm_expression_plot(atlas.max = all.atlas.max, 
                           atlas.cat = all.atlas.category.4, 
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues fold 4')

make_swarm_expression_plot(atlas.max = all.atlas.max, 
                           atlas.cat = all.atlas.category.5, 
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues fold 5')


make_swarm_expression_plot(atlas.max = filter(all.atlas.max, 
                                              ensg_id %in% filter(all.atlas.category.4, 
                                                                  elevated.category == "tissue enriched")$ensg_id), 
                           atlas.cat = all.atlas.category.4, 
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues fold 4 only tissue enriched')

make_swarm_expression_plot(atlas.max = filter(all.atlas.max, 
                                              ensg_id %in% filter(all.atlas.category.5, 
                                                                  elevated.category == "tissue enriched")$ensg_id), 
                           atlas.cat = all.atlas.category.5, 
                           maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                           tissue_column = "consensus_content_name",
                           outpath = result_folder,
                           prefix = 'all_tissues fold 5 only tissue enriched')

make_swarm_expression_plot_difference <- function(atlas.max, atlas.cat1, atlas.cat2, maxEx_column, tissue_column, outpath, prefix, label.all = F) {
  genes <- 
    left_join(atlas.cat1, atlas.cat2, by = "ensg_id", suffix = c(".4", ".5")) %>%
    filter(elevated.category.4 == "tissue enriched" & elevated.category.5 != "tissue enriched") %$% ensg_id
  
  plot.data <- 
    atlas.max %>%
    ungroup() %>%
    mutate(expression = eval(parse(text = maxEx_column)),
           Grouping = eval(parse(text = tissue_column))) %>%
    dplyr::select(Grouping, ensg_id, expression) %>%
    #Plot 1 % highest
    filter(ensg_id %in% genes) %>%
    group_by(Grouping) %>%
    # Write gene name if highest 2 % per tissue and/or highest 1 % in total
    mutate(highest = expression >= quantile(expression, probs = 0.99)) %>%
    ungroup() %>%
    mutate(highest = ifelse(highest, T, expression >= quantile(expression, probs = 0.98)),
           highest = ifelse(rep(label.all, length(highest)), T, highest)) %>%
    left_join(dplyr::select(atlas.cat1, ensg_id, `enriched tissues`, express.category.2, elevated.category), by = "ensg_id") %>%
    filter(`enriched tissues` == Grouping) %>%
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

make_swarm_expression_plot_difference(atlas.max = all.atlas.max,
                                      atlas.cat1 = all.atlas.category.4, 
                                      atlas.cat2 = all.atlas.category.5,
                                      maxEx_column = "limma_gene_dstmm.zero.impute.expression_maxEx", 
                                      tissue_column = "consensus_content_name",
                                      outpath = result_folder,
                                      prefix = 'all_tissues fold 4 to 5 only tissue enriched all', label.all = T)

left_join(all.atlas.category.4, all.atlas.category.5, by = "ensg_id", suffix = c(".4", ".5")) %>%
  filter(elevated.category.4 == "tissue enriched" & elevated.category.5 != "tissue enriched", 
         elevated.category.4 == "tissue enriched" & elevated.category.5 != "tissue enriched", 
         elevated.category.4 == "tissue enriched" & elevated.category.5 != "tissue enriched")


plot_classification_chord <- function(cat1, cat2, savename = "") {
  png(paste(result_folder, paste0("chord not norm vs norm.png ", savename, ".png"), sep = "/"), width = 12, height = 12, units = "in", res = 300)
  unique_categories <- unique(c(cat1$category.text, cat2$category.text))
  plot.colors <- cat.cols[which(names(cat.cols) %in% unique_categories)]
  
  left_join(cat1, cat2, by = "ensg_id") %>%
    group_by(category.x, category.text.x, category.y, category.text.y) %>%
    summarise(number = length(category.x)) %>%
    ungroup() %>%
    mutate(category.text.x = paste("from", category.text.x),
           category.text.y = paste("to", category.text.y)) %$%
    chord_classification(from = category.text.x, 
                         to = category.text.y, 
                         sizes = number,
                         grid.col = c(setNames(plot.colors, paste("from", names(plot.colors))), setNames(plot.colors, paste("to", names(plot.colors)))),
                         groups = c(rep(1, length(plot.colors)),
                                    rep(2, length(plot.colors))), 
                         plot.order = c(paste("from", names(plot.colors)),
                                        rev(paste("to", names(plot.colors)))))
  
  dev.off()
}

plot_express_elevated_classification_chord <- function(cat1, cat2, savename = "") {
  png(paste(result_folder, paste0("chord vs expressed ", savename, ".png"), sep = "/"), width = 12, height = 12, units = "in", res = 300)
  plot.data <- 
    left_join(cat1, cat2, by = "ensg_id") %>%
    group_by(express.category.2.x, express.category.2.y) %>%
    summarise(number = length(express.category.2.x)) %>%
    ungroup() %>%
    mutate(express.category.2.text.x = paste("from", express.category.2.x),
           express.category.2.text.y = paste("to", express.category.2.y)) 
  
  plot.colors <- c(setNames(expressed.cat.cols, paste("from", names(expressed.cat.cols))),
                   setNames(expressed.cat.cols, paste("to", names(expressed.cat.cols))))
  
  plot.groups <- c(rep(1, length(expressed.cat.cols)),
                   rep(2, length(expressed.cat.cols)))
  
  plot.order <- c(paste("from", names(expressed.cat.cols)),
                  rev(paste("to", names(expressed.cat.cols))))
  
  plot.data %$%
    chord_classification(from = express.category.2.text.x, 
                         to = express.category.2.text.y, 
                         sizes = number,
                         grid.col = plot.colors,
                         groups = plot.groups, 
                         plot.order = plot.order)
  
  dev.off()
  
  png(paste(result_folder, paste0("chord vs elevated ", savename, ".png"), sep = "/"), width = 12, height = 12, units = "in", res = 300)
  plot.data <- 
    left_join(cat1, cat2, by = "ensg_id") %>%
    group_by(elevated.category.x, elevated.category.y) %>%
    summarise(number = length(elevated.category.x)) %>%
    ungroup() %>%
    mutate(elevated.category.text.x = paste("from", elevated.category.x),
           elevated.category.text.y = paste("to", elevated.category.y)) 
  
  
  cats <- c('tissue enriched', 'group enriched', 'tissue enhanced', 'low tissue specificity', 'not detected')
  cats.cols <- elevated.cat.cols[match(cats, names(elevated.cat.cols))]
  plot.colors <- c(setNames(cats.cols, paste("from", cats)),
                   setNames(cats.cols, paste("to", cats)))
  
  
  
  plot.groups <- c(rep(1, length(cats)),
                   rep(2, length(cats)))
  
  plot.order <- c(paste("from", cats),
                  rev(paste("to", cats)))
  
  plot.data %$%
    chord_classification(from = elevated.category.text.x, 
                         to = elevated.category.text.y, 
                         sizes = number,
                         grid.col = plot.colors,
                         groups = plot.groups, 
                         plot.order = plot.order)
  
  dev.off()
}

plot_classification_chord(all.atlas.category.5, all.atlas.category.4, savename = "fold 4 vs fold 4")

plot_express_elevated_classification_chord(all.atlas.category.5, all.atlas.category.4, savename = "fold 4 vs fold 4")


all.atlas.max.wide <- generate_wide(all.atlas.max, ensg_column='ensg_id', 
                                    group_column='consensus_content_name', 
                                    max_column="limma_gene_dstmm.zero.impute.expression_maxEx")

all.atlas.max.wide %>%
  calc_elevated.table(atlas.categories = all.atlas.category.4) %>%
  apply(2, FUN = function(x) table(factor(x, levels = c(1,1.5,2,3,4,5,6)))) %>%
  t() %>%
  set_colnames(setNames(c("Not detected in this tissue","Not detected in any tissues","Tissue enriched",
                          "Group enriched","Tissue enhanced","Mixed in this tissue","Expressed in all tissues"), 
                        c(1, 1.5, 2, 3, 4, 5, 6))) %>%
  make_elevated_bar_plot(outpath = result_folder,
                         prefix = 'all_tissues fold 4')

all.atlas.max.wide %>%
  calc_elevated.table(atlas.categories = all.atlas.category.5) %>%
  apply(2, FUN = function(x) table(factor(x, levels = c(1,1.5,2,3,4,5,6)))) %>%
  t() %>%
  set_colnames(setNames(c("Not detected in this tissue","Not detected in any tissues","Tissue enriched",
                          "Group enriched","Tissue enhanced","Mixed in this tissue","Expressed in all tissues"), 
                        c(1, 1.5, 2, 3, 4, 5, 6))) %>%
  make_elevated_bar_plot(outpath = result_folder,
                         prefix = 'all_tissues fold 5')
###################################################################
## Checking genes <1 for HPA
###################################################################

HPA_below_1 <-
  all.atlas %>%
  filter(expression < 1 & method == "HPA") %>%
  select(ensg_id, content_name) %>%
  unique() %>%
  left_join(all.atlas, by = c("ensg_id", "content_name")) %>% 
  left_join(all.atlas.category, by = "ensg_id") %>%
  rename(X = 5, NX = 15) %>%
  gather(key = "Type", value = "Value", X, NX) %>%
  
  mutate(Type = factor(Type, levels = c("X", "NX")),
         #tissue = factor(tissue, levels = unique(tissue[order(group)])),
         label = paste(express.category.2, 
                       elevated.category, sep = "\n"))
genes_not_detected <- 
  unique(HPA_below_1$ensg_id)

for(i in seq(1, length(genes_not_detected), 5)) {
  j <- i + 5
  if(j > length(genes_not_detected)) j <- length(genes_not_detected)
  HPA_below_1 %>%
    filter(ensg_id %in% genes_not_detected[i:j]) %>%
    ggplot(aes(content_name, Value, fill = method, group = method))+
    geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
    geom_hline(yintercept = 2, color = NA)+
    geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
    geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
    #annotate("text", y= max(plot.data$Value), x =1,label="Custom Title",hjust=1) +
    
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.text.x = element_blank()) +
    scale_fill_manual(values = dataset.colors)+
    facet_grid(ensg_id ~ Type, scales = "free") 
  ggsave(paste(result_folder, paste0("X NZX bars genes ", i, "-", j, ".png"), sep = "/"), width = 40, height = 20, dpi = 300, units = "cm")
  
}

all_atlas_cat <- 
  all.atlas %>%
  left_join(all.atlas.category, by = "ensg_id") %>%
  rename(X = expression, TMM = dstmm.zero.expression, NGX = gene_dstmm.zero.impute.expression, NX = limma_gene_dstmm.zero.impute.expression) %>%
  mutate(X_low = X < 1)



genes_in_GTEXFANTOM <- 
  all.atlas %>%
  group_by(ensg_id) %>%
  summarise(imp = any(is.na(expression))) %>%
  filter(!imp)
# 
# all.atlas %>%
#   filter(ensg_id %in% genes_not_in_GTEXFANTOM$ensg_id) %>%
#   filter(expression < 1 & method == "HPA") %>%
#   select(ensg_id, content_name) %>%
#   unique() %>%
#   left_join(all.atlas, by = c("ensg_id", "content_name")) %>% 
#   left_join(all.atlas.category, by = "ensg_id") %>%
#   rename(X = 5, NX = 15) %>%
#   
#   ggplot(aes(X + 1, NX + 1, color = method, group = method))+
#   geom_hline(yintercept = 2, color = "red", linetype = "dashed")+
#   geom_vline(xintercept = 2, color = "red", linetype = "dashed")+
#   geom_point(alpha = 0.1)+
#   
#   simple_theme+
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         axis.text.x = element_blank()) +
#   scale_x_log10() + 
#   scale_y_log10() + 
#   scale_color_manual(values = dataset.colors)

all_atlas_cat %>%
  filter(ensg_id %in% genes_not_in_GTEXFANTOM) %>%
  filter(X < 1 & method == "HPA") %>%
  ggplot(aes(X, NX, color = method, group = method))+
  geom_hline(yintercept = c(1, 2, 3), color = "red", linetype = "dashed")+
  geom_vline(xintercept = 1, color = "red", linetype = "dashed")+
  geom_point(alpha = 0.1)+
  simple_theme+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.text.x = element_blank()) +
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values = dataset.colors)

genes_not_in_GTEXFANTOM <- 
  all.atlas %>%
  filter(method %in% c("FANTOM", "GTEx")) %>%
  filter(is.na(expression)) %>%
  group_by(ensg_id) %>%
  summarise(n = length(unique(method))) %>%
  filter(n == 2) %$% 
  unique(ensg_id)

top_values <- 
  all_atlas_cat %>%
  filter(ensg_id %in% genes_not_in_GTEXFANTOM) %>% 
  filter(X < 1 & method == "HPA") %>% filter(NX>3) %$% 
  ensg_id[order(NX, decreasing = T)]  %>% unique() #%>%
head(20)

for( i in seq(1, length(top_values), 5)) {
  j = i + 4
  if(j > length(top_values)) j = length(top_values)
  
  all_atlas_cat  %>%
    #filter(X < 1 & method == "HPA") %>% 
    gather(key = "Type", value = "Value", X, TMM, NGX, NX)%>%
    filter(ensg_id %in% top_values[i:j]) %>%
    #filter(ensg_id == "ENSG00000276017") %>% View
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
           #tissue = factor(tissue, levels = unique(tissue[order(group)])),
           label = paste(express.category.2, 
                         elevated.category, sep = "\n")) %>%
    ggplot(aes(content_name, Value, fill = method, group = method))+
    geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
    geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
    geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
    #annotate("text", y= max(plot.data$Value), x =1,label="Custom Title",hjust=1) +
    
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = dataset.colors)+
    facet_grid(ensg_id ~ Type, scales = "free") 
  
  ggsave(paste(result_folder, paste0(i, "-", j, " top TMM.png"), sep = "/"), width = 20, height = 10)
  
}

randomgenes <- sample(unique(all.atlas$ensg_id), 100)
for( i in seq(1, length(randomgenes), 5)) {
  j = i + 4
  if(j > length(randomgenes)) j = length(randomgenes)
  
  all_atlas_cat  %>%
    #filter(X < 1 & method == "HPA") %>% 
    gather(key = "Type", value = "Value", X, TMM, NGX, NX)%>%
    filter(ensg_id %in% randomgenes[i:j]) %>%
    #filter(ensg_id == "ENSG00000276017") %>% View
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
           #tissue = factor(tissue, levels = unique(tissue[order(group)])),
           label = paste(express.category.2, 
                         elevated.category, sep = "\n")) %>%
    ggplot(aes(content_name, Value, fill = method, group = method))+
    geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
    geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
    geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
    #annotate("text", y= max(plot.data$Value), x =1,label="Custom Title",hjust=1) +
    
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = dataset.colors)+
    facet_grid(ensg_id ~ Type, scales = "free") 
  
  ggsave(paste(result_folder, paste0(i, "-", j, " random TMM.png"), sep = "/"), width = 20, height = 10)
  
}


AsaGener <- 
  c("ENSG00000185686","ENSG00000006016", "ENSG00000000419","ENSG00000006744","ENSG00000006757",
    "ENSG00000007047","ENSG00000007520","ENSG00000008441","ENSG00000240747",
    "ENSG00000124702","ENSG00000130695","ENSG00000078295","ENSG00000274897",
    "ENSG00000233954","ENSG00000278599","ENSG00000260300","ENSG00000126266",
    "ENSG00000160948")

for( i in seq(1, length(AsaGener), 5)) {
  j = i + 4
  if(j > length(AsaGener)) j = length(AsaGener)
  
  all_atlas_cat  %>%
    #filter(X < 1 & method == "HPA") %>% 
    gather(key = "Type", value = "Value", X, NX)%>%
    filter(ensg_id %in% AsaGener[i:j]) %>%
    #filter(ensg_id == "ENSG00000276017") %>% View
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
           #tissue = factor(tissue, levels = unique(tissue[order(group)])),
           label = paste(express.category.2, 
                         elevated.category, sep = "\n")) %>%
    ggplot(aes(content_name, Value, fill = method, group = method))+
    geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
    geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
    geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
    #annotate("text", y= max(plot.data$Value), x =1,label="Custom Title",hjust=1) +
    
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
    scale_fill_manual(values = dataset.colors)+
    facet_grid(ensg_id ~ Type, scales = "free") 
  
  ggsave(paste(result_folder, paste0(i, "-", j, " Åsa TMM.png"), sep = "/"), width = 20, height = 10)
  
}

essential_genes <- readr::read_delim("./doc/essential genes.txt", delim = "\t") %$% ENSG

for( i in seq(1, length(essential_genes), 4)) {
  j = i + 3
  if(j > length(essential_genes)) j = length(essential_genes)
  
  all_atlas_cat  %>%
    #filter(X < 1 & method == "HPA") %>% 
    gather(key = "Type", value = "Value", X, NX)%>%
    filter(ensg_id %in% essential_genes[i:j]) %>%
    #filter(ensg_id == "ENSG00000276017") %>% View
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
           #tissue = factor(tissue, levels = unique(tissue[order(group)])),
           label = paste(express.category.2, 
                         elevated.category, sep = "\n")) %>%
    ggplot(aes(content_name, Value, fill = method, group = method))+
    
    #geom_vline(xintercept = c(5), size = 5, alpha = 0.3, color = "pink")+
    geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
    
    geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
    annotate(geom = "rect", xmin=c(56, 71)-0.5,xmax=c(56, 71)+0.5,ymin=-Inf,ymax=Inf, alpha=0.1,fill="green")+
    geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
    #annotate("text", y= max(plot.data$Value), x =1,label="Custom Title",hjust=1) +
    
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
    scale_fill_manual(values = dataset.colors)+
    facet_grid(ensg_id ~ Type, scales = "free") 
  
  ggsave(paste(result_folder, paste0("essential genes", i, "-", j, ".png"), sep = "/"), width = 20, height = 10)
  
}

for(meth in c("HPA", "GTEx", "FANTOM", "Blood")) {
  all_atlas_cat %>%
    filter(ensg_id %in% essential_genes) %>%
    filter(method == meth) %>%
    gather(key = "Type", value = "Value", X, NX, TMM)%>%
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX"))) %>%
    ggplot(aes(content_name, Value + 1, fill = method, color = method))+
    geom_violin(draw_quantiles = 0.5, alpha = 0.5)+
    simple_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
    scale_fill_manual(values = dataset.colors)+
    scale_color_manual(values = dataset.colors)+
    facet_wrap( ~ Type, scales = "free", ncol = 1) +
    scale_y_log10()
  
  ggsave(paste(result_folder, paste0("essential genes boxplot ", meth, ".png"), sep = "/"), width = 20, height = 10)
}


readr::write_delim(tibble(top_values), paste(result_folder, "Topgener.txt", sep = "/"), delim  = "\t")


all.atlas %>%
  filter(expression < 1 & method == "HPA") %>%
  select(ensg_id, content_name) %>%
  unique() %>%
  left_join(all.atlas, by = c("ensg_id", "content_name")) %>% 
  left_join(all.atlas.category, by = "ensg_id") %>%
  rename(X = 5, NX = 15) %>%
  
  ggplot(aes(X + 1, NX + 1, fill = method, alpha = ..count..))+
  geom_hline(yintercept = 2, color = "red", linetype = "dashed")+
  geom_vline(xintercept = 2, color = "red", linetype = "dashed")+
  geom_hex(bins = 100)+
  
  simple_theme+
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values = dataset.colors)


AsaGener <- 
  c("ENSG00000185686","ENSG00000006016", "ENSG00000000419","ENSG00000006744","ENSG00000006757",
    "ENSG00000007047","ENSG00000007520","ENSG00000008441","ENSG00000240747",
    "ENSG00000124702","ENSG00000130695","ENSG00000078295","ENSG00000274897",
    "ENSG00000233954","ENSG00000278599","ENSG00000260300","ENSG00000126266",
    "ENSG00000160948")
#### Generate pdf with ALL genes

genes <- unique(all_atlas_cat$ensg_id)

pdf(paste(result_folder, "All genes bar.pdf", sep = "/"), width = 16, height = 8)

pb <- timerProgressBar(min = 1, max = length(genes))
on.exit(close(pb))


p <- 
  all_atlas_cat  %>%
  gather(key = "Type", value = "Value", X, NX) %>%
  mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
         #tissue = factor(tissue, levels = unique(tissue[order(group)])),
         label = paste(express.category.2, 
                       elevated.category, sep = "\n")) %>%
  {lapply(seq(1, length(genes), 3), length(genes),
          FUN = function(i, jmax) {
            setTimerProgressBar(pb, i)
            j <- i + 2
            if(j>jmax) j <- jmax
            filter(., ensg_id %in% genes[i:j]) %>%
              ggplot(aes(content_name, Value, fill = method, group = method))+
              geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
              
              geom_bar(stat = "identity", show.legend = F, color = "black", position = "dodge")+
              annotate(geom = "rect", xmin=c(56, 71)-0.5,xmax=c(56, 71)+0.5,ymin=-Inf,ymax=Inf, alpha=0.1,fill="green")+
              geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
              
              simple_theme+
              theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
              scale_fill_manual(values = dataset.colors)+
              facet_grid(ensg_id ~ Type, scales = "free") 
            }) 
  } 

print(p)


dev.off()


#### Generate pdf with ÅSA genes

pdf(paste(result_folder, "Åsa genes bar.pdf", sep = "/"), width = 32, height = 8)

pb <- timerProgressBar(min = 1, max = length(AsaGener))
on.exit(close(pb))


p <- 
  all_atlas_cat  %>%
  gather(key = "Type", value = "Value", X, NX, TMM) %>%
  mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX")),
         #tissue = factor(tissue, levels = unique(tissue[order(group)])),
         label = paste(express.category.2, 
                       elevated.category, sep = "\n")) %>%
                       {lapply(seq(1, length(AsaGener), 3), length(AsaGener),
                               FUN = function(i, jmax) {
                                 setTimerProgressBar(pb, i)
                                 j <- i + 2
                                 if(j>jmax) j <- jmax
                                 filter(., ensg_id %in% AsaGener[i:j]) %>%
                                   ggplot(aes(content_name, Value, fill = method, group = method))+
                                   geom_hline(yintercept = 1:4, color = "red", linetype = "dashed")+
                                   
                                   geom_bar(stat = "identity", show.legend = F, size = 1, color = "black", position = "dodge")+
                                   annotate(geom = "rect", xmin=c(56, 71)-0.5,xmax=c(56, 71)+0.5,ymin=-Inf,ymax=Inf, alpha=0.1,fill="green")+
                                   geom_text(aes(3, 1, label = label), vjust = -1, hjust = 0, size = 2)+
                                   
                                   simple_theme+
                                   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8)) +
                                   scale_fill_manual(values = dataset.colors)+
                                   facet_grid(ensg_id ~ Type, scales = "free") 
                               }) 
                       } 

print(p)


dev.off()



##----- generate ÅSA layout pdf
tissue_colors <- readr::read_delim("ref/colors_92.tsv", delim = "\t")


print_gene_bar_pdf <- function(genes, savename) {
  
  pdf(paste(result_folder, savename, sep = "/"), width = 32, height = 8)
  plot.col <- with(tissue_colors, setNames(c(color, color, elevated.cat.cols), c(tissue_name, organ_name, names(elevated.cat.cols))))
  
  pb <- timerProgressBar(min = 1, max = length(genes))
  on.exit(close(pb))
  
  plot.data <- 
    all_atlas_cat  %>%
    filter(!imputed) %>%
    filter(method != "Blood") %>%
    
    gather(key = "Type", value = "Value", X, NX, TMM) %>%
    select(consensus_content_name, content_name, ensg_id, method, Type, Value, express.category.2, elevated.category) %>%
    
    rbind(left_join(all_atlas_cat  %>%
                      filter(!imputed) %>%
                      filter(method != "Blood") %>%
                      select(consensus_content_name, content_name, ensg_id, express.category.2, elevated.category),
                    all.atlas.max %>% 
                      select(consensus_content_name, ensg_id, 
                             expression_maxEx, 
                             dstmm.zero.expression_maxEx,
                             limma_gene_dstmm.zero.impute.expression_maxEx) %>% 
                      
                      rename(X = expression_maxEx,
                             TMM = dstmm.zero.expression_maxEx,
                             NX = limma_gene_dstmm.zero.impute.expression_maxEx) %>%
                      
                      mutate(method = "consensus"), 
                    by = c("ensg_id", "consensus_content_name")) %>%
            
            gather(key = "Type", value = "Value", X, NX, TMM)) %>%
    
    mutate(Type = factor(Type, levels = c("X", "TMM", "NGX", "NX", "consensus")),
           content_name = factor(content_name, 
                                 levels = select(., content_name, consensus_content_name) %>% 
                                   unique() %$% content_name[order(consensus_content_name)]),
           label = paste(express.category.2, 
                         elevated.category, sep = "\n"),
           label = ifelse(Type == "X" & method == "HPA", label, NA)#,consensus_Value = ifelse(Type == "NX", consensus_Value, NA)
    ) 
  
  p <- 
    plot.data %>%
    {lapply(1:length(genes),
            FUN = function(i) {
              setTimerProgressBar(pb, i)
              filter(., ensg_id %in% genes[i]) %>%
                ggplot(aes(content_name, Value, fill = consensus_content_name))+
                geom_hline(yintercept = 1:3, color = "red", linetype = "dashed")+
                
                geom_bar(stat = "identity", show.legend = F, size = 0.5, color = "black", position = "dodge")+
                
                geom_text(aes(5, 1, label = label), vjust = -1, hjust = 0, size = 5)+
                ggtitle(genes[i])+
                simple_theme+
                theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.8),
                      strip.text.x = element_text(size = 15),
                      strip.text.y = element_text(size = 15)) +
                scale_fill_manual(values = plot.col)+#dataset.colors)+
                facet_grid(method ~ Type)  
            })} 
  
  pb <- timerProgressBar(min = 1, max = length(genes))
  on.exit(close(pb))
  
  for(i in 1:length(genes)){
    setTimerProgressBar(pb, i)
    print(p[[i]])
  }
  
  
  
  dev.off()
}



print_gene_bar_pdf(AsaGener, "Åsas gener.pdf")
genes <- unique(all_atlas_cat$ensg_id)
print_gene_bar_pdf(genes, "Alla gener.pdf")


################################################################
# Plots for checking transfer of categories
################################################################



# Create classification based on only HPA

all.atlas.max.HPA <-
  all.atlas %>%
  # Remove genes that are imputed
  filter(!imputed) %>%
  # Keep only HPA
  filter(method == "HPA") %>%
  group_by(consensus_content_name, ensg_id) %>% 
  mutate(method = as.character(method)) %>%
  dplyr::summarise_at(.funs = funs(maxEx = max(., na.rm = T),
                                   method = get_method(method, ., max(., na.rm = T))),
                      .vars = grep("expression$", colnames(.), value = T)) 
# 

## Step 4. Category
all.atlas.category.HPA <- 
  get.categories.with.num.expressed(all.atlas.max.HPA,
                                    max_column = "expression_maxEx",
                                    cat_column = "consensus_content_name",
                                    enrich.fold = 5,
                                    under.lim = 1,
                                    group.num = 6)


plot.data <- 
  all.atlas.category %>%
  left_join(all.atlas.category.HPA, by = "ensg_id", suffix = c("", ".HPA")) %>%
  filter(elevated.category == "tissue enriched" | elevated.category.HPA == "tissue enriched") %>%
  mutate(tissue = ifelse(elevated.category == "tissue enriched", `enriched tissues`, elevated.category),
         tissue.HPA = ifelse(elevated.category.HPA == "tissue enriched", `enriched tissues.HPA`, elevated.category.HPA)) %>%
  group_by(tissue, tissue.HPA) %>%
  summarise(n = length(tissue)) 

plot.col <- with(tissue_colors, setNames(c(color, color, elevated.cat.cols), c(tissue_name, organ_name, names(elevated.cat.cols))))
plot.data %>%
  ungroup() %>%
  rename(from = tissue, to = tissue.HPA, sizes = n) %>%
  chordDiagram(annotationTrack = "grid", 
               preAllocateTracks = 1, 
               grid.col = c(plot.col[match(gsub(" 1", "", unique(c(.$from, .$to))), 
                                           names(plot.col))], "blood" = "dark red"))

circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, cex=0.6,adj = c(0, 0.3))
  #  circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
}, bg.border = NA)
