get.categories.with.num.expressed <- function(expMax.tb, max_column, cat_column, enrich.fold, group.num, under.lim = 1) {
  expMax.tb <- mutate(expMax.tb, 
                      max.norm.expression = eval(parse(text = max_column)),
                      categories = eval(parse(text = cat_column)))
  num_groups <- length(unique(expMax.tb$categories))
  
  category.mat <- 
    expMax.tb %>%
    mutate(max.norm.expression = round(max.norm.expression, 4)) %>%
    group_by(ensg_id) %>%
    summarise(mean.exp = mean(max.norm.expression, na.rm = T),
              min.exp = min(max.norm.expression, na.rm = T),
              max.exp = max(max.norm.expression, na.rm = T), 
              max.2nd = sort(max.norm.expression)[length(max.norm.expression)-1],
              num.exp = length(which(max.norm.expression >= under.lim)),
              frac.exp = num.exp/length(max.norm.expression[!is.na(max.norm.expression)])*100,
              lim = max.exp/enrich.fold, #5
              num_over = length(max.norm.expression[which(max.norm.expression > lim & max.norm.expression > under.lim)]), 
              mean_over = mean(max.norm.expression[which(max.norm.expression > lim & max.norm.expression > under.lim)]),
              min_over = min(max.norm.expression[which(max.norm.expression > lim & max.norm.expression > under.lim)]),
              max_under_lim = max(max.norm.expression[which(max.norm.expression < min_over)], under.lim*0.1),
              expressed.in = paste(sort(categories[which(max.norm.expression > lim & max.norm.expression > under.lim)]), collapse=", "),
              expressed.in.index = paste(sort(which(max.norm.expression > lim & max.norm.expression > under.lim)), collapse=", "),
              num.enh = length(which(max.norm.expression/mean.exp>=enrich.fold & max.norm.expression >= under.lim)), # 5
              enhanced.in = paste(sort(categories[which(max.norm.expression/mean.exp>=enrich.fold & max.norm.expression >= under.lim)]), collapse=", "), # 5
              enhanced.in.index =  paste(which(max.norm.expression/mean.exp>=enrich.fold & max.norm.expression >= under.lim), collapse=", "), # 5
              num.na = num_groups - length(max.norm.expression),
              max.2nd.or.lim = max(max.2nd, under.lim*0.1),
              tissues.under.lim = paste(sort(categories[which(max.norm.expression < under.lim)]), collapse=", "),
              tissues.under.lim.index =  paste(which(max.norm.expression < under.lim), collapse=", ")) %>%
    mutate(category = ifelse(# Genes not expressed above limit:
      num.exp == 0, 1,
      # Genes with expression 5 times more than anything else are tissue enriched:
      ifelse(max.exp/max.2nd.or.lim >= enrich.fold, 2, # 5
             # Genes with expression 5 times more than other tissues in groups of max group_num - 1 are group enriched:
             ifelse(max.exp >= lim & 
                      num_over < group.num & num_over > 1 &
                      mean_over/max_under_lim >= enrich.fold, 3, #5
                    # Genes with expression in tissues 5 times more than the mean are tissue enhance:
                    ifelse(num.enh > 0, 4, 
                           # Genes expressed in all tissues, mixed:
                           ifelse(frac.exp == 100, 6, 5))))),
      tissue.spec.score = case_when(category == 2 ~ max.exp/max.2nd.or.lim,
                                    category == 3 ~ mean_over/max_under_lim, 
                                    category == 4 ~ max.exp/mean.exp)) %>%
    left_join(tibble(category = 1:6, 
                     category.text = c("not detected", "tissue enriched", "group enriched", 
                                       "tissue enhanced", "mixed", "expressed in all tissues")), 
              by = "category") %$%
    tibble(ensg_id = ensg_id,
           category = category, 
           category.text = category.text,
           num.expressed = num.exp, 
           fraction.expressed = frac.exp,
           max.exp = max.exp,
           `enriched tissues` = ifelse(category %in% 4, enhanced.in, 
                                       ifelse(category %in% 2:3, expressed.in, "")),
           `enriched tissues index` = ifelse(category %in% 4, enhanced.in.index, 
                                             ifelse(category %in% 2:3, expressed.in.index, "")),
           `tissue/group specific score` = ifelse(category %in% 2:4, tissue.spec.score, ""),
           num.na = num.na,
           `tissues under lim` = tissues.under.lim,
           `tissues under lim index` = tissues.under.lim.index)
  
  category.mat <- 
    category.mat%>%
    mutate(express.category = case_when(num.expressed == max(num.expressed)                  ~ "expressed in all",
                                        num.expressed > 0 & num.expressed<max(num.expressed) ~ "expressed in some",
                                        num.expressed == 0                                   ~ "not expressed"),
           
           express.category.2 = case_when(fraction.expressed == 100                          ~ "expressed in all",
                                          fraction.expressed >= 31 & fraction.expressed<100  ~ "expressed in many",
                                          num.expressed > 1 & fraction.expressed<31          ~ "expressed in some",
                                          num.expressed == 1                                 ~ "expressed in single",
                                          num.expressed == 0                                 ~ "not expressed"),
           
           elevated.category = case_when(category.text == "tissue enriched"                  ~ "tissue enriched",
                                         category.text == "tissue enhanced"                  ~ "tissue enhanced",
                                         category.text == "group enriched"                   ~ "group enriched",
                                         category.text == "not detected"                     ~ "not detected",
                                         TRUE                                                ~ "low tissue specificity"),
           elevated.category.2 = case_when(category.text == "tissue enriched" | 
                                             category.text == "tissue enhanced" | 
                                             category.text == "group enriched"               ~ "general elevated",
                                           category.text == "not detected"                   ~ "not detected",
                                           TRUE                                              ~ "low tissue specificity")) 
  return(category.mat)
}	


calc_elevated.table <- function(tb.wide, atlas.categories, under.lim = 1, cat.colum = "category", colum = "enriched tissues") {
  elevated.table <- 
    matrix(0,nrow=nrow(atlas.categories), ncol=ncol(tb.wide), dimnames = dimnames(tb.wide)) 
  
  enriched.pos <- which(atlas.categories[, cat.colum][[1]]==2)
  group.pos <- which(atlas.categories[, cat.colum][[1]]==3)
  enhanced.pos <- which(atlas.categories[, cat.colum][[1]]==4)
  
  for(j in 1:ncol(tb.wide)) {
    tissue<-colnames(tb.wide)[j]
    pos<-grep(tissue,atlas.categories[, colum][[1]])
    
    table(atlas.categories[, cat.colum][[1]][pos])
    elevated.table[intersect(enriched.pos, pos),j]<-2
    elevated.table[intersect(group.pos, pos),j]<-3
    elevated.table[intersect(enhanced.pos, pos),j]<-4
    elevated.table[which(tb.wide[,j]<under.lim),j]<-1
    elevated.table[which(tb.wide[,j]<under.lim & atlas.categories[, cat.colum][[1]]==1),j]<-1.5
    elevated.table[which(elevated.table[,j]==0 & atlas.categories[, cat.colum][[1]]==6),j]<-6
    elevated.table[which(elevated.table[,j]==0),j] <-5
  }
  return(elevated.table)
}


calc_elevated.summary.table <- function(elevated.table, celltype = F) { 
  cat_names <- setNames(c("Not detected in this tissue","Not detected in any tissues","Tissue enriched",
                          "Group enriched","Tissue enhanced","Mixed in this tissue","Expressed in all tissues"), 
                        c(1, 1.5, 2, 3, 4, 5, 6))
  
  if(celltype) cat_names <- gsub("Tissue", "Celltype", gsub("tissue", "celltype", gsub(pattern = "tissues", "celltypes", cat_names)))
  
  elevated.table %>% 
    apply(2, FUN = function(x) table(factor(x, levels = c(1,1.5,2,3,4,5,6)))) %>%
    t() %>%
    set_colnames(cat_names)
}
