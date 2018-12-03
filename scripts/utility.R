## utilize
under_limit <- function(new.expression, expression, method) {
  tibble(new.expression, expression, method) %>%
    filter(method == "GTEx") %$%
    sort(new.expression)[{a <- sort(c(expression, 1)); which(a == 1)[1]}]
}


get_method <- function(method, norm.expression, max.norm.expression){
  met <- method[which(near(norm.expression,max.norm.expression))]
  
  if(length(met) > 1) {
    return(paste0(sort(unique(met), decreasing = T), collapse = "; "))
  } 
  return(met)
}

generate_wide <- function(expMax.tb, ensg_column, group_column, max_column) { #all.atlas.max, max_column
  expMax.tb %>%
    dplyr::select(group_column, ensg_column, max_column) %>%
    spread(key = group_column, value = max_column) %>%
    column_to_rownames(var=ensg_column) %>%
    round(4)
}
