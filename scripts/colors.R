### hello world! ###

cat.cols <- c("expressed in all tissues" = "#4daf4a",
              "mixed" = "#377eb8"   ,
              "tissue enhanced" = "#984ea3",
              "group enriched" = "#FF9D00",
              "tissue enriched" = "#e41a1c",
              "not detected" = "gray",
              "expressed in all celltypes" = "#4daf4a",
              "celltype enhanced" = "#984ea3",
              "celltype enriched" = "#e41a1c")

cat2.cols <- c("Expressed in all tissues" = "#377eb8",
               "Expressed in all celltypes" = "#377eb8",
               "Mixed in this tissue" = "#4daf4a",
               "Mixed in this celltypes" = "#4daf4a",
               "Tissue enhanced" = "#984ea3",
               "Celltype enhanced" = "#984ea3",
               "Group enriched" = "#FF9D00",
               "Tissue enriched" = "#e41a1c",
               "Celltype enriched" = "#e41a1c",
               "Not detected in any celltypes" = "grey",
               "Not detected in any celltypes" = "grey",
               "Not detected in this celltype" = "dark gray",
               "Not detected in this celltype" = "dark gray")

dataset.colors <- setNames(c("dark red", "red","blue","green3"), c("Blood", "HPA","GTEx","FANTOM"))

elevated.cat.cols <- rev(c("not detected" = "grey",
                           "low tissue specificity" = "grey40",
                           "low celltype specificity" = "grey40",
                           "tissue enhanced" = "#984ea3",
                           "celltype enhanced" = "#984ea3",
                           "group enriched" = "#FF9D00",
                           "tissue enriched" = "#e41a1c",
                           "celltype enriched" = "#e41a1c"))

expressed.cat.cols <- c("expressed in all" = "#253494",
                        "expressed in many" = "#2c7fb8",
                        "expressed in some" = "#41b6c4",
                        "expressed in single" = "#a1dab4",
                        "not expressed" = "grey", 
                        "detected in all" = "#253494",
                        "detected in many" = "#2c7fb8",
                        "detected in some" = "#41b6c4",
                        "detected in single" = "#a1dab4",
                        "not detected " = "grey")

protein.class.palette <- c("secreted" = '#911D51',
                           "membrane" = '#6D4BAA', 
                           "other" = '#008490', 
                           "cd_marker" = '#318F1E', 
                           "transcription_factors" = '#B8801B', 
                           "mitochondrial" = '#E371B4', 
                           "ribosomal" = '#89A0F3', 
                           "none" = "black",
                           '#00C9BC', '#97C542', '#FFA05E')

protein.localization.palette2 <- c("membrane" = "#CE70A4",
                                  "secreted" = 	"#FCAC3B",
                                  "membrane and secreted isoforms" = "#755A85")

protein.localization.palette <- c("intracellular and membrane isoforms" = "#858141",
                                  "membrane" = "#6DB9C6",
                                  "intracellular" = "#FCAC3B",
                                  "secreted" = 	"#CE70A4",
                                  "intracellular and secreted isoforms" = "#CF5734",
                                  "membrane and secreted isoforms" = "#755A85",
                                  "intracellular, membrane, secreted isoforms" = "#794A39")

monaco_schmiedel_lineage <- 
  c('basophil' = 'granulocytes',
    'Central memory CD8 T-cell' = 'T-cells',
    'classical monocyte' = 'monocytes',
    'Effector memory CD8 T-cell' = 'T-cells',
    'Exhausted memory B-cell' = 'B-cells',
    'intermediate monocyte' = 'monocytes',
    'MAIT T-cell' = 'T-cells',
    'Memory CD4 T-cell TFH' = 'T-cells',
    'Memory CD4 T-cell Th1' = 'T-cells',
    'Memory CD4 T-cell Th1/Th17' = 'T-cells',
    'Memory CD4 T-cell Th17' = 'T-cells',
    'Memory CD4 T-cell Th2' = 'T-cells',
    'myeloid DC' = 'dendritic cells',
    'naive B-cell' = 'B-cells',
    'naive CD4 T-cell' = 'T-cells',
    'naive CD8 T-cell' = 'T-cells',
    'neutrophil' = 'granulocytes',
    'NK-cell' = 'NK-cell',
    'non-classical monocyte' = 'monocytes',
    'Non-switched memory B-cell' = 'B-cells',
    'Non-Vd2 gdTCR' = 'T-cells',
    'Plasmablast' = 'B-cells',
    'plasmacytoid DC' = 'dendritic cells',
    'Progenitor cell' = 'progenitor',
    'Switched memory B-cell' = 'B-cells',
    'T-reg' = 'T-cells',
    'Terminal effector memory CD4 T-cell' = 'T-cells',
    'Terminal effector memory CD8 T-cell' = 'T-cells',
    'total PBMC' = 'total PBMC',
    'Vd2 gdTCR' = 'T-cells',
    'Memory T-reg' = 'T-cells',
    'Naive CD4 T-cell activated' = 'T-cells',
    'Naive CD8 T-cell activated' = 'T-cells',
    'Naive T-reg' = 'T-cells',
    'eosinophil' = 'granulocytes',
    'gdTCR' = 'T-cells',
    'memory B-cell' = 'B-cells',
    'memory CD4 T-cell' = 'T-cells',
    'memory CD8 T-cell' = 'T-cells') %>%
  {tibble(content_name = names(.),
          lineage = .)}

  
