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
                        "not expressed" = "grey"  )

protein.class.palette <- c("secreted" = '#911D51',
                           "membrane" = '#6D4BAA', 
                           "other" = '#008490', 
                           "cd_marker" = '#318F1E', 
                           "transcription_factors" = '#B8801B', 
                           "mitochondrial" = '#E371B4', 
                           "ribosomal" = '#89A0F3', 
                           '#00C9BC', '#97C542', '#FFA05E')
