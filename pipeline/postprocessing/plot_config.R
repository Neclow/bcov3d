library(ape)
library(dplyr)
library(stringr)

source("src/config.R")

mast_analysis <- function(mast_path, criterion = "BIC") {
  mast_scores <- read.csv(paste0(mast_path, ".csv"), row.names = 1)
  mast_scores <- mast_scores[order(as.numeric(rownames(mast_scores))), ]
  argbest <- rownames(mast_scores)[which.min(mast_scores[[criterion]])]
  return(argbest)
}

virusinfo <- read.csv(file.path(metadir, "receptors.csv"))
virusinfo <- virusinfo %>%
  mutate(
    Variant = ifelse(grepl("\\(", Virus),
      str_extract(Virus, "\\(.*\\)"),
      ""
    ),
    Virus = str_replace(Virus, "\\s*\\(.*\\)", ""),
    Virus = str_replace_all(Virus, " ", "")
  ) %>%
  mutate(Variant = str_replace_all(Variant, "[\\(\\)]", ""))

head(virusinfo)

argbest_aa <- mast_analysis(mast_path_aa, "BIC")
argbest_3di <- mast_analysis(mast_path_3di, "BIC")

trs_aa <- read.tree(
  file = file.path(
    mast_path_aa, paste0(argbest_aa, ".treefile.annotated")
  )
)

trs_3di <- read.tree(
  file = file.path(
    mast_path_3di, paste0(argbest_3di, ".treefile.annotated")
  )
)

for (i in seq_along(trs_aa)) {
  for (j in seq_along(trs_aa[[i]]$tip.label)) {
    label_aa <- trs_aa[[i]]$tip.label[j]
    label_aa <- gsub("'", "", label_aa)
    label_aa <- gsub("Cov", "CoV", label_aa)
    label_aa <- gsub("Bat CoV", "BatCoV", label_aa)
    trs_aa[[i]]$tip.label[j] <- label_aa
    label_3di <- trs_3di[[i]]$tip.label[j]
    label_3di <- gsub("'", "", label_3di)
    label_3di <- gsub("Cov", "CoV", label_3di)
    label_3di <- gsub("Bat CoV", "BatCoV", label_3di)
    trs_3di[[i]]$tip.label[j] <- label_3di
  }
}

# NTD
ntd_class_aa <- 1
ntd_class_3di <- 2
# RBD
rbd_class_aa <- 3
rbd_class_3di <- 1
# S2
s2_class_aa <- 2
s2_class_3di <- 3

ntd_tree_aa <- trs_aa[[ntd_class_aa]]
ntd_tree_3di <- trs_3di[[ntd_class_3di]]
rbd_tree_aa <- trs_aa[[rbd_class_aa]]
rbd_tree_3di <- trs_3di[[rbd_class_3di]]
s2_tree_aa <- trs_aa[[s2_class_aa]]
s2_tree_3di <- trs_3di[[s2_class_3di]]
