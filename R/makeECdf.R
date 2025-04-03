library(readxl)
library(dplyr)
library(stringr)

source('R/convertGenes_updateGenes.R')

organs <- c("brain", "testis", "liver", "spleen", "lung", "kidney", "colon", "small_intestine", "heart", "EDL", "soleus", "pancreas", "organ")

is_any_capitalized <- function(vec) {
  any(grepl("^[A-Z]+$", vec))
}

all_genes <- c()

ec_df <- tibble(EC = character(),
                Genes = character())

for (i in 1:13) {

  markers <- read_excel('rutils/kalucka_murine_ecatlas_allvascularmarkers.xlsx', sheet = i)

  for (j in 1:ncol(markers)) {
    genes <- as.vector(na.omit(markers[,j][[1]]))

    if(i == 12 | (i == 13 & j == 11)) {

      genes <- convertGenes(genes, convert_to = 'mouse')
    }

    all_genes <- unique(c(all_genes,genes))

    name <- gsub(" ", "_", names(markers[,j]))

    name <- str_to_lower(paste0(organs[i], "_", name))

    df = tibble(rep(name, length(genes)), genes)

    ec_df = rbind(ec_df, df)
  }

}

colnames(ec_df) <- c("label", "genes")

ec_df <- ec_df %>%
  mutate(vascular_subtype = case_when(
    str_detect(label, "large") & str_detect(label, "artery") ~ "Large Artery",
    str_detect(label, "artery") & !str_detect(label, "large") ~ "Artery",
    str_detect(label, "arterial") | str_detect(label, "arteriole") ~ "Arteriole",
    str_detect(label, "vein") ~ "Vein",
    str_detect(label, "venous") | str_detect(label, "venule") ~ "Venule",
    str_detect(label, "proliferating") ~ "Proliferating",
    str_detect(label, "capillary") & !str_detect(label, "arterial") &
      !str_detect(label, "venous") & !str_detect(label, "arteriole") &
      !str_detect(label, "venule") ~ "Capillary",
    str_detect(label, "lymphatic") | str_detect(label, "lypmhatic") ~ "Lymphatic",
    str_detect(label, "angiogenic") ~ "Angiogenic",
    str_detect(label, "choroidplexus") | str_detect(label, "glomeruli") |
      str_detect(label, "colon_...7") ~ "Capillary",
    str_detect(label, "interferon") ~ "Immunogenic",
  ))


write.csv(ec_df, file = 'rutils/ec_vascular_subtype_signatures_dataframe.csv')
