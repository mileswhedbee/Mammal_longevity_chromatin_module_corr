rm(list = ls())

library(dplyr)
library(purrr)
library(tibble)
library(Biobase)
library(readr)

gse_list <- getGEO(
  GEO = "GSE190756",
  GSEMatrix = TRUE,
  AnnotGPL = FALSE
)

all_geo_metadata <- purrr::imap_dfr(
  gse_list,
  function(eset, platform_file) {
    
    Biobase::pData(eset) %>%
      tibble::rownames_to_column("GSM") %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        geo_matrix_file = platform_file,
        geo_platform = Biobase::annotation(eset)
      )
  }
) %>%
  dplyr::distinct(GSM, .keep_all = TRUE) %>%
  dplyr::arrange(GSM)

dim(all_geo_metadata)

all_geo_metadata %>%
  dplyr::count(geo_platform, sort = TRUE)
