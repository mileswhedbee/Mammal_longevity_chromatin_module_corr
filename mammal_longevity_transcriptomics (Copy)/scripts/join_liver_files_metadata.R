library(dplyr)
library(stringr)
library(tibble)
library(readr)

metadata_objects <- readRDS(
  "metadata/GSE190756_metadata_objects.rds"
)

all_geo_metadata <- metadata_objects$all_geo_metadata
liver_file_map <- metadata_objects$liver_file_map
liver_geo_metadata <- metadata_objects$liver_geo_metadata

input_dir <- "data/GSE190756_processed"

getwd()

dir.create("metadata", recursive = TRUE, showWarnings = FALSE)

dir.exists("metadata")
file.access("metadata", mode = 2)

all_expression_files <- list.files(
  path = input_dir,
  pattern = "\\.exp\\.salmon\\.txt\\.gz$",
  full.names = TRUE
)

liver_expression_files <- all_expression_files[
  stringr::str_detect(
    basename(all_expression_files),
    stringr::regex("Liver", ignore_case = TRUE)
  )
]

liver_file_map <- tibble::tibble(
  GSM = stringr::str_extract(
    basename(liver_expression_files),
    "GSM\\d+"
  ),
  filename = basename(liver_expression_files)
) %>%
  dplyr::arrange(GSM)

nrow(liver_file_map)

liver_geo_metadata <- liver_file_map %>%
  dplyr::left_join(
    all_geo_metadata,
    by = "GSM"
  )

readr::write_csv(
  liver_geo_metadata,
  "metadata/GSE190756_liver_GEO_metadata_joined.csv"
)

liver_geo_metadata %>%
  dplyr::filter(is.na(geo_platform)) %>%
  dplyr::select(GSM, filename)

nrow(liver_geo_metadata)

liver_geo_metadata %>%
  dplyr::count(
    geo_platform,
    organism_ch1,
    sort = TRUE
  ) %>%
  print(n = Inf, width = Inf)

saveRDS(
  list(
    all_geo_metadata = all_geo_metadata,
    liver_file_map = liver_file_map,
    liver_geo_metadata = liver_geo_metadata
  ),
  "metadata/GSE190756_metadata_objects.rds"
)

library(dplyr)
library(readr)

liver_analysis_metadata <- liver_geo_metadata %>%
  dplyr::transmute(
    GSM,
    filename,
    species = organism_ch1,
    tissue = source_name_ch1,
    replicate = NA_character_,
    sex = NA_character_,
    age = NA_character_,
    notes = dplyr::coalesce(title, source_name_ch1),
    geo_platform
  ) %>%
  dplyr::arrange(species, GSM)

readr::write_csv(
  liver_analysis_metadata,
  "results/liver_chromatin_module/GSE190756_liver_sample_metadata.csv"
)

