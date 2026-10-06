rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
  library(tibble)
})

input_dir <- "data/GSE190756_processed"
output_dir <- "results/liver_chromatin_module"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

module_genes <- c(
  "Atrx",
  "Banf1",
  "Cbx1",
  "Cbx3",
  "Cbx5",
  "Dnmt1",
  "Dnmt3a",
  "Ehmt1",
  "Ehmt2",
  "Lmna",
  "Morc3",
  "Setdb1",
  "Suv39h1",
  "Suv39h2",
  "Trim28",
  "Uhrf1"
)

writeLines(
  module_genes,
  con = file.path(output_dir, "heterochromatin_core_genes.txt")
)

all_files <- list.files(
  path = input_dir,
  pattern = "\\.exp\\.salmon\\.txt\\.gz$",
  full.names = TRUE
)

liver_files <- all_files[
  str_detect(basename(all_files), regex("Liver", ignore_case = TRUE))
]

if (length(liver_files) == 0) {
  stop(
    "No liver files found. Check input_dir and make sure the expression files ",
    "end in .exp.salmon.txt.gz."
  )
}

message("Number of liver files found: ", length(liver_files))

get_gsm <- function(path) {
  gsm <- str_extract(basename(path), "GSM\\d+")
  
  if (is.na(gsm)) {
    stop("Could not extract GSM accession from: ", basename(path))
  }
  
  return(gsm)
}

read_expression_file <- function(path, module_genes) {
  gsm <- get_gsm(path)
  
  x <- read_tsv(
    file = path,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  
  if (ncol(x) < 2) {
    stop("File has fewer than two columns: ", basename(path))
  }
  
  first_value <- as.character(x[[1]][1])
  
  has_header <- first_value %in% c(
    "Name",
    "name",
    "Gene",
    "gene",
    "GeneID",
    "gene_id"
  )
  
  if (has_header) {
    header <- as.character(x[1, ])
    x <- x[-1, , drop = FALSE]
    colnames(x) <- header
    
    if ("TPM" %in% colnames(x)) {
      gene_values <- as.character(x[[1]])
      expr_values <- suppressWarnings(as.numeric(x[["TPM"]]))
    } else {
      gene_values <- as.character(x[[1]])
      expr_values <- suppressWarnings(as.numeric(x[[2]]))
    }
  } else {
    gene_values <- as.character(x[[1]])
    expr_values <- suppressWarnings(as.numeric(x[[2]]))
  }
  
  result <- tibble(
    GSM = gsm,
    filename = basename(path),
    gene = gene_values,
    expression = expr_values
  ) %>%
    filter(gene %in% module_genes)
  
  return(result)
}

module_long <- map_dfr(
  liver_files,
  ~ read_expression_file(.x, module_genes)
)

if (nrow(module_long) == 0) {
  stop(
    "No module genes were extracted. Check that gene symbols are in column 1."
  )
}

write_csv(
  module_long,
  file.path(output_dir, "liver_module_expression_long_by_sample.csv")
)

all_gsms <- sort(unique(vapply(liver_files, get_gsm, character(1))))

expected_sample_gene <- expand_grid(
  GSM = all_gsms,
  gene = module_genes
)

gene_counts <- module_long %>%
  count(GSM, gene, name = "n_rows")

detection_table <- expected_sample_gene %>%
  left_join(gene_counts, by = c("GSM", "gene")) %>%
  mutate(
    n_rows = replace_na(n_rows, 0L),
    detection_status = case_when(
      n_rows == 0 ~ "missing",
      n_rows == 1 ~ "present_once",
      n_rows > 1 ~ "duplicated"
    )
  ) %>%
  arrange(GSM, gene)

write_csv(
  detection_table,
  file.path(output_dir, "liver_module_gene_detection_by_sample.csv")
)

detection_summary <- detection_table %>%
  count(gene, detection_status, name = "n_samples") %>%
  complete(
    gene = module_genes,
    detection_status = c("missing", "present_once", "duplicated"),
    fill = list(n_samples = 0L)
  ) %>%
  pivot_wider(
    names_from = detection_status,
    values_from = n_samples,
    values_fill = 0
  ) %>%
  mutate(
    total_samples = missing + present_once + duplicated,
    prop_present_once = present_once / total_samples
  ) %>%
  arrange(desc(prop_present_once), gene)

write_csv(
  detection_summary,
  file.path(output_dir, "liver_module_gene_detection_summary.csv")
)

detection_problems <- detection_table %>%
  filter(detection_status != "present_once") %>%
  arrange(GSM, gene)

write_csv(
  detection_problems,
  file.path(output_dir, "liver_module_detection_problems.csv")
)

print(detection_summary)

if (any(detection_table$n_rows > 1)) {
  stop(
    "At least one gene occurs more than once in one or more samples. ",
    "Inspect liver_module_detection_problems.csv before continuing."
  )
}

module_wide_by_sample <- module_long %>%
  dplyr::distinct(GSM, filename, gene, .keep_all = TRUE) %>%
  dplyr::select(GSM, filename, gene, expression) %>%
  tidyr::pivot_wider(
    id_cols = c(GSM, filename),
    names_from = gene,
    values_from = expression
  ) %>%
  dplyr::arrange(GSM)

write_csv(
  module_wide_by_sample,
  file.path(output_dir, "liver_module_expression_wide_by_sample.csv")
)

expression_summary <- module_long %>%
  group_by(gene) %>%
  summarise(
    n_samples = n(),
    n_nonmissing = sum(!is.na(expression)),
    n_positive = sum(expression > 0, na.rm = TRUE),
    min_expression = min(expression, na.rm = TRUE),
    median_expression = median(expression, na.rm = TRUE),
    mean_expression = mean(expression, na.rm = TRUE),
    max_expression = max(expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene)

write_csv(
  expression_summary,
  file.path(output_dir, "liver_module_expression_summary_by_gene.csv")
)

print(expression_summary)

metadata_template <- tibble(
  GSM = all_gsms,
  species = NA_character_,
  tissue = "Liver",
  replicate = NA_character_,
  sex = NA_character_,
  age = NA_character_,
  notes = NA_character_
)

template_file <- file.path(
  output_dir,
  "GSE190756_liver_sample_metadata_TEMPLATE.csv"
)

if (!file.exists(template_file)) {
  write_csv(metadata_template, template_file)
}

metadata_file <- file.path(
  output_dir,
  "GSE190756_liver_sample_metadata.csv"
)

if (!file.exists(metadata_file)) {
  message(
    "\nStage 1 completed successfully.\n\n",
    "Inspect these output files:\n",
    "  liver_module_gene_detection_summary.csv\n",
    "  liver_module_detection_problems.csv\n",
    "  liver_module_expression_wide_by_sample.csv\n\n",
    "Next, fill in the species column in:\n  ",
    template_file,
    "\n\nThen save it as:\n  ",
    metadata_file,
    "\n\nRe-run this script to create species-level expression values and ",
    "chromatin module scores."
  )
} else {
  sample_metadata <- read_csv(
    metadata_file,
    show_col_types = FALSE
  ) %>%
    mutate(GSM = as.character(GSM))
  
  if (!all(c("GSM", "species") %in% names(sample_metadata))) {
    stop(
      "The metadata file must contain columns named GSM and species."
    )
  }
  
  missing_species <- sample_metadata %>%
    filter(is.na(species) | species == "")
  
  if (nrow(missing_species) > 0) {
    stop(
      "Some samples have missing species values. Fill in every species value ",
      "in the metadata file before rerunning."
    )
  }
  
  module_with_metadata <- module_long %>%
    left_join(sample_metadata, by = "GSM")
  
  unmatched_samples <- module_with_metadata %>%
    filter(is.na(species)) %>%
    distinct(GSM)
  
  if (nrow(unmatched_samples) > 0) {
    stop(
      "These samples were not found in the metadata file:\n",
      paste(unmatched_samples$GSM, collapse = "\n")
    )
  }
  
  module_species_long <- module_with_metadata %>%
    mutate(log2_expression = log2(expression + 1)) %>%
    group_by(species, gene) %>%
    summarise(
      n_replicates = sum(!is.na(log2_expression)),
      mean_log2_expression = mean(log2_expression, na.rm = TRUE),
      sd_log2_expression = sd(log2_expression, na.rm = TRUE),
      .groups = "drop"
    )
  
  write_csv(
    module_species_long,
    file.path(output_dir, "liver_module_expression_long_by_species.csv")
  )
  
  module_species_z <- module_species_long %>%
    group_by(gene) %>%
    mutate(
      gene_mean = mean(mean_log2_expression, na.rm = TRUE),
      gene_sd = sd(mean_log2_expression, na.rm = TRUE),
      gene_z = if_else(
        is.na(gene_sd) | gene_sd == 0,
        NA_real_,
        (mean_log2_expression - gene_mean) / gene_sd
      )
    ) %>%
    ungroup()
  
  write_csv(
    module_species_z,
    file.path(output_dir, "liver_module_expression_zscores_by_species.csv")
  )
  
  module_score <- module_species_z %>%
    group_by(species) %>%
    summarise(
      n_module_genes = sum(!is.na(gene_z)),
      chromatin_module_score = mean(gene_z, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(chromatin_module_score))
  
  write_csv(
    module_score,
    file.path(output_dir, "liver_chromatin_module_score_by_species.csv")
  )
  
  module_species_wide <- module_species_long %>%
    dplyr::select(species, gene, mean_log2_expression) %>%
    pivot_wider(
      names_from = gene,
      values_from = mean_log2_expression
    ) %>%
    arrange(species)
  
  write_csv(
    module_species_wide,
    file.path(output_dir, "liver_module_expression_wide_by_species.csv")
  )
  
  saveRDS(
    list(
      module_genes = module_genes,
      module_long = module_long,
      detection_table = detection_table,
      detection_summary = detection_summary,
      sample_metadata = sample_metadata,
      module_species_long = module_species_long,
      module_species_z = module_species_z,
      module_score = module_score
    ),
    file.path(output_dir, "liver_chromatin_module_analysis_objects.rds")
  )
  
  print(module_score)
  
  message(
    "\nStage 2 completed successfully.\n\n",
    "Your species-level score file is:\n  ",
    file.path(output_dir, "liver_chromatin_module_score_by_species.csv"),
    "\n\nThe next step is to join this file to your AnAge adult-mass and ",
    "maximum-longevity data, match species to the mammalian phylogeny, and fit PGLS."
  )
}

