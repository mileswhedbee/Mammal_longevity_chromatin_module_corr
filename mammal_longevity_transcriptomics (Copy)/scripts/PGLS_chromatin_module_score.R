rm(list = ls())

module_score_file <- "results/liver_chromatin_module/liver_chromatin_module_score_by_species.csv"

anage_file <- "data/Anage_mamms_for_analysis.csv"

tree_file <- "data/Species_for_TreeTime_upload.nwk"

library(readr)
library(dplyr)
library(stringr)
library(ape)

module_scores <- readr::read_csv(
  module_score_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    species = stringr::str_squish(species)
  )

anage_raw <- readr::read_csv(
  anage_file,
  show_col_types = FALSE
)

mammal_tree <- ape::read.tree(tree_file)

names(module_scores)
names(anage_raw)
glimpse(module_scores)
glimpse(anage_raw)

length(mammal_tree$tip.label)
mammal_tree$tip.label
is.null(mammal_tree$edge.length)
summary(mammal_tree$edge.length)

mammal_tree$tip.label <- mammal_tree$tip.label %>%
  stringr::str_replace_all("_", " ") %>%
  stringr::str_squish()

anage_traits <- anage_raw %>%
  dplyr::transmute(
    anage_species = stringr::str_squish(
      paste(Genus, Species)
    ),
    
    # Replace these two fields with the exact column names in anage_raw:
    adult_mass_g = as.numeric(`Adult.weight..g.`),
    max_longevity_y = as.numeric(`Maximum.longevity..yrs.`)
  ) %>%
  dplyr::mutate(
    anage_species = stringr::str_replace_all(anage_species, "_", " "),
    anage_species = stringr::str_squish(anage_species)
  ) %>%
  dplyr::filter(
    !is.na(anage_species),
    anage_species != "",
    is.finite(adult_mass_g),
    adult_mass_g > 0,
    is.finite(max_longevity_y),
    max_longevity_y > 0
  ) %>%
  dplyr::distinct(anage_species, .keep_all = TRUE)

anage_traits <- anage_traits %>%
  dplyr::mutate(
    species_key = normalize_species(anage_species)
  )

normalize_species <- function(x) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_lower()
}


tree_taxa <- tibble::tibble(
  tree_species = mammal_tree$tip.label
) %>%
  dplyr::mutate(
    tree_species = stringr::str_replace_all(tree_species, "_", " "),
    tree_species = stringr::str_squish(tree_species),
    species_key = normalize_species(tree_species)
  )

#####################################3
# PGLS join
########################################

module_scores <- module_scores %>%
  dplyr::mutate(
    species_key = normalize_species(.data[["species"]])
  )

anage_traits <- anage_traits %>%
  dplyr::mutate(
    species_key = normalize_species(.data[["anage_species"]])
  )


pgls_data <- tree_taxa %>%
  dplyr::left_join(
    module_scores %>%
      dplyr::select(
        species_key,
        expression_species = species,
        chromatin_module_score,
        n_module_genes
      ),
    by = "species_key"
  ) %>%
  dplyr::left_join(
    anage_traits %>%
      dplyr::select(
        species_key,
        anage_species,
        adult_mass_g,
        max_longevity_y
      ),
    by = "species_key"
  ) %>%
  dplyr::mutate(
    log10_adult_mass_g = log10(adult_mass_g),
    log10_max_longevity_y = log10(max_longevity_y)
  )


if (!requireNamespace("caper", quietly = TRUE)) {
  install.packages("caper")
}

################# run PGLS
library(caper)

comp_data <- caper::comparative.data(
  phy = mammal_tree,
  data = as.data.frame(pgls_data),
  names.col = "tree_species",
  vcv = TRUE,
  vcv.dim = 2,
  na.omit = FALSE,
  warn.dropped = TRUE
)

comp_data$dropped

pgls_primary <- caper::pgls(
  log10_max_longevity_y ~
    log10_adult_mass_g +
    chromatin_module_score,
  data = comp_data,
  lambda = "ML"
)

summary(pgls_primary)

pgls_mass_only <- caper::pgls(
  log10_max_longevity_y ~ log10_adult_mass_g,
  data = comp_data,
  lambda = "ML"
)

AIC(pgls_mass_only, pgls_primary)


################################
# drop tree tips one by one and rerun analysis

loo_results <- purrr::map_dfr(
  pgls_data$tree_species,
  function(drop_species) {
    dat_i <- pgls_data %>%
      dplyr::filter(tree_species != drop_species)
    
    tree_i <- ape::drop.tip(mammal_tree, drop_species)
    
    comp_i <- caper::comparative.data(
      phy = tree_i,
      data = as.data.frame(dat_i),
      names.col = "tree_species",
      vcv = TRUE,
      vcv.dim = 2,
      na.omit = FALSE,
      warn.dropped = FALSE
    )
    
    fit_i <- caper::pgls(
      log10_max_longevity_y ~
        log10_adult_mass_g +
        chromatin_module_score,
      data = comp_i,
      lambda = "ML"
    )
    
    coefs_i <- summary(fit_i)$coefficients
    
    tibble::tibble(
      excluded_species = drop_species,
      n_species = nrow(dat_i),
      lambda_ml = unname(fit_i$param["lambda"]),
      module_beta = coefs_i["chromatin_module_score", "Estimate"],
      module_se = coefs_i["chromatin_module_score", "Std. Error"],
      module_p = coefs_i["chromatin_module_score", "Pr(>|t|)"],
      AIC = AIC(fit_i)
    )
  }
)

loo_results %>%
  dplyr::arrange(module_beta) %>%
  print(n = Inf, width = Inf)

readr::write_csv(
  loo_results,
  "results/liver_chromatin_module/pgls_19_species_leave_one_out.csv"
)


##########################################################
# Check Sciurus carolinensis, is it real or technical?
###########################################################

read_expression_file <- function(file_path) {
  x <- readr::read_tsv(
    file_path,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  if (ncol(x) < 2) {
    stop(
      "Expected at least 2 columns in: ", file_path,
      "\nFound ", ncol(x), " column(s)."
    )
  }
  
  gene_col <- names(x)[1]
  expression_col <- names(x)[2]
  
  x %>%
    dplyr::transmute(
      gene = as.character(.data[[gene_col]]),
      expression = suppressWarnings(
        as.numeric(.data[[expression_col]])
      )
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      is.finite(expression),
      expression >= 0
    )
}

test_file <- expression_manifest$file_path[1]

test_expression <- read_expression_file(test_file)

test_expression %>%
  dplyr::slice_head(n = 10) %>%
  print(n = Inf, width = Inf)

test_expression %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    min_expression = min(expression),
    median_expression = median(expression),
    max_expression = max(expression),
    n_zero = sum(expression == 0)
  )


sample_global_qc <- expression_manifest %>%
  dplyr::mutate(
    global_qc = purrr::map(
      file_path,
      function(path) {
        dat <- read_expression_file(path) %>%
          dplyr::mutate(
            log2_expression = log2(expression + 1)
          )
        
        tibble::tibble(
          n_quantified_genes = nrow(dat),
          n_expressed_genes = sum(dat$expression > 0, na.rm = TRUE),
          fraction_expressed = mean(dat$expression > 0, na.rm = TRUE),
          median_log2_expression = median(
            dat$log2_expression,
            na.rm = TRUE
          ),
          mean_log2_expression = mean(
            dat$log2_expression,
            na.rm = TRUE
          ),
          p25_log2_expression = unname(
            stats::quantile(
              dat$log2_expression,
              0.25,
              na.rm = TRUE
            )
          ),
          p75_log2_expression = unname(
            stats::quantile(
              dat$log2_expression,
              0.75,
              na.rm = TRUE
            )
          )
        )
      }
    )
  ) %>%
  dplyr::select(-file_path) %>%
  tidyr::unnest(global_qc)

species_global_qc <- sample_global_qc %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    median_of_sample_medians = median(
      median_log2_expression,
      na.rm = TRUE
    ),
    mean_of_sample_medians = mean(
      median_log2_expression,
      na.rm = TRUE
    ),
    sd_sample_medians = sd(
      median_log2_expression,
      na.rm = TRUE
    ),
    mean_expressed_genes = mean(
      n_expressed_genes,
      na.rm = TRUE
    ),
    mean_fraction_expressed = mean(
      fraction_expressed,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    global_median_z = as.numeric(
      scale(median_of_sample_medians)
    )
  ) %>%
  dplyr::arrange(median_of_sample_medians)

species_global_qc %>%
  dplyr::select(
    species,
    n_samples,
    median_of_sample_medians,
    global_median_z,
    mean_expressed_genes,
    mean_fraction_expressed
  ) %>%
  print(n = Inf, width = Inf)


library(dplyr)
library(ggplot2)
install.packages("ggrepel")
library(ggrepel)
library(caper)

# PGLS used to remove adult-mass effects from longevity
pgls_longevity_mass <- caper::pgls(
  log10_max_longevity_y ~ log10_adult_mass_g,
  data = comp_data,
  lambda = "ML"
)

# PGLS used to remove adult-mass effects from module score
pgls_module_mass <- caper::pgls(
  chromatin_module_score ~ log10_adult_mass_g,
  data = comp_data,
  lambda = "ML"
)

# Attach residuals in the same order as pgls_data
plot_data <- pgls_data %>%
  dplyr::mutate(
    longevity_mass_residual = residuals(pgls_longevity_mass),
    module_mass_residual = residuals(pgls_module_mass),
    label_species = dplyr::if_else(
      tree_species == "Sciurus carolinensis",
      tree_species,
      NA_character_
    )
  )

# A simple visual line through the residualized data.
# Do not use its OLS p-value as the PGLS inference.
plot_line <- stats::lm(
  longevity_mass_residual ~ module_mass_residual,
  data = plot_data
)

p_module <- summary(pgls_primary)$coefficients[
  "chromatin_module_score",
  "Pr(>|t|)"
]

beta_module <- summary(pgls_primary)$coefficients[
  "chromatin_module_score",
  "Estimate"
]

lambda_primary <- unname(pgls_primary$param["lambda"])

p_partial <- ggplot(
  plot_data,
  aes(
    x = module_mass_residual,
    y = longevity_mass_residual
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "grey75"
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "grey75"
  ) +
  geom_point(
    aes(
      colour = tree_species == "Sciurus carolinensis",
      shape = tree_species == "Sciurus carolinensis"
    ),
    size = 3.2,
    alpha = 0.9
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.8,
    colour = "#1F4E79",
    fill = "#9ECAE1"
  ) +
  ggrepel::geom_text_repel(
    aes(label = label_species),
    na.rm = TRUE,
    size = 3.6,
    min.segment.length = 0,
    seed = 1
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "#4D4D4D",
      "TRUE" = "#D55E00"
    ),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c(
      "FALSE" = 16,
      "TRUE" = 17
    ),
    guide = "none"
  ) +
  labs(
    x = "Chromatin-module score residual\n(after PGLS adjustment for adult mass)",
    y = "log10(maximum longevity) residual\n(after PGLS adjustment for adult mass)",
    title = "Mass-adjusted association of liver chromatin-module expression with longevity",
    subtitle = paste0(
      "PGLS: βmodule = ",
      sprintf("%.3f", beta_module),
      "; p = ",
      sprintf("%.3f", p_module),
      "; Pagel's λ = ",
      sprintf("%.2f", lambda_primary),
      "; n = ",
      nrow(plot_data)
    ),
    caption = paste(
      "Each point is one species. The fitted line and ribbon are visual aids;",
      "formal inference comes from the PGLS model."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(colour = "grey25"),
    plot.caption = element_text(
      hjust = 0,
      colour = "grey35",
      size = 9
    )
  )

p_partial
