rm(list = ls())
library(dplyr)
library(stringr)
library(tidyr)
library(ape)
library(caper)
library(ggplot2)

read.csv(file = "data/GSE190756_family.soft") -> organism_list
read.delim(file = "data/mamms.species.txt", col.names = F) -> long_orgs
read.delim(file = "data/anage_data.txt") -> anage

colnames(organism_list) = "Organism"

# filter out rows == "!Platform_organism"
platform_rows <- organism_list %>%
  filter(str_starts(Organism, "!Platform_organism"))

df_split <- platform_rows %>%
  separate_wider_delim(
    Organism,
    delim = regex("\\s+"),
    names = c("field", "equals_sign", "Genus","Species")
  ) %>% select("Genus","Species")

paste(df_split$Genus, df_split$Species, sep = "_") -> df_split$Organism
df_split %>% select("Organism") -> df_split

colnames(long_orgs) = "Org"

long_orgs <- long_orgs %>%
  separate_wider_delim(
    Org,
    delim = regex("\\s+"),
    names = c("Genus","Species")
  ) 

paste(long_orgs$Genus, long_orgs$Species, sep = "_") -> long_orgs$Genus_species

long_orgs %>% select("Genus_species") -> long_orgs

# Join matching species together into new df
df1_intersection <- df_split %>%
  semi_join(long_orgs, by = c("Organism" = "Genus_species"))

paste(anage$Genus, anage$Species, sep = "_") -> anage$Species_lookup

anage <- anage %>%
  semi_join(df1_intersection, by = c("Species_lookup" = "Organism"))

#write.csv(anage, file = "output/Anage_mamms_for_analysis.csv", quote = F)

tree <- read.tree("data/Species_for_TreeTime_upload.nwk")

tree
length(tree$tip.label)
head(tree$tip.label, 20)

plot(tree, cex = 0.4)

# clean steps
tree$tip.label <- tree$tip.label |>
  gsub(" ", "_", x = _) |>
  gsub("'", "", x = _) |>
  gsub('"', "", x = _)

anage$Species_lookup <- anage$Species_lookup |>
  gsub(" ", "_", x = _) |>
  gsub("'", "", x = _) |>
  gsub('"', "", x = _)

# make dataframe
anage_pgls <- anage %>%
  transmute(
    Species_lookup,
    Adult_weight_g = suppressWarnings(as.numeric(Adult.weight..g.)),
    Maximum_longevity_yrs = suppressWarnings(
      as.numeric(Maximum.longevity..yrs.)
    )
  ) %>%
  filter(
    !is.na(Species_lookup),
    !is.na(Adult_weight_g),
    !is.na(Maximum_longevity_yrs),
    Adult_weight_g > 0,
    Maximum_longevity_yrs > 0
  ) %>%
  distinct(Species_lookup, .keep_all = TRUE) %>%
  mutate(
    log_adult_weight = log10(Adult_weight_g),
    log_max_longevity = log10(Maximum_longevity_yrs)
  )

shared_species <- intersect(
  anage_pgls$Species_lookup,
  tree$tip.label
)

length(shared_species)

setdiff(anage_pgls$Species_lookup, tree$tip.label) |> head(30)

setdiff(tree$tip.label, anage_pgls$Species_lookup) |> head(30)

anage_pgls_matched <- anage_pgls %>%
  filter(Species_lookup %in% shared_species)

tree_matched <- drop.tip(
  tree,
  setdiff(tree$tip.label, shared_species)
)

nrow(anage_pgls_matched)
length(tree_matched$tip.label)

all(anage_pgls_matched$Species_lookup %in% tree_matched$tip.label)

comparative_anage <- comparative.data(
  phy = tree_matched,
  data = anage_pgls_matched,
  names.col = "Species_lookup",
  vcv = TRUE,
  na.omit = TRUE,
  warn.dropped = TRUE
)

pgls_fit <- pgls(
  log_max_longevity ~ log_adult_weight,
  data = comparative_anage,
  lambda = "ML"
)

summary(pgls_fit)

# Coefficients: estimate, SE, t, and p value
summary(pgls_fit)$coefficients

# R-squared values
summary(pgls_fit)$r.squared

# Includes the maximum-likelihood estimate of Pagel's lambda
summary(pgls_fit)$param



pgls_style <- list(
  point_col = grDevices::adjustcolor("gray25", alpha.f = 0.55),
  point_pch = 16,
  point_cex = 1.2,
  line_col = "firebrick",
  line_lwd = 2.5,
  axis_col = "black",
  label_col = "black",
  bg_col = "white"
)

plot(
  anage_pgls_matched$log_adult_weight,
  anage_pgls_matched$log_max_longevity,
  pch = pgls_style$point_pch,
  cex = pgls_style$point_cex,
  col = pgls_style$point_col,
  xlab = expression(log[10]("Adult weight (g)")),
  ylab = expression(log[10]("Maximum longevity (years)")),
  las = 1,
  bty = "l",
  col.axis = pgls_style$axis_col,
  col.lab = pgls_style$label_col,
  bg = pgls_style$bg_col,
  main = "Mammals (N=19) for transcriptomic analysis"
)

abline(
  a = coef(pgls_fit)[1],
  b = coef(pgls_fit)[2],
  col = pgls_style$line_col,
  lwd = pgls_style$line_lwd
)
