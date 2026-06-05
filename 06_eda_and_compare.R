# =====================================================================
# 06_eda_and_compare.R
# Reads output/features.csv and produces the exploratory plots and
# species comparisons the report needs, plus a small classifier used
# ONLY as a sanity check ("do these features separate species?").
#
# Run after 05_run_all.R. Saves figures to output/figs/.
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

FIG_DIR <- file.path(OUTPUT_DIR, "figs")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# load features, with a clear error if 05 hasn't run properly
.features_path <- file.path(OUTPUT_DIR, "features.csv")
if (!file.exists(.features_path)) {
  stop("features.csv not found at:\n  ", .features_path,
       "\nRun R/05_run_all.R first.", call. = FALSE)
}
feat <- readr::read_csv(.features_path, show_col_types = FALSE)
if (nrow(feat) == 0 || !"species" %in% names(feat)) {
  stop("features.csv exists but is empty or missing the 'species' column.\n",
       "This means R/05_run_all.R processed 0 images -- usually a path/template\n",
       "mismatch in R/config.R. Re-run R/05_run_all.R; its preflight check will\n",
       "tell you exactly which paths it tried and what is actually on disk.",
       call. = FALSE)
}
feat <- feat |>
  mutate(species = factor(species),
         species_name = factor(species_name))

# Distribution of a feature across species (histogram + ECDF) 
plot_feature_by_species <- function(df, var) {
  p_hist <- ggplot(df, aes(.data[[var]], fill = species_name)) +
    geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
    facet_wrap(~ species_name, ncol = 1, strip.position = "right") +
    labs(title = paste("Distribution of", var), x = var, y = "count") +
    theme_minimal() + theme(legend.position = "none")
  p_ecdf <- ggplot(df, aes(.data[[var]], color = species_name)) +
    stat_ecdf(linewidth = 0.8) +
    labs(title = paste("ECDF of", var), x = var, y = "F(x)", color = "species") +
    theme_minimal()
  p_hist | p_ecdf
}

# Save a few of the most informative outline features
for (v in c("aspect_ratio_blade", "circularity", "solidity",
            "eccentricity", "asymmetry_index")) {
  ggsave(file.path(FIG_DIR, paste0("dist_", v, ".png")),
         plot_feature_by_species(feat, v), width = 11, height = 7, dpi = 120)
}

# Petiole effect: raw vs blade aspect ratio
petiole_plot <- feat |>
  select(species_name, aspect_ratio_raw, aspect_ratio_blade) |>
  pivot_longer(-species_name, names_to = "type", values_to = "aspect") |>
  mutate(type = recode(type,
                       aspect_ratio_raw = "with petiole",
                       aspect_ratio_blade = "blade only")) |>
  ggplot(aes(species_name, aspect, fill = type)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  coord_flip() +
  labs(title = "Length/width ratio: effect of trimming the petiole",
       x = NULL, y = "aspect ratio", fill = NULL) +
  theme_minimal()
ggsave(file.path(FIG_DIR, "petiole_effect.png"), petiole_plot,
       width = 9, height = 6, dpi = 120)

# Color: hue & greenness by species 
color_plot <- (
  ggplot(feat, aes(species_name, H_mean, fill = species_name)) +
    geom_boxplot(alpha = 0.7) + coord_flip() +
    labs(title = "Mean hue by species", x = NULL, y = "H (0-1)") +
    theme_minimal() + theme(legend.position = "none")
) | (
  ggplot(feat, aes(species_name, green_frac_mean, fill = species_name)) +
    geom_boxplot(alpha = 0.7) + coord_flip() +
    labs(title = "Mean green fraction by species", x = NULL, y = "G/(R+G+B)") +
    theme_minimal() + theme(legend.position = "none")
)
ggsave(file.path(FIG_DIR, "color_by_species.png"), color_plot,
       width = 11, height = 5, dpi = 120)

# Vein density by species 
vein_plot <- ggplot(feat, aes(species_name, vein_len_per_area, fill = species_name)) +
  geom_boxplot(alpha = 0.7) + coord_flip() +
  labs(title = "Vein length per unit leaf area by species",
       x = NULL, y = "skeleton px / leaf px") +
  theme_minimal() + theme(legend.position = "none")
ggsave(file.path(FIG_DIR, "vein_by_species.png"), vein_plot,
       width = 9, height = 5, dpi = 120)

# Per-species summary table 
summary_table <- feat |>
  group_by(species_name) |>
  summarise(
    n = n(),
    area_mean = mean(area),
    aspect_blade = mean(aspect_ratio_blade),
    circularity = mean(circularity),
    solidity = mean(solidity),
    asymmetry = mean(asymmetry_index),
    hue = mean(H_mean),
    green = mean(green_frac_mean),
    vein = mean(vein_len_per_area),
    .groups = "drop"
  )
write.csv(summary_table, file.path(OUTPUT_DIR, "species_summary.csv"), row.names = FALSE)
print(summary_table)

message("EDA done. Figures in ", FIG_DIR)

