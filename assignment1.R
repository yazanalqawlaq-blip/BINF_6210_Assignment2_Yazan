# BINF6210 — Assignment 1 (Storyboard)
# Yazan Alqawlaq
# Date: 2025-10-13
# Question: Do Bombus BIN communities show comparable distance–decay and β-partition structure within the Nearctic vs Western Palearctic?
# Project type: Hypothesis-testing
# Hypothesis: Within each region, BIN dissimilarity increases with geographic distance, and turnover dominates β-diversity (nestedness is minor).

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan) # distances, Mantel, betadisper, adonis2
  library(betapart) # beta.pair (Jaccard partition)
  library(geosphere) # distHaversine
})

theme_set(theme_minimal(base_size = 12))
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold"),
    strip.background = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(6, 6, 6, 6)
  )

# directory
fig_dir <- "../figs"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# 1) Data
data_path <- "../data/result.tsv"
dfBOLD <- readr::read_tsv(data_path, show_col_types = FALSE)

dfBOLD <- dfBOLD %>%
  rename(country = `country/ocean`, province_state = `province/state`) %>%
  mutate(
    province_state = na_if(province_state, ""),
    site_id = paste(country_iso, province_state, sep = ":")
  ) %>%
  filter(!is.na(province_state), !country %in% c("NA", "Unrecoverable"))

# setting regions
realm_by_site <- dfBOLD %>%
  filter(!is.na(realm)) %>%
  count(site_id, realm, name = "n_realm") %>%
  group_by(site_id) %>%
  slice_max(n_realm, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(region = case_when(
    str_detect(realm, regex("Nearctic", TRUE)) ~ "Nearctic",
    str_detect(realm, regex("Palearctic", TRUE)) ~ "Western Palearctic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(region)) %>%
  select(site_id, region)

# 2) balancing sites
build_with_threshold <- function(bin_threshold = 10) {
  site_stats <- dfBOLD %>%
    filter(!is.na(bin_uri)) %>%
    summarise(n_BINs = n_distinct(bin_uri), records = dplyr::n(), .by = site_id)

  site_keep0 <- site_stats %>%
    inner_join(realm_by_site, by = "site_id") %>%
    filter(n_BINs >= bin_threshold) %>%
    select(site_id, region, n_BINs, records)

  # median calculation
  coords_candidates <- dfBOLD %>%
    filter(site_id %in% site_keep0$site_id, !is.na(coord)) %>%
    mutate(coord_clean = stringr::str_squish(stringr::str_replace_all(coord, "[\\[\\]()]", ""))) %>%
    tidyr::separate_wider_delim(coord_clean, ",",
      names = c("lat_chr", "lon_chr"),
      too_many = "drop", too_few = "align_start"
    ) %>%
    mutate(
      lat_num = readr::parse_number(lat_chr),
      lon_num = readr::parse_number(lon_chr),
      lat_num = case_when(
        str_detect(lat_chr, "[Ss]") ~ -abs(lat_num),
        str_detect(lat_chr, "[Nn]") ~ abs(lat_num),
        TRUE ~ lat_num
      ),
      lon_num = case_when(
        str_detect(lon_chr, "[Ww]") ~ -abs(lon_num),
        str_detect(lon_chr, "[Ee]") ~ abs(lon_num),
        TRUE ~ lon_num
      )
    ) %>%
    filter(
      !is.na(lat_num), !is.na(lon_num),
      dplyr::between(lat_num, -90, 90),
      dplyr::between(lon_num, -180, 180),
      !(lat_num == 0 & lon_num == 0)
    ) %>%
    summarise(
      lat = median(lat_num, na.rm = TRUE),
      lon = median(lon_num, na.rm = TRUE),
      .by = site_id
    )

  site_keep <- site_keep0 %>% inner_join(coords_candidates, by = "site_id")

  sites_NA <- site_keep %>%
    filter(region == "Nearctic") %>%
    pull(site_id)
  sites_EU <- site_keep %>%
    filter(region == "Western Palearctic") %>%
    pull(site_id)
  S <- min(length(sites_NA), length(sites_EU))
  if (is.na(S) || S < 3) {
    return(NULL)
  }

  set.seed(20251012)
  pick_NA <- sample(sites_NA, S)
  pick_EU <- sample(sites_EU, S)
  sites_sel <- c(pick_NA, pick_EU)

  site_meta <- site_keep %>%
    filter(site_id %in% sites_sel) %>%
    select(site_id, region, n_BINs, lat, lon)

  df_bin_site <- dfBOLD %>%
    filter(site_id %in% sites_sel, !is.na(bin_uri)) %>%
    count(site_id, bin_uri, name = "n")

  comm_wide <- df_bin_site %>%
    tidyr::pivot_wider(names_from = bin_uri, values_from = n, values_fill = 0)

  if (ncol(comm_wide) < 3) {
    return(NULL)
  }

  comm_mat <- comm_wide %>%
    column_to_rownames("site_id") %>%
    as.matrix()
  comm_pa <- ifelse(comm_mat > 0, 1L, 0L)

  # align data
  site_meta <- site_meta %>% filter(site_id %in% rownames(comm_pa))
  coord_order <- site_meta %>% arrange(match(site_id, rownames(comm_pa)))
  stopifnot(identical(coord_order$site_id, rownames(comm_pa)))

  # pairwise geo distance (km)
  coords_mat <- as.matrix(coord_order[, c("lon", "lat")])
  geo_mat <- geosphere::distm(coords_mat, fun = geosphere::distHaversine) / 1000
  rownames(geo_mat) <- colnames(geo_mat) <- rownames(comm_pa)

  # Jaccard dissimilarities
  jac_dist <- vegan::vegdist(comm_pa, method = "jaccard", binary = TRUE)
  jac_mat <- as.matrix(jac_dist)

  list(
    S = S, site_meta = site_meta, comm_pa = comm_pa,
    jac_dist = jac_dist, jac_mat = jac_mat, geo_mat = geo_mat
  )
}

built <- build_with_threshold(10)
if (is.null(built)) {
  message("Not enough balanced sites at ≥10 BINs; retrying with ≥5 BINs.")
  built <- build_with_threshold(5)
  if (is.null(built)) stop("Still not enough balanced sites at ≥5 BINs.")
}

S <- built$S
site_meta <- built$site_meta
comm_pa <- built$comm_pa
jac_dist <- built$jac_dist
jac_mat <- built$jac_mat
geo_mat <- built$geo_mat

# 3) Figure 1 — Distance–decay within each region
sites_NA_sel <- site_meta$site_id[site_meta$region == "Nearctic"]
sites_EU_sel <- site_meta$site_id[site_meta$region == "Western Palearctic"]

jac_na <- jac_mat[sites_NA_sel, sites_NA_sel, drop = FALSE]
geo_na <- geo_mat[sites_NA_sel, sites_NA_sel, drop = FALSE]
jac_eu <- jac_mat[sites_EU_sel, sites_EU_sel, drop = FALSE]
geo_eu <- geo_mat[sites_EU_sel, sites_EU_sel, drop = FALSE]

mantel_na <- vegan::mantel(as.dist(jac_na), as.dist(geo_na), method = "spearman", permutations = 999)
mantel_eu <- vegan::mantel(as.dist(jac_eu), as.dist(geo_eu), method = "spearman", permutations = 999)

dd_df <- bind_rows(
  tibble(
    region = "Nearctic",
    distance_km = as.numeric(as.dist(geo_na)),
    jaccard = as.numeric(as.dist(jac_na))
  ),
  tibble(
    region = "Western Palearctic",
    distance_km = as.numeric(as.dist(geo_eu)),
    jaccard = as.numeric(as.dist(jac_eu))
  )
) %>% tidyr::drop_na(distance_km, jaccard)

ggplot(dd_df, aes(distance_km, jaccard)) +
  geom_point(alpha = 0.35, size = 1.6) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~region, ncol = 2, scales = "free_x") +
  scale_x_continuous(labels = scales::label_number(big.mark = ",")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Distance–decay of Bombus BIN composition",
    x = "Distance (km)",
    y = "Jaccard dissimilarity"
  ) +
  theme_paper

ggsave(file.path(fig_dir, "A_distance_decay.png"),
  
  width = 8, height = 4.8, dpi = 300
)

# Reference (outside course content):
# vegan::mantel — https://search.r-project.org/CRAN/refmans/vegan/html/mantel.html
# geosphere::distHaversine — https://search.r-project.org/CRAN/refmans/geosphere/html/distHaversine.html

# 4) Figure 2 — Jaccard
bp <- betapart::beta.pair(comm_pa, index.family = "jaccard")
b_total <- as.matrix(bp$beta.jac) 
b_turn <- as.matrix(bp$beta.jtu) 
b_nest <- as.matrix(bp$beta.jne) 

pairs_vec <- function(mat, sites) as.numeric(as.dist(mat[sites, sites, drop = FALSE]))

df_beta <- bind_rows(
  tibble(
    region = "Nearctic",
    total = pairs_vec(b_total, sites_NA_sel),
    turn = pairs_vec(b_turn, sites_NA_sel),
    nest = pairs_vec(b_nest, sites_NA_sel)
  ),
  tibble(
    region = "Western Palearctic",
    total = pairs_vec(b_total, sites_EU_sel),
    turn = pairs_vec(b_turn, sites_EU_sel),
    nest = pairs_vec(b_nest, sites_EU_sel)
  )
) %>%
  tidyr::pivot_longer(c(total, turn, nest), names_to = "component", values_to = "value") %>%
  mutate(
    component = factor(
      recode(component,
        total = "Total β (Jaccard)",
        turn  = "Turnover",
        nest  = "Nestedness"
      ),
      levels = c("Total β (Jaccard)", "Turnover", "Nestedness")
    ),
    region = factor(region, levels = c("Nearctic", "Western Palearctic"))
  )

ggplot(df_beta, aes(region, value)) +
  geom_violin(fill = "grey95", color = "grey70", scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.16, outlier.shape = NA) +
  facet_wrap(~component, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "β-diversity partition within regions",
    x = NULL, y = "Dissimilarity"
  ) +
  theme_paper

ggsave(file.path(fig_dir, "B_beta_partition.png"),
  
  width = 10, height = 4, dpi = 300
)

# Reference (outside course content):
# betapart::beta.pair — https://search.r-project.org/CRAN/refmans/betapart/html/beta.pair.html

# Figure 3: Unique vs shared BINs 
stopifnot(exists("comm_pa"), exists("site_meta"), exists("fig_dir"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

near_sites <- site_meta$site_id[site_meta$region == "Nearctic"]
west_sites <- site_meta$site_id[site_meta$region == "Western Palearctic"]

near_bins <- colnames(comm_pa)[colSums(comm_pa[rownames(comm_pa) %in% near_sites, , drop = FALSE] > 0) > 0]
west_bins <- colnames(comm_pa)[colSums(comm_pa[rownames(comm_pa) %in% west_sites, , drop = FALSE] > 0) > 0]

shared <- intersect(near_bins, west_bins)
na_only <- setdiff(near_bins, west_bins)
wp_only <- setdiff(west_bins, near_bins)
total_union <- length(union(near_bins, west_bins)) 

df_bar <- tibble::tibble(
  group = c("Nearctic only", "Shared", "Western Palearctic only"),
  count = c(length(na_only), length(shared), length(wp_only))
) |>
  dplyr::mutate(
    pct         = 100 * count / total_union,
    label_count = format(count, big.mark = ","),
    label_pct   = sprintf("%.1f%%", pct)
  ) |>
  dplyr::arrange(dplyr::desc(count)) |>
  dplyr::mutate(group = factor(group, levels = group))

ymax <- max(df_bar$count)

ggplot2::ggplot(df_bar, ggplot2::aes(group, count)) +
  ggplot2::geom_col(width = 0.65, fill = "grey70", colour = "black") +
  ggplot2::geom_text(ggplot2::aes(y = count + 0.05 * ymax, label = label_count),
    size = 5, fontface = "bold"
  ) +
  ggplot2::annotate("text",
    x = 1:3, y = -0.07 * ymax,
    label = paste0("of all BINs: ", df_bar$label_pct), size = 3.8
  ) +
  ggplot2::coord_cartesian(ylim = c(-0.10 * ymax, 1.20 * ymax), clip = "off") +
  ggplot2::labs(
    title = "Unique vs shared Bombus BINs by region",
    x = NULL, y = "Number of BINs",
    caption = paste0("Total BINs = ", total_union)
  ) +
  ggplot2::theme_minimal(base_size = 14) +
  ggplot2::theme(
    plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold"),
    plot.caption     = ggplot2::element_text(hjust = 0.5, margin = ggplot2::margin(t = 10)),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x      = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    plot.margin      = ggplot2::margin(10, 20, 40, 20)
  )

ggplot2::ggsave(file.path(fig_dir, "C_unique_shared_BINs.png"),
  
  width = 6.5, height = 4.8, dpi = 300
)
# References:
# Base set ops: intersect(), union(), setdiff()
#   https://stat.ethz.ch/R-manual/R-patched/library/base/html/sets.html
# ggplot2 layers used here:
#   geom_col / geom_text / annotate / coord_cartesian / theme_minimal
#   https://ggplot2.tidyverse.org/reference/

# 6) For Results and discussion
cat("\n=== Balanced site summary ===\n")
site_meta %>%
  count(region, name = "n_sites") %>%
  print()
cat(sprintf("Sites per region used (balanced): %d\n", S))
cat(sprintf("Pairs per region (within): %d\n", S * (S - 1) / 2))
cat(sprintf("Mantel (Nearctic): r=%.3f, p=%.3f\n", mantel_na$statistic, mantel_na$signif))
cat(sprintf("Mantel (Western Palearctic): r=%.3f, p=%.3f\n", mantel_eu$statistic, mantel_eu$signif))
