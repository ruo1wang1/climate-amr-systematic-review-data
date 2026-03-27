# Public release script for the climate-AMR manuscript

# Panel A — WB region × income (FY26) background + dots (unweighted studies)

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  library(grid)
})

# 0) Paths
project_root <- Sys.getenv("CLIMATE_AMR_GITHUB_ROOT", unset = ".")
project_root <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
data_file <- file.path(project_root, "figure_source_data", "source_data_figure1a_map_country_counts.csv")
shp_file  <- file.path(project_root, "external_data", "world_boundary_country.shp")
out_dir   <- file.path(project_root, "generated_figures", "figure1a_map")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_pdf   <- file.path(out_dir, "figure1a_map_panelA.pdf")

# 1) Palette (unclassified -> light grey via na.value)
fill_values <- c(
  "East Asia & Pacific (HIC)"              = "#4E99AF",
  "East Asia & Pacific (LMIC)"             = "#C0E3F0",
  "Europe & Central Asia (HIC)"            = "#B18FC0",
  "Europe & Central Asia (LMIC)"           = "#D8CDDE",
  "Latin America & Caribbean (HIC)"        = "#5066AE",
  "Latin America & Caribbean (LMIC)"       = "#D7DDEC",
  "Middle East & North Africa (HIC)"       = "#FFE850",
  "Middle East & North Africa (LMIC)"      = "#EDD8A9",
  "North America (HIC)"                    = "#C85256",
  "North America (LMIC)"                   = "#D5B0AB",
  "South Asia (HIC)"                       = "#559E76",
  "South Asia (LMIC)"                      = "#D5B0AB",
  "Sub-Saharan Africa (HIC)"               = "#559E76",
  "Sub-Saharan Africa (LMIC)"              = "#DEECCA"
)

fill_breaks <- c(
  "East Asia & Pacific (HIC)", "East Asia & Pacific (LMIC)",
  "Europe & Central Asia (HIC)", "Europe & Central Asia (LMIC)",
  "Latin America & Caribbean (HIC)", "Latin America & Caribbean (LMIC)",
  "Middle East & North Africa (HIC)", "Middle East & North Africa (LMIC)",
  "North America (HIC)", "North America (LMIC)",
  "South Asia (HIC)", "South Asia (LMIC)",
  "Sub-Saharan Africa (HIC)", "Sub-Saharan Africa (LMIC)"
)

NA_FILL <- "#E6E6E6"  # light grey for unclassified

# 2) Read plotting data
dat <- read_csv(data_file, show_col_types = FALSE) %>%
  mutate(
    iso2 = toupper(str_trim(iso2)),
    iso3 = toupper(str_trim(iso3)),
    world_bank_region = str_trim(as.character(world_bank_region)),
    world_bank_income_hic_lmic = str_trim(as.character(world_bank_income_hic_lmic)),
    study_count_unweighted = as.numeric(study_count_unweighted)
  )

to_na_if_unclassified <- function(x){
  x2 <- str_trim(as.character(x))
  x2[x2 == ""] <- NA_character_
  x2[str_detect(tolower(x2), "not classified by world bank")] <- NA_character_
  x2
}

dat <- dat %>%
  mutate(
    WB_region_clean = to_na_if_unclassified(world_bank_region),
    income_clean    = to_na_if_unclassified(world_bank_income_hic_lmic),
    income_clean    = ifelse(income_clean %in% c("HIC","LMIC"), income_clean, NA_character_),
    fill_key = ifelse(
      is.na(WB_region_clean) | is.na(income_clean),
      NA_character_,
      paste0(WB_region_clean, " (", income_clean, ")")
    )
  )

# 3) Read shapefile & build join keys
world <- st_read(shp_file, quiet = TRUE) %>%
  st_make_valid() %>%
  st_transform(4326)

world_keyed <- world %>%
  mutate(
    shp_iso2 = toupper(str_trim(as.character(ISO_1))),   # iso2 candidate
    shp_iso3 = toupper(str_trim(as.character(ISO))),     # iso3 candidate
    shp_name = as.character(COUNTRY)
  )

# 4) Two-step join (iso2 first, then iso3 fallback)
world_A <- world_keyed %>%
  left_join(dat %>% select(iso2, study_count_unweighted, fill_key),
            by = c("shp_iso2" = "iso2")) %>%
  mutate(join_used = ifelse(!is.na(fill_key) | !is.na(study_count_unweighted), "iso2", NA_character_))

world_B_fill <- world_A %>%
  st_drop_geometry() %>%
  filter(is.na(fill_key) & is.na(study_count_unweighted)) %>%
  select(shp_iso3) %>%
  distinct() %>%
  left_join(dat %>% select(iso3, study_count_unweighted, fill_key),
            by = c("shp_iso3" = "iso3"))

world_joined <- world_A %>%
  left_join(world_B_fill, by = "shp_iso3", suffix = c("", "_iso3")) %>%
  mutate(
    fill_key = coalesce(fill_key, fill_key_iso3),
    study_count_unweighted = coalesce(study_count_unweighted, study_count_unweighted_iso3),
    join_used = coalesce(join_used, ifelse(!is.na(fill_key_iso3) | !is.na(study_count_unweighted_iso3), "iso3", NA_character_))
  ) %>%
  select(-ends_with("_iso3"))

world_joined <- world_joined %>%
  mutate(study_count_unweighted = ifelse(is.na(study_count_unweighted), 0, study_count_unweighted))

# 5) Dots
dots <- world_joined %>%
  filter(study_count_unweighted > 0) %>%
  mutate(geom_pt = st_point_on_surface(geometry)) %>%
  st_as_sf() %>%
  st_set_geometry("geom_pt")

# 6) Europe bbox + locator box
europe_bbox <- st_bbox(c(xmin = -12, ymin = 34, xmax = 45, ymax = 72), crs = st_crs(4326))
europe_box_sf <- st_as_sfc(europe_bbox) %>% st_sf(geometry = .)

world_eu <- st_crop(world_joined, europe_bbox)
dots_eu  <- st_crop(dots, europe_bbox)

# 7) Scales
size_breaks <- c(1, 2, 5, 10, 20, 50)

fill_scale <- scale_fill_manual(
  name     = "WB region × income (FY26)",
  values   = fill_values,
  breaks   = fill_breaks,
  drop     = FALSE,
  na.value = NA_FILL
)

size_scale <- scale_size_continuous(
  name   = "Number of studies (unweighted)",
  breaks = size_breaks,
  range  = c(1.5, 9),
  limits = c(1, max(size_breaks, max(dots$study_count_unweighted, na.rm = TRUE)))
)

# 8) Main map
p_main <- ggplot() +
  geom_sf(data = world_joined, aes(fill = fill_key), color = "black", linewidth = 0.20) +
  geom_sf(data = europe_box_sf, fill = NA, color = "black", linewidth = 0.35) +
  geom_sf(
    data = dots,
    aes(size = study_count_unweighted),
    shape = 21, fill = "grey20", color = "black",
    alpha = 0.70, stroke = 0.25
  ) +
  fill_scale +
  size_scale +
  coord_sf(xlim = c(-180, 180), ylim = c(-58, 85), expand = FALSE) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.key.height = unit(5, "mm"),
    legend.key.width  = unit(10, "mm"),
    plot.margin = margin(6, 6, 6, 6, unit = "mm"),
    plot.caption = element_text(size = 10, hjust = 0)
  ) +
  guides(
    fill = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(color = "black")),
    size = guide_legend(nrow = 1, byrow = TRUE)
  ) +
  labs(
    caption = paste0(
      "Panel A. Global distribution of included studies by study setting / sampling location. ",
      "Background colours represent World Bank analytical regions (FY26) stratified by income group (HIC vs LMIC). ",
      "Dots denote countries where included studies were conducted; dot size reflects the unweighted number of studies. ",
      "Countries/economies without World Bank region or income classification are shown in light grey."
    )
  )

# 9) Inset (legend OFF)
p_inset <- ggplot() +
  geom_sf(data = world_eu, aes(fill = fill_key), color = "black", linewidth = 0.20) +
  geom_sf(
    data = dots_eu,
    aes(size = study_count_unweighted),
    shape = 21, fill = "grey20", color = "black",
    alpha = 0.55, stroke = 0.30
  ) +
  fill_scale +
  size_scale +
  coord_sf(
    xlim = c(europe_bbox["xmin"], europe_bbox["xmax"]),
    ylim = c(europe_bbox["ymin"], europe_bbox["ymax"]),
    expand = FALSE
  ) +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  )

# 10) Compose (ONLY change: inset y lowered)
p_final <- ggdraw() +
  draw_plot(p_main,  x = 0.14, y = 0.0, width = 0.86, height = 1.0) +
  draw_plot(p_inset, x = 0.02, y = 0.42, width = 0.26, height = 0.30)  # <- moved down

# 11) Save
ggsave(
  filename = out_pdf,
  plot = p_final,
  device = cairo_pdf,
  width = 12,
  height = 7.2,
  units = "in",
  dpi = 300
)

message("Saved: ", out_pdf)
