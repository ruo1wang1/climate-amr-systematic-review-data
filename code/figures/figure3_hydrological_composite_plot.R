# Public release script for the climate-AMR manuscript

rm(list = ls())

packages <- c("tidyverse", "grid", "scales")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

project_root <- Sys.getenv("CLIMATE_AMR_GITHUB_ROOT", unset = ".")
base_path <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
data_path <- file.path(base_path, "figure_source_data")
output_path <- file.path(base_path, "generated_figures", "figure3_hydrological")

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

df_csv <- read.csv(
  file.path(data_path, "source_data_figure3a_hydrological_meta_summary.csv"),
  stringsAsFactors = FALSE
)

data_from_csv <- df_csv %>%
  filter(analysis_code %in% c("Precip_Beta", "Precip_Corr", "Humid_Beta", "Water_Corr")) %>%
  transmute(
    Analysis_Code = analysis_code,
    Hydrological_Factor = hydrological_factor,
    Effect_Type = effect_type,
    k_effects = k_effects,
    n_studies = n_studies,
    n_total_isolates = total_isolates,
    Pooled_Effect = pooled_effect,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper,
    P_Value = p_value,
    Original_Effect = pooled_effect,
    Original_CI_Lower = ci_lower,
    Original_CI_Upper = ci_upper
  )

runoff_data <- data.frame(
  Analysis_Code = c("Runoff_Surface_OR", "Runoff_Subsurface_OR"),
  Hydrological_Factor = c("Surface Runoff", "Sub-surface Runoff"),
  Effect_Type = c("OR", "OR"),
  k_effects = c(1, 1),
  n_studies = c(1, 1),
  n_total_isolates = c(32002684, 32002684),
  Pooled_Effect = c(1.22, 0.84),
  CI_Lower = c(1.05, 0.76),
  CI_Upper = c(1.4, 0.93),
  P_Value = c(0.04, 0.04), # Assuming <0.05 for both
  Original_Effect = c(1.22, 0.84),
  Original_CI_Lower = c(1.05, 0.76),
  Original_CI_Upper = c(1.4, 0.93)
)

combined_data <- bind_rows(data_from_csv, runoff_data)

processed_data <- combined_data %>%
  mutate(

    Transformed_Effect = ifelse(Effect_Type == "OR", log(Pooled_Effect), Pooled_Effect),
    Transformed_CI_Lower = ifelse(Effect_Type == "OR", log(CI_Lower), CI_Lower),
    Transformed_CI_Upper = ifelse(Effect_Type == "OR", log(CI_Upper), CI_Upper),

    is_significant = P_Value < 0.05,

    Shape_Type = case_when(
      Effect_Type == "OR" ~ "Log Odds Ratio",
      Effect_Type == "Correlation" ~ "Correlation (r)",
      Effect_Type == "Beta" ~ "Beta (β)",
      TRUE ~ "Other"
    ),

    Original_Value_CI = sprintf("%.3f (%.3f, %.3f)%s",
                                Original_Effect, Original_CI_Lower, Original_CI_Upper,
                                ifelse(is_significant, "*", "")),

    k_effects_display = as.character(k_effects),
    n_studies_display = as.character(n_studies),
    n_total_isolates_display = format(n_total_isolates, big.mark = ",", scientific = FALSE, trim = TRUE)
  ) %>%

  mutate(
    Factor_Label = case_when(
      Analysis_Code == "Precip_Corr" ~ "Precipitation (r)",
      Analysis_Code == "Runoff_Surface_OR" ~ "Runoff (Surface ln OR)",
      Analysis_Code == "Runoff_Subsurface_OR" ~ "Runoff (Sub-surface ln OR)",
      Analysis_Code == "Humid_Beta" ~ "Humidity (β)",
      Analysis_Code == "Precip_Beta" ~ "Precipitation (β)",
      Analysis_Code == "Water_Corr" ~ "Water Temperature (r)",
      TRUE ~ Hydrological_Factor
    )
  )

factor_order <- c(
  "Precipitation (r)",
  "Runoff (Surface ln OR)",
  "Runoff (Sub-surface ln OR)",
  "Humidity (β)",
  "Precipitation (β)",
  "Water Temperature (r)"
)

processed_data$Factor_Label <- factor(processed_data$Factor_Label, levels = rev(factor_order)) # rev() for plotting from top to bottom

processed_data <- processed_data %>% arrange(Factor_Label)

COLORS <- list(
  text_main = "gray20",
  text_sub = "gray50",
  line = "gray30",
  background_alt = "gray98",

  significant = "#E41A1C",
  non_significant = "gray70",

  "Precipitation (r)" = "#1f78b4",
  "Runoff (Surface ln OR)" = "#33a02c",
  "Runoff (Sub-surface ln OR)" = "#b2df8a",
  "Humidity (β)" = "#ff7f00",
  "Precipitation (β)" = "#6a3d9a",
  "Water Temperature (r)" = "#e31a1c"
)

FONT_SIZES <- list(
  title = 16,
  header = 10,
  data = 9,
  axis_label = 10,
  axis_tick = 9,
  footer = 8,
  legend = 9
)

SHAPES <- list(
  "Log Odds Ratio" = 23, # Diamond
  "Correlation (r)" = 21, # Circle
  "Beta (β)" = 22 # Square
)

draw_hydrological_figure <- function(data) {

  pushViewport(viewport(width = 0.95, height = 0.95))

  grid.text(
    "Comprehensive Comparison of Hydrological Factors on AMR",
    x = 0.02, y = 0.98, just = "left",
    gp = gpar(fontsize = FONT_SIZES$title, fontface = "bold", col = COLORS$text_main)
  )

  y_header <- 0.90

  grid.text("Variable", x = 0.10, y = y_header, just = "left",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  grid.text("k", x = 0.25, y = y_header, just = "center",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  grid.text("n", x = 0.32, y = y_header, just = "center",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  grid.text("Isolates", x = 0.40, y = y_header, just = "center",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  # Forest plot header (empty space, will be filled by axis label later)
  # The forest plot will be between x = 0.50 and x = 0.70

  grid.text("Value", x = 0.80, y = y_header, just = "center",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  grid.text("95% CI", x = 0.90, y = y_header, just = "center",
            gp = gpar(fontsize = FONT_SIZES$header, fontface = "bold", col = COLORS$text_main))

  grid.lines(
    x = c(0.02, 0.98),
    y = c(y_header - 0.02, y_header - 0.02),
    gp = gpar(col = COLORS$line, lwd = 1.5)
  )

  n_rows <- nrow(data)
  row_height <- 0.06
  y_start_data <- y_header - 0.05

  plot_x_min_val <- min(data$Transformed_CI_Lower, na.rm = TRUE) * 1.1
  plot_x_max_val <- max(data$Transformed_CI_Upper, na.rm = TRUE) * 1.1
  plot_x_range <- c(min(-0.7, plot_x_min_val), max(0.7, plot_x_max_val))

  plot_area_x_start <- 0.50
  plot_area_x_end <- 0.70

  for(i in 1:n_rows) {
    row <- data[i, ]
    y_pos <- y_start_data - ((i - 1) * row_height)

    if (i %% 2 == 0) {
      grid.rect(x = 0.5, y = y_pos, width = 0.96, height = row_height,
                gp = gpar(fill = COLORS$background_alt, col = NA))
    }

    # Variable Label (colored, bold italic)
    grid.text(
      as.character(row$Factor_Label),
      x = 0.10, y = y_pos,
      just = "left",
      gp = gpar(fontsize = FONT_SIZES$data, col = COLORS[[as.character(row$Factor_Label)]], fontface = "bold.italic")
    )

    # k_effects
    grid.text(
      row$k_effects_display,
      x = 0.25, y = y_pos,
      just = "center",
      gp = gpar(fontsize = FONT_SIZES$data, col = COLORS$text_main)
    )

    # n_studies
    grid.text(
      row$n_studies_display,
      x = 0.32, y = y_pos,
      just = "center",
      gp = gpar(fontsize = FONT_SIZES$data, col = COLORS$text_main)
    )

    # Isolates (with thousand separators)
    grid.text(
      row$n_total_isolates_display,
      x = 0.40, y = y_pos,
      just = "center",
      gp = gpar(fontsize = FONT_SIZES$data, col = COLORS$text_main)
    )

    # FOREST PLOT AREA (x: 0.50 to 0.70)

    effect_val <- row$Transformed_Effect
    ci_lower <- row$Transformed_CI_Lower
    ci_upper <- row$Transformed_CI_Upper

    # Map to coordinates
    effect_x <- plot_area_x_start + ((effect_val - plot_x_range[1]) / (plot_x_range[2] - plot_x_range[1])) * (plot_area_x_end - plot_area_x_start)
    ci_lower_x <- plot_area_x_start + ((ci_lower - plot_x_range[1]) / (plot_x_range[2] - plot_x_range[1])) * (plot_area_x_end - plot_area_x_start)
    ci_upper_x <- plot_area_x_start + ((ci_upper - plot_x_range[1]) / (plot_x_range[2] - plot_x_range[1])) * (plot_area_x_end - plot_area_x_start)

    ci_lower_x <- max(plot_area_x_start, min(plot_area_x_end, ci_lower_x))
    ci_upper_x <- max(plot_area_x_start, min(plot_area_x_end, ci_upper_x))
    effect_x <- max(plot_area_x_start, min(plot_area_x_end, effect_x))

    # CI line (thicker)
    grid.lines(
      x = c(ci_lower_x, ci_upper_x),
      y = c(y_pos, y_pos),
      gp = gpar(col = COLORS$text_main, lwd = 2)
    )

    # Point estimate (shape based on Effect_Type, fill based on significance)
    point_size <- 0.008

    fill_color <- ifelse(row$is_significant, COLORS$significant, COLORS$non_significant)

    if (row$Shape_Type == "Log Odds Ratio") { # Diamond
      grid.polygon(
        x = c(effect_x, effect_x + point_size, effect_x, effect_x - point_size),
        y = c(y_pos + point_size, y_pos, y_pos - point_size, y_pos),
        gp = gpar(fill = fill_color, col = COLORS$text_main, lwd = 0.5)
      )
    } else if (row$Shape_Type == "Correlation (r)") { # Circle
      grid.circle(
        x = effect_x, y = y_pos, r = point_size,
        gp = gpar(fill = fill_color, col = COLORS$text_main, lwd = 0.5)
      )
    } else if (row$Shape_Type == "Beta (β)") { # Square
      grid.rect(
        x = effect_x, y = y_pos, width = 2 * point_size, height = 2 * point_size,
        gp = gpar(fill = fill_color, col = COLORS$text_main, lwd = 0.5)
      )
    }

    # Value (95% CI)
    grid.text(
      row$Original_Value_CI,
      x = 0.85, y = y_pos,
      just = "center",
      gp = gpar(fontsize = FONT_SIZES$data, col = COLORS$text_main)
    )
  }

  # FOREST PLOT AXIS

  y_axis_pos <- y_start_data - (n_rows * row_height) - 0.03

  # Null line (0 = No Association, vertical dashed line)
  null_x_coord <- plot_area_x_start + ((0 - plot_x_range[1]) / (plot_x_range[2] - plot_x_range[1])) * (plot_area_x_end - plot_area_x_start)

  grid.lines(
    x = c(null_x_coord, null_x_coord),
    y = c(y_header - 0.02, y_axis_pos + 0.01),
    gp = gpar(col = COLORS$line, lwd = 1.5, lty = 2)
  )

  # X-axis line
  grid.lines(
    x = c(plot_area_x_start, plot_area_x_end),
    y = c(y_axis_pos, y_axis_pos),
    gp = gpar(col = COLORS$text_main, lwd = 1)
  )

  # Tick marks and labels
  tick_values <- pretty(plot_x_range, n = 5)
  tick_values <- tick_values[tick_values >= plot_x_range[1] & tick_values <= plot_x_range[2]]

  for(tick in tick_values) {
    tick_x <- plot_area_x_start + ((tick - plot_x_range[1]) / (plot_x_range[2] - plot_x_range[1])) * (plot_area_x_end - plot_area_x_start)

    # Tick mark
    grid.lines(
      x = c(tick_x, tick_x),
      y = c(y_axis_pos, y_axis_pos - 0.008),
      gp = gpar(col = COLORS$text_main, lwd = 1)
    )

    # Tick label
    grid.text(
      sprintf("%.2f", tick),
      x = tick_x, y = y_axis_pos - 0.02,
      just = "center",
      gp = gpar(fontsize = FONT_SIZES$axis_tick, col = COLORS$text_main)
    )
  }

  # X-axis label
  grid.text(
    "Standardized Effect (Direction)",
    x = (plot_area_x_start + plot_area_x_end) / 2, y = y_axis_pos - 0.04,
    just = "center",
    gp = gpar(fontsize = FONT_SIZES$axis_label, col = COLORS$text_main, fontface = "bold")
  )

  # Direction labels
  grid.text(
    "Reduced AMR Risk",
    x = plot_area_x_start, y = y_axis_pos - 0.06,
    just = "left",
    gp = gpar(fontsize = FONT_SIZES$axis_tick, col = COLORS$text_sub, fontface = "italic")
  )
  grid.text(
    "Increased AMR Risk",
    x = plot_area_x_end, y = y_axis_pos - 0.06,
    just = "right",
    gp = gpar(fontsize = FONT_SIZES$axis_tick, col = COLORS$text_sub, fontface = "italic")
  )

  legend_y_start <- y_axis_pos - 0.12
  legend_x_start <- 0.02

  grid.text("Legend:", x = legend_x_start, y = legend_y_start, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, fontface = "bold", col = COLORS$text_main))

  legend_y_pos <- legend_y_start - 0.025
  legend_x_offset_base <- 0.015

  # Log Odds Ratio (Diamond)
  grid.polygon(
    x = c(legend_x_start + legend_x_offset_base, legend_x_start + legend_x_offset_base + 0.004, legend_x_start + legend_x_offset_base, legend_x_start + legend_x_offset_base - 0.004),
    y = c(legend_y_pos + 0.004, legend_y_pos, legend_y_pos - 0.004, legend_y_pos),
    gp = gpar(fill = COLORS$significant, col = COLORS$text_main, lwd = 0.5)
  )
  grid.text("Log Odds Ratio", x = legend_x_start + legend_x_offset_base + 0.01, y = legend_y_pos, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, col = COLORS$text_main))

  # Correlation (r) (Circle)
  legend_x_offset_corr <- legend_x_offset_base + 0.08
  grid.circle(
    x = legend_x_start + legend_x_offset_corr, y = legend_y_pos, r = 0.004,
    gp = gpar(fill = COLORS$significant, col = COLORS$text_main, lwd = 0.5)
  )
  grid.text("Correlation (r)", x = legend_x_start + legend_x_offset_corr + 0.01, y = legend_y_pos, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, col = COLORS$text_main))

  # Beta (β) (Square)
  legend_x_offset_beta <- legend_x_offset_corr + 0.08
  grid.rect(
    x = legend_x_start + legend_x_offset_beta, y = legend_y_pos, width = 0.008, height = 0.008,
    gp = gpar(fill = COLORS$significant, col = COLORS$text_main, lwd = 0.5)
  )
  grid.text("Beta (β)", x = legend_x_start + legend_x_offset_beta + 0.01, y = legend_y_pos, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, col = COLORS$text_main))

  legend_y_pos_sig <- legend_y_pos - 0.02

  # P < 0.05 (Significant)
  grid.rect(
    x = legend_x_start + legend_x_offset_base, y = legend_y_pos_sig, width = 0.008, height = 0.008,
    gp = gpar(fill = COLORS$significant, col = COLORS$text_main, lwd = 0.5)
  )
  grid.text("P < 0.05", x = legend_x_start + legend_x_offset_base + 0.01, y = legend_y_pos_sig, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, col = COLORS$text_main))

  # P >= 0.05 (Non-significant)
  legend_x_offset_nonsig <- legend_x_offset_base + 0.08
  grid.rect(
    x = legend_x_start + legend_x_offset_nonsig, y = legend_y_pos_sig, width = 0.008, height = 0.008,
    gp = gpar(fill = COLORS$non_significant, col = COLORS$text_main, lwd = 0.5)
  )
  grid.text("P ≥ 0.05", x = legend_x_start + legend_x_offset_nonsig + 0.01, y = legend_y_pos_sig, just = "left",
            gp = gpar(fontsize = FONT_SIZES$legend, col = COLORS$text_main))

  # FOOTER NOTE

  grid.text(
    "Note: Point shapes indicate effect type; filled colors indicate statistical significance (P < 0.05).",
    x = 0.5, y = 0.03,
    just = "center",
    gp = gpar(fontsize = FONT_SIZES$footer, col = COLORS$text_sub, fontface = "italic")
  )

  popViewport()
}

# SAVE FIGURE

# PDF (vector format for publication)
pdf(file.path(output_path, "Fig_Hydrological_Effect_Estimates_Comprehensive.pdf"),
    width = 12, height = 8)
draw_hydrological_figure(processed_data)
dev.off()

# PNG (high resolution for presentations)
png(file.path(output_path, "Fig_Hydrological_Effect_Estimates_Comprehensive.png"),
    width = 12, height = 8, units = "in", res = 600)
draw_hydrological_figure(processed_data)
dev.off()

# CREATE DATA TABLE (for supplementary materials)

data_table_output <- processed_data %>%
  select(
    `Hydrological Factor` = Factor_Label,
    `Effect Type` = Shape_Type,
    `No. of effect sizes` = k_effects,
    `No. of studies` = n_studies,
    `No. of isolates` = n_total_isolates_display,
    `Pooled Effect` = Original_Effect,
    `CI Lower` = Original_CI_Lower,
    `CI Upper` = Original_CI_Upper,
    `P-value` = P_Value,
    `Significant (P<0.05)` = is_significant
  )

write_csv(
  data_table_output,
  file.path(output_path, "Table_Hydrological_Effect_Estimates.csv")
)

# SUMMARY STATISTICS

for(i in 1:nrow(processed_data)) {
  row <- processed_data[i, ]
  cat(sprintf("    │ %-30s | Type: %-18s | k: %-2s | n: %-2s | Isolates: %-15s | Value (95%% CI): %-25s | P: %.4f %s │\n",
              as.character(row$Factor_Label),
              as.character(row$Shape_Type),
              row$k_effects_display,
              row$n_studies_display,
              row$n_total_isolates_display,
              row$Original_Value_CI,
              row$P_Value,
              ifelse(row$is_significant, "(Significant)", "(Non-significant)")))
  if(i < nrow(processed_data)) {
    cat("    ├───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤\n")
  }
}
