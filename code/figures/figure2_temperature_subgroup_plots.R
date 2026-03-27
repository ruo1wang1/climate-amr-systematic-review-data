# Public release script for the climate-AMR manuscript

#     UNIFIED Single-Page Multi-Subgroup Forest Plot
#

# SETUP
rm(list = ls())
gc()
while(dev.cur() > 1) dev.off()

cat("\n================================================================\n")
cat("   Temperature & AMR: UNIFIED Single-Page Forest Plot\n")
cat("================================================================\n\n")

# Load packages
library(tidyverse)
library(grid)

# DIRECTORIES
project_root <- Sys.getenv("CLIMATE_AMR_GITHUB_ROOT", unset = ".")
base_dir <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
input_dir <- file.path(base_dir, "figure_source_data")
output_dir <- file.path(base_dir, "generated_figures", "figure2_temperature_subgroups")

if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# MORANDI COLOR PALETTE
MORANDI <- list(
  pathogen1 = "#D4A59A", pathogen2 = "#C9A997",
  pathogen3 = "#B8A390", pathogen4 = "#A89B8E",
  income1 = "#9ABACC", income2 = "#7A9BB0",
  region1 = "#A8B5A0", region2 = "#8FA589",
  bg = "#FAF9F7", text = "#2C2C2C",
  text_sub = "#5A5A5A", line = "#CCCCCC"
)

# READ AND PROCESS DATA
cat(">>> Reading data...\n")

data_raw <- read_csv(
  file.path(input_dir, "source_data_figure2b_temperature_annual_subgroup_summary.csv"),
  show_col_types = FALSE
)

data <- data_raw %>%
  rename(
    Analysis = 1, Subgroup_Variable = 2, Subgroup_Level = 3,
    k_effects = 4, n_studies = 5, Study_List = 6,
    Total_Isolates = 7, Isolates_Display = 8,
    Pooled_RR = 9, CI_Lower = 10, CI_Upper = 11,
    P_Value = 12, I2 = 13
  ) %>%
  filter(Analysis == "Annual Temperature", !is.na(Pooled_RR))

# Format data
data_formatted <- data %>%
  mutate(
    Subgroup_Order = case_when(
      Subgroup_Variable == "Pathogen" ~ 1,
      Subgroup_Variable == "Income Level" ~ 2,
      Subgroup_Variable == "Region" ~ 3,
      TRUE ~ 99
    ),

    Subgroup_Label = case_when(
      str_detect(Subgroup_Level, "baumannii") ~ "  A. baumannii",
      str_detect(Subgroup_Level, "coli") ~ "  E. coli",
      str_detect(Subgroup_Level, "pneumoniae") ~ "  K. pneumoniae",
      str_detect(Subgroup_Level, "aeruginosa") ~ "  P. aeruginosa",
      str_detect(Subgroup_Level, "High-incon") ~ "  High-income (HICs)",
      str_detect(Subgroup_Level, "Upper-midd") ~ "  Upper-middle-income (UMICs)",
      str_detect(Subgroup_Level, "East Asia") ~ "  East Asia & Pacific",
      str_detect(Subgroup_Level, "Europe") ~ "  Europe & Central Asia",
      TRUE ~ paste0("  ", Subgroup_Level)
    ),

    Isolates_Formatted = format(Total_Isolates, big.mark = ",", scientific = FALSE),
    RR_Text = sprintf("%.2f", Pooled_RR),
    CI_Text = sprintf("[%.2f, %.2f]", CI_Lower, CI_Upper),
    I2_Text = sprintf("%.1f%%", I2),

    Color = case_when(
      Subgroup_Variable == "Pathogen" & str_detect(Subgroup_Label, "baumannii") ~ MORANDI$pathogen1,
      Subgroup_Variable == "Pathogen" & str_detect(Subgroup_Label, "coli") ~ MORANDI$pathogen2,
      Subgroup_Variable == "Pathogen" & str_detect(Subgroup_Label, "pneumoniae") ~ MORANDI$pathogen3,
      Subgroup_Variable == "Pathogen" & str_detect(Subgroup_Label, "aeruginosa") ~ MORANDI$pathogen4,
      Subgroup_Variable == "Income Level" & str_detect(Subgroup_Label, "High") ~ MORANDI$income1,
      Subgroup_Variable == "Income Level" & str_detect(Subgroup_Label, "Upper") ~ MORANDI$income2,
      Subgroup_Variable == "Region" & str_detect(Subgroup_Label, "East") ~ MORANDI$region1,
      Subgroup_Variable == "Region" & str_detect(Subgroup_Label, "Europe") ~ MORANDI$region2,
      TRUE ~ "#999999"
    ),

    Het_P_Text = case_when(
      Subgroup_Variable == "Pathogen" ~ "P_heterogeneity = 0.0019",
      Subgroup_Variable == "Income Level" ~ "P_heterogeneity < 0.0001",
      Subgroup_Variable == "Region" ~ "P_heterogeneity = 0.7144",
      TRUE ~ ""
    ),

    Row_Type = "data"
  ) %>%
  arrange(Subgroup_Order, Subgroup_Label)

# Add section headers and spacing
plot_data <- tibble()

for(subgrp in c("Pathogen", "Income Level", "Region")) {

  # Section header
  header_row <- tibble(
    Subgroup_Variable = subgrp,
    Subgroup_Label = subgrp,
    k_effects = NA, n_studies = NA,
    Isolates_Formatted = "", RR_Text = "", CI_Text = "", I2_Text = "",
    Pooled_RR = NA, CI_Lower = NA, CI_Upper = NA,
    Color = MORANDI$text,
    Het_P_Text = unique(data_formatted$Het_P_Text[data_formatted$Subgroup_Variable == subgrp]),
    Row_Type = "header",
    Subgroup_Order = unique(data_formatted$Subgroup_Order[data_formatted$Subgroup_Variable == subgrp])
  )

  # Data rows
  data_rows <- data_formatted %>% filter(Subgroup_Variable == subgrp)

  # Heterogeneity p-value row
  het_row <- tibble(
    Subgroup_Variable = subgrp,
    Subgroup_Label = paste0("    ", unique(data_rows$Het_P_Text)),
    k_effects = NA, n_studies = NA,
    Isolates_Formatted = "", RR_Text = "", CI_Text = "", I2_Text = "",
    Pooled_RR = NA, CI_Lower = NA, CI_Upper = NA,
    Color = MORANDI$text_sub,
    Het_P_Text = "",
    Row_Type = "het_p",
    Subgroup_Order = unique(data_rows$Subgroup_Order)
  )

  # Spacer row
  spacer_row <- tibble(
    Subgroup_Variable = subgrp,
    Subgroup_Label = "",
    k_effects = NA, n_studies = NA,
    Isolates_Formatted = "", RR_Text = "", CI_Text = "", I2_Text = "",
    Pooled_RR = NA, CI_Lower = NA, CI_Upper = NA,
    Color = MORANDI$bg,
    Het_P_Text = "",
    Row_Type = "spacer",
    Subgroup_Order = unique(data_rows$Subgroup_Order)
  )

  plot_data <- bind_rows(plot_data, header_row, data_rows, het_row, spacer_row)
}

n_total_rows <- nrow(plot_data)

cat("    Total rows:", n_total_rows, "\n\n")

# DRAW UNIFIED FOREST PLOT

draw_unified_forest <- function() {

  grid.newpage()

  # Main viewport with margins
  pushViewport(viewport(
    x = 0.5, y = 0.5,
    width = 0.96, height = 0.96,
    just = c("center", "center")
  ))

  # Title
  grid.text(
    "Temperature and AMR: Multi-Subgroup Meta-Analysis",
    x = 0.5, y = 0.98,
    just = c("center", "top"),
    gp = gpar(fontsize = 14, fontface = "bold", col = MORANDI$text)
  )

  # Column headers
  y_header <- 0.94

  grid.text("Subgroup", x = 0.08, y = y_header, just = "left",
            gp = gpar(fontsize = 10, fontface = "bold", col = MORANDI$text))

  grid.text("No. of\neffect sizes", x = 0.30, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = MORANDI$text, lineheight = 0.9))

  grid.text("No. of\nstudies", x = 0.38, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = MORANDI$text, lineheight = 0.9))

  grid.text("No. of isolates", x = 0.48, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = MORANDI$text))

  grid.text("RR", x = 0.78, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = MORANDI$text))

  grid.text("95% CI", x = 0.86, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = MORANDI$text))

  grid.text("I²", x = 0.94, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = MORANDI$text))

  # Header underline
  grid.lines(
    x = c(0.04, 0.98),
    y = c(y_header - 0.015, y_header - 0.015),
    gp = gpar(col = MORANDI$line, lwd = 1.5)
  )

  # Data rows
  y_start <- y_header - 0.03
  row_height <- 0.88 / n_total_rows

  for(i in 1:n_total_rows) {

    row <- plot_data[i, ]
    y_pos <- y_start - (i * row_height)

    # Skip spacer rows (just empty space)
    if(row$Row_Type == "spacer") next

    # Section headers (bold)
    if(row$Row_Type == "header") {

      grid.text(
        row$Subgroup_Label,
        x = 0.08, y = y_pos,
        just = "left",
        gp = gpar(fontsize = 10, fontface = "bold", col = MORANDI$text)
      )

      # Horizontal line under section header
      grid.lines(
        x = c(0.08, 0.98),
        y = c(y_pos - 0.005, y_pos - 0.005),
        gp = gpar(col = MORANDI$line, lwd = 1, lty = 3)
      )

      next
    }

    # Heterogeneity p-value rows (italic)
    if(row$Row_Type == "het_p") {

      grid.text(
        row$Subgroup_Label,
        x = 0.08, y = y_pos,
        just = "left",
        gp = gpar(fontsize = 9, fontface = "italic", col = MORANDI$text_sub)
      )

      next
    }

    # Data rows
    if(row$Row_Type == "data") {

      # Subgroup label
      grid.text(
        row$Subgroup_Label,
        x = 0.08, y = y_pos,
        just = "left",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # Number of effect sizes
      grid.text(
        as.character(row$k_effects),
        x = 0.30, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # Number of studies
      grid.text(
        as.character(row$n_studies),
        x = 0.38, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # Number of isolates
      grid.text(
        row$Isolates_Formatted,
        x = 0.48, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # Forest plot
      plot_min <- 0.95
      plot_max <- 1.12
      plot_x_start <- 0.58
      plot_x_end <- 0.72

      rr_val <- row$Pooled_RR
      ci_lower <- row$CI_Lower
      ci_upper <- row$CI_Upper

      rr_x <- plot_x_start + ((rr_val - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)
      ci_lower_x <- plot_x_start + ((ci_lower - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)
      ci_upper_x <- plot_x_start + ((ci_upper - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

      ci_lower_x <- max(plot_x_start, min(plot_x_end, ci_lower_x))
      ci_upper_x <- max(plot_x_start, min(plot_x_end, ci_upper_x))
      rr_x <- max(plot_x_start, min(plot_x_end, rr_x))

      # CI line
      grid.lines(
        x = c(ci_lower_x, ci_upper_x),
        y = c(y_pos, y_pos),
        gp = gpar(col = row$Color, lwd = 2.5)
      )

      # Point estimate
      grid.rect(
        x = rr_x, y = y_pos,
        width = unit(0.012, "npc"),
        height = unit(0.012, "npc"),
        gp = gpar(fill = row$Color, col = row$Color)
      )

      # RR
      grid.text(
        row$RR_Text,
        x = 0.78, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # CI
      grid.text(
        row$CI_Text,
        x = 0.86, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )

      # I²
      grid.text(
        row$I2_Text,
        x = 0.94, y = y_pos,
        just = "center",
        gp = gpar(fontsize = 9, col = MORANDI$text)
      )
    }
  }

  # Forest plot axis
  y_axis <- 0.04

  # Null line (RR = 1.0)
  null_x <- plot_x_start + ((1.0 - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

  grid.lines(
    x = c(null_x, null_x),
    y = c(y_header - 0.02, y_axis + 0.02),
    gp = gpar(col = MORANDI$line, lwd = 1.5, lty = 2)
  )

  # X-axis
  grid.lines(
    x = c(plot_x_start, plot_x_end),
    y = c(y_axis, y_axis),
    gp = gpar(col = MORANDI$text, lwd = 1)
  )

  # Tick marks
  tick_values <- c(0.95, 1.00, 1.05, 1.10)
  for(tick in tick_values) {
    tick_x <- plot_x_start + ((tick - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

    grid.lines(
      x = c(tick_x, tick_x),
      y = c(y_axis, y_axis - 0.008),
      gp = gpar(col = MORANDI$text, lwd = 1)
    )

    grid.text(
      sprintf("%.2f", tick),
      x = tick_x, y = y_axis - 0.015,
      just = "center",
      gp = gpar(fontsize = 8, col = MORANDI$text)
    )
  }

  # X-axis label
  grid.text(
    "Relative Risk per 1°C increase in temperature",
    x = 0.65, y = y_axis - 0.03,
    just = "center",
    gp = gpar(fontsize = 10, col = MORANDI$text)
  )

  popViewport()
}

# SAVE PLOT

cat(">>> Saving unified forest plot...\n")

# PDF
pdf(file.path(output_dir, "FigC_Temperature_Unified_Forest.pdf"),
    width = 12, height = 11)
draw_unified_forest()
dev.off()

# PNG
png(file.path(output_dir, "FigC_Temperature_Unified_Forest.png"),
    width = 12, height = 11, units = "in", res = 600)
draw_unified_forest()
dev.off()

# EXPORT DATA TABLE

export_table <- plot_data %>%
  filter(Row_Type == "data") %>%
  select(
    Subgroup_Variable,
    Subgroup = Subgroup_Label,
    `No. of effect sizes` = k_effects,
    `No. of studies` = n_studies,
    `No. of isolates` = Isolates_Formatted,
    RR = RR_Text,
    `95% CI` = CI_Text,
    `I²` = I2_Text
  )

write_csv(
  export_table,
  file.path(output_dir, "Unified_Forest_Data_Table.csv")
)

# FINAL OUTPUT

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║     UNIFIED SINGLE-PAGE FOREST PLOT COMPLETE                ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("📊 Output Files:\n")

cat("🎨 Design Features:\n")

cat("📁 Location:", output_dir, "\n\n")
cat("🎉 Ready for journal submission!\n\n")

# END

#     Temperature and AMR: Seasonal Sensitivity Analysis
#     SIMPLIFIED - Effect Estimates Only (with I² on the right)
#

# SETUP
rm(list = ls())
gc()
while(dev.cur() > 1) dev.off()

cat("\n================================================================\n")
cat("   Temperature & AMR: Simplified Effect Estimates Figure\n")
cat("================================================================\n\n")

# Load packages
library(tidyverse)
library(grid)

# DIRECTORIES
project_root <- Sys.getenv("CLIMATE_AMR_GITHUB_ROOT", unset = ".")
base_dir <- normalizePath(project_root, winslash = "/", mustWork = FALSE)
input_dir <- file.path(base_dir, "Temperature_AMR_v28.9_Enhanced/Tables")
output_dir <- file.path(base_dir, "Temperature_AMR_v28.9_Enhanced/Multisubgroup_Figures")

if(!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# COLOR PALETTE
COLORS <- list(
  annual = "#4E99AF",
  summer = "#B18FC0",
  winter = "#FFE850",

  bg_light = "#FAFAFA",
  text_main = "#2C2C2C",
  text_sub = "#5A5A5A",
  line = "#D0D0D0",
  grid = "#E8E8E8"
)

# READ DATA
cat(">>> Reading data...\n")

overall_data <- read_csv(
  file.path(input_dir, "Overall_Analysis_Summary_Enhanced.csv"),
  show_col_types = FALSE
)

cat("    Data loaded:", nrow(overall_data), "rows\n\n")

# PROCESS DATA
cat(">>> Processing data...\n")

temp_data <- overall_data %>%
  rename(
    Analysis = 1, Effect_Type = 2, k_effects = 3, n_studies = 4,
    Total_Isolates = 5, Isolates_Display = 6,
    Pooled_Effect = 7, CI_Lower = 8, CI_Upper = 9, P_Value = 10,
    PI_Lower = 11, PI_Upper = 12, I2 = 13, Tau2 = 14, Tau = 15,
    H = 16, Q = 17, Q_pval = 18
  ) %>%
  filter(Analysis %in% c("Annual Temperature", "Summer Temperature", "Winter Temperature")) %>%
  mutate(
    Order = case_when(
      Analysis == "Annual Temperature" ~ 1,
      Analysis == "Summer Temperature" ~ 2,
      Analysis == "Winter Temperature" ~ 3,
      TRUE ~ 99
    ),

    # Scientific labels
    Label = case_when(
      Analysis == "Annual Temperature" ~ "Annual mean temperature",
      Analysis == "Summer Temperature" ~ "Summer temperature",
      Analysis == "Winter Temperature" ~ "Winter temperature",
      TRUE ~ Analysis
    ),

    # Color scheme
    Color = case_when(
      Analysis == "Annual Temperature" ~ COLORS$annual,
      Analysis == "Summer Temperature" ~ COLORS$summer,
      Analysis == "Winter Temperature" ~ COLORS$winter,
      TRUE ~ "#999999"
    ),

    # Format numbers
    Isolates_Formatted = format(Total_Isolates, big.mark = ",", scientific = FALSE),
    RR_Text = sprintf("%.2f", Pooled_Effect),
    CI_Text = sprintf("[%.2f, %.2f]", CI_Lower, CI_Upper),
    I2_Text = sprintf("%.1f%%", I2),

    # P-values
    P_Text = case_when(
      P_Value < 0.001 ~ "< 0.001",
      P_Value < 0.01 ~ sprintf("%.4f", P_Value),
      TRUE ~ sprintf("%.3f", P_Value)
    )
  ) %>%
  arrange(Order)

cat("    Processed:", nrow(temp_data), "temperature analyses\n\n")

# DRAW SIMPLIFIED FIGURE - EFFECT ESTIMATES ONLY

draw_simplified_figure <- function() {

  grid.newpage()

  # Main viewport with generous margins
  pushViewport(viewport(
    x = 0.5, y = 0.5,
    width = 0.96, height = 0.94,
    just = c("center", "center")
  ))

  # Main title
  grid.text(
    "Temperature and Antimicrobial Resistance: Effect Estimates",
    x = 0.5, y = 0.96,
    just = c("center", "top"),
    gp = gpar(fontsize = 14, fontface = "bold", col = COLORS$text_main)
  )

  # EFFECT ESTIMATES TABLE WITH I² ON THE RIGHT

  # Column headers
  y_header <- 0.88

  grid.text("Temperature\nMetric", x = 0.12, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = COLORS$text_main, lineheight = 0.9))

  grid.text("No. of\neffect sizes", x = 0.26, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = COLORS$text_main, lineheight = 0.9))

  grid.text("No. of\nstudies", x = 0.35, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = COLORS$text_main, lineheight = 0.9))

  grid.text("No. of\nisolates", x = 0.45, y = y_header, just = "center",
            gp = gpar(fontsize = 9, fontface = "bold", col = COLORS$text_main, lineheight = 0.9))

  # Forest plot header (empty space)
  # The forest plot will be between x = 0.54 and x = 0.74

  grid.text("RR", x = 0.80, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = COLORS$text_main))

  grid.text("95% CI", x = 0.88, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = COLORS$text_main))

  grid.text("I²", x = 0.96, y = y_header, just = "center",
            gp = gpar(fontsize = 10, fontface = "bold", col = COLORS$text_main))

  # Header underline
  grid.lines(
    x = c(0.02, 0.99),
    y = c(y_header - 0.025, y_header - 0.025),
    gp = gpar(col = COLORS$line, lwd = 2)
  )

  # Data rows
  n_rows <- nrow(temp_data)
  row_height <- 0.20  # Generous spacing
  y_start <- y_header - 0.06

  for(i in 1:n_rows) {

    row <- temp_data[i, ]
    y_pos <- y_start - ((i - 1) * row_height)

    # Temperature label (colored, bold)
    grid.text(
      row$Label,
      x = 0.12, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = row$Color, fontface = "bold")
    )

    # k_effects
    grid.text(
      as.character(row$k_effects),
      x = 0.26, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main)
    )

    # n_studies
    grid.text(
      as.character(row$n_studies),
      x = 0.35, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main)
    )

    # Isolates (with thousand separators)
    grid.text(
      row$Isolates_Formatted,
      x = 0.45, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main)
    )

    # FOREST PLOT AREA (x: 0.54 to 0.74)

    plot_min <- 0.99
    plot_max <- 1.09
    plot_x_start <- 0.54
    plot_x_end <- 0.74

    rr_val <- row$Pooled_Effect
    ci_lower <- row$CI_Lower
    ci_upper <- row$CI_Upper

    # Map to coordinates
    rr_x <- plot_x_start + ((rr_val - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)
    ci_lower_x <- plot_x_start + ((ci_lower - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)
    ci_upper_x <- plot_x_start + ((ci_upper - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

    # Clip to plotting area
    ci_lower_x <- max(plot_x_start, min(plot_x_end, ci_lower_x))
    ci_upper_x <- max(plot_x_start, min(plot_x_end, ci_upper_x))
    rr_x <- max(plot_x_start, min(plot_x_end, rr_x))

    # CI line (thicker)
    grid.lines(
      x = c(ci_lower_x, ci_upper_x),
      y = c(y_pos, y_pos),
      gp = gpar(col = row$Color, lwd = 3.5)
    )

    # Point estimate (diamond shape, larger)
    diamond_size <- 0.012
    grid.polygon(
      x = c(rr_x, rr_x + diamond_size, rr_x, rr_x - diamond_size),
      y = c(y_pos + diamond_size, y_pos, y_pos - diamond_size, y_pos),
      gp = gpar(fill = row$Color, col = row$Color)
    )

    # TEXT COLUMNS

    # RR
    grid.text(
      row$RR_Text,
      x = 0.80, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main)
    )

    # 95% CI
    grid.text(
      row$CI_Text,
      x = 0.88, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main)
    )

    # I² (rightmost column)
    grid.text(
      row$I2_Text,
      x = 0.96, y = y_pos,
      just = "center",
      gp = gpar(fontsize = 10, col = COLORS$text_main, fontface = "bold")
    )
  }

  # FOREST PLOT AXIS

  y_axis <- y_start - (n_rows * row_height) - 0.04

  # Null line (RR = 1.0, vertical dashed line)
  null_x <- plot_x_start + ((1.0 - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

  grid.lines(
    x = c(null_x, null_x),
    y = c(y_header - 0.03, y_axis + 0.01),
    gp = gpar(col = COLORS$line, lwd = 2, lty = 2)
  )

  # X-axis line
  grid.lines(
    x = c(plot_x_start, plot_x_end),
    y = c(y_axis, y_axis),
    gp = gpar(col = COLORS$text_main, lwd = 1.5)
  )

  # Tick marks and labels
  tick_values <- c(1.00, 1.02, 1.04, 1.06, 1.08)
  for(tick in tick_values) {
    if(tick >= plot_min && tick <= plot_max) {
      tick_x <- plot_x_start + ((tick - plot_min) / (plot_max - plot_min)) * (plot_x_end - plot_x_start)

      # Tick mark
      grid.lines(
        x = c(tick_x, tick_x),
        y = c(y_axis, y_axis - 0.012),
        gp = gpar(col = COLORS$text_main, lwd = 1.5)
      )

      # Tick label
      grid.text(
        sprintf("%.2f", tick),
        x = tick_x, y = y_axis - 0.028,
        just = "center",
        gp = gpar(fontsize = 9, col = COLORS$text_main)
      )
    }
  }

  # X-axis label
  grid.text(
    "Relative Risk per 1°C increase in temperature",
    x = 0.64, y = y_axis - 0.06,
    just = "center",
    gp = gpar(fontsize = 10, col = COLORS$text_main, fontface = "bold")
  )

  # FOOTER NOTE

  grid.text(
    "Note: Diamond size represents point estimate; horizontal lines indicate 95% confidence intervals; I² quantifies heterogeneity.",
    x = 0.5, y = 0.03,
    just = "center",
    gp = gpar(fontsize = 8, col = COLORS$text_sub, fontface = "italic")
  )

  popViewport()
}

# SAVE FIGURE

cat(">>> Saving simplified figure...\n")

# PDF (vector format for publication)
pdf(file.path(output_dir, "Fig_Temperature_Effect_Estimates_Simplified.pdf"),
    width = 12, height = 6)
draw_simplified_figure()
dev.off()

# PNG (high resolution for presentations)
png(file.path(output_dir, "Fig_Temperature_Effect_Estimates_Simplified.png"),
    width = 12, height = 6, units = "in", res = 600)
draw_simplified_figure()
dev.off()

# CREATE DATA TABLE

cat(">>> Creating data table...\n")

data_table <- temp_data %>%
  select(
    `Temperature Metric` = Label,
    `No. of effect sizes` = k_effects,
    `No. of studies` = n_studies,
    `No. of isolates` = Isolates_Formatted,
    `Pooled RR` = RR_Text,
    `95% CI` = CI_Text,
    `I²` = I2_Text,
    `P-value` = P_Text
  )

write_csv(
  data_table,
  file.path(output_dir, "Table_Temperature_Effect_Estimates.csv")
)

# CREATE ALTERNATIVE LAYOUT (HORIZONTAL)

cat(">>> Creating alternative horizontal layout...\n")

pdf(file.path(output_dir, "Fig_Temperature_Effect_Estimates_Horizontal.pdf"),
    width = 14, height = 5)

draw_simplified_figure()

dev.off()

png(file.path(output_dir, "Fig_Temperature_Effect_Estimates_Horizontal.png"),
    width = 14, height = 5, units = "in", res = 600)

draw_simplified_figure()

dev.off()

# SUMMARY STATISTICS

cat("\n╔══════════════════════════════════════════════════════════════╗\n")
cat("║     SIMPLIFIED EFFECT ESTIMATES FIGURE COMPLETE             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n\n")

cat("📊 Output Files:\n")
cat("   1️⃣  Fig_Temperature_Effect_Estimates_Simplified.pdf (12×6\")\n")
cat("      → Main figure (recommended for publication)\n\n")
cat("   2️⃣  Fig_Temperature_Effect_Estimates_Simplified.png (600 DPI)\n")
cat("      → High-resolution image\n\n")
cat("   3️⃣  Fig_Temperature_Effect_Estimates_Horizontal.pdf (14×5\")\n")
cat("      → Alternative horizontal layout\n\n")
cat("   4️⃣  Table_Temperature_Effect_Estimates.csv\n")
cat("      → Data table for supplementary materials\n\n")

cat("🎨 Design Features:\n")

cat("📈 Key Results Summary:\n")
cat("   ┌─────────────────────────────────────────────────────┐\n")
for(i in 1:nrow(temp_data)) {
  row <- temp_data[i, ]
  cat(sprintf("   │ %-25s RR = %.2f [%.2f, %.2f]  │\n",
              row$Label, row$Pooled_Effect, row$CI_Lower, row$CI_Upper))
  cat(sprintf("   │   Studies: %2d | Isolates: %15s │\n",
              row$n_studies, row$Isolates_Formatted))
  cat(sprintf("   │   I² = %5s | P-value: %8s        │\n",
              row$I2_Text, row$P_Text))
  if(i < nrow(temp_data)) {
    cat("   ├─────────────────────────────────────────────────────┤\n")
  }
}
cat("   └─────────────────────────────────────────────────────┘\n\n")

cat("📐 Table Layout:\n")
cat("   Temperature Metric | k_effects | n_studies | Isolates | [Forest] | RR | 95%CI | I²\n")
cat("   ──────────────────────────────────────────────────────────────────────────────\n")
cat("   Position (x):   0.12      0.26       0.35       0.45    0.54-0.74  0.80  0.88  0.96\n\n")

cat("📁 Location:", output_dir, "\n\n")
cat("🎉 Clean, focused, publication-ready!\n\n")

# END
