#!/usr/bin/env Rscript
# Public release script for the climate-AMR manuscript

if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("scales")) install.packages("scales")

library(tidyverse)
library(ggplot2)
library(scales)

project_root <- Sys.getenv("CLIMATE_AMR_GITHUB_ROOT", unset = ".")
setwd(normalizePath(project_root, winslash = "/", mustWork = FALSE))

data <- read.csv(
  file.path(project_root, "figure_source_data", "source_data_figure3b_runoff_effect_estimates.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8"
)

runoff_data <- data %>%
  filter(!grepl("Average prevalence|Multiple", pathogen, ignore.case = TRUE)) %>%
  mutate(
    Runoff_Type = runoff_type,
    Pathogen_Clean = pathogen_display,
    AMR_Type = amr_type,
    Region_Type = region_type,
    Significant = is_significant,
    Pathogen_AMR = paste0(Pathogen_Clean, "\n(", AMR_Type, ")"),
    OR_Value = estimate,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  ) %>%
  filter(!is.na(Runoff_Type))

surface_data <- runoff_data %>%
  filter(Runoff_Type == "Surface") %>%
  select(Pathogen_Clean, AMR_Type, Region_Type,
         OR_Surface = OR_Value, CI_Lower_Surface = CI_Lower,
         CI_Upper_Surface = CI_Upper, Sig_Surface = Significant,
         Pathogen_AMR)

subsurface_data <- runoff_data %>%
  filter(Runoff_Type == "Sub-surface") %>%
  select(Pathogen_Clean, AMR_Type, Region_Type,
         OR_Subsurface = OR_Value, CI_Lower_Subsurface = CI_Lower,
         CI_Upper_Subsurface = CI_Upper, Sig_Subsurface = Significant,
         Pathogen_AMR)

paired_data <- surface_data %>%
  inner_join(subsurface_data,
             by = c("Pathogen_Clean", "AMR_Type", "Region_Type", "Pathogen_AMR")) %>%
  mutate(
    AMR_Order = ifelse(AMR_Type == "CR", 1, 2),
    Mean_OR = (OR_Surface + OR_Subsurface) / 2
  ) %>%
  arrange(AMR_Order, desc(Mean_OR)) %>%
  mutate(
    Pathogen_AMR = factor(Pathogen_AMR, levels = unique(Pathogen_AMR))
  )

color_surface <- "#E74C3C"
color_subsurface <- "#3498DB"
color_line <- "#BDC3C7"
color_ref <- "#7F8C8D"
color_sig_surface <- "#C0392B"
color_sig_subsurface <- "#2874A6"

p_final <- ggplot(paired_data) +

  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = color_ref,
             linewidth = 0.8,
             alpha = 0.7) +

  geom_segment(aes(x = OR_Subsurface, xend = OR_Surface,
                   y = Pathogen_AMR, yend = Pathogen_AMR),
               color = color_line,
               linewidth = 2.5,
               alpha = 0.5) +

  geom_errorbarh(aes(xmin = CI_Lower_Subsurface,
                     xmax = CI_Upper_Subsurface,
                     y = Pathogen_AMR),
                 height = 0,
                 color = color_subsurface,
                 linewidth = 2,
                 alpha = 0.75) +

  geom_errorbarh(aes(xmin = CI_Lower_Surface,
                     xmax = CI_Upper_Surface,
                     y = Pathogen_AMR),
                 height = 0,
                 color = color_surface,
                 linewidth = 2,
                 alpha = 0.75) +

  geom_point(aes(x = OR_Subsurface, y = Pathogen_AMR, shape = Region_Type),
             fill = color_subsurface,
             color = "white",
             size = 6.5,
             stroke = 1.5) +

  geom_point(aes(x = OR_Surface, y = Pathogen_AMR, shape = Region_Type),
             fill = color_surface,
             color = "white",
             size = 6.5,
             stroke = 1.5) +

  geom_text(data = paired_data %>% filter(Sig_Surface),
            aes(x = OR_Surface, y = Pathogen_AMR),
            label = "★",
            color = color_sig_surface,
            size = 6,
            vjust = -1.3,
            fontface = "bold") +

  geom_text(data = paired_data %>% filter(Sig_Subsurface),
            aes(x = OR_Subsurface, y = Pathogen_AMR),
            label = "★",
            color = color_sig_subsurface,
            size = 6,
            vjust = -1.3,
            fontface = "bold") +

  scale_x_continuous(
    name = "Relative Risk (RR) with 95% Confidence Interval",
    breaks = seq(0.5, 2.5, 0.5),
    limits = c(0.45, 2.55),
    expand = c(0.01, 0.01)
  ) +

  scale_y_discrete(name = NULL) +

  scale_shape_manual(
    name = "Data Source",
    values = c("Global" = 21, "Europe" = 24),
    labels = c("Global" = "Global", "Europe" = "Europe")
  ) +

  theme_minimal(base_size = 14) +
  theme(

    axis.title.x = element_text(
      face = "bold",
      size = 14,
      color = "#2C3E50",
      margin = margin(t = 15, b = 5)
    ),

    axis.text.y = element_text(
      size = 12.5,
      color = "#34495E",
      lineheight = 1.2,
      face = "plain"
    ),

    axis.text.x = element_text(
      size = 12,
      color = "#34495E"
    ),

    axis.line.x = element_line(
      color = "#2C3E50",
      linewidth = 0.8
    ),
    axis.line.y = element_line(
      color = "#95A5A6",
      linewidth = 0.5
    ),

    axis.ticks.x = element_line(
      color = "#2C3E50",
      linewidth = 0.6
    ),
    axis.ticks.y = element_line(
      color = "#95A5A6",
      linewidth = 0.4
    ),
    axis.ticks.length = unit(0.2, "cm"),

    panel.grid.major.x = element_line(
      color = "#ECF0F1",
      linewidth = 0.5
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(
      color = "#F8F9FA",
      linewidth = 0.4,
      linetype = "solid"
    ),

    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(
      face = "bold",
      size = 12,
      color = "#2C3E50"
    ),
    legend.text = element_text(
      size = 11,
      color = "#34495E"
    ),
    legend.key.size = unit(1.1, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = 20, b = 5),
    legend.key = element_rect(fill = "white", color = NA),

    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#BDC3C7", fill = NA, linewidth = 0.8),

    plot.margin = margin(15, 20, 15, 15)
  )

legend_data <- data.frame(
  x = rep(0.5, 4),
  y = rep(1, 4),
  label = c("Surface Runoff", "Sub-surface Runoff", "p < 0.05 (Surface)", "p < 0.05 (Sub-surface)"),
  color = c(color_surface, color_subsurface, color_sig_surface, color_sig_subsurface)
)

p_final <- p_final +

  annotate("point",
           x = 0.52, y = as.numeric(paired_data$Pathogen_AMR[1]) + 0.3,
           fill = color_surface, color = "white",
           shape = 21, size = 5.5, stroke = 1.2) +
  annotate("text",
           x = 0.57, y = as.numeric(paired_data$Pathogen_AMR[1]) + 0.3,
           label = "Surface Runoff",
           hjust = 0, size = 4, color = "#2C3E50", fontface = "bold") +

  annotate("point",
           x = 0.52, y = as.numeric(paired_data$Pathogen_AMR[1]) - 0.4,
           fill = color_subsurface, color = "white",
           shape = 21, size = 5.5, stroke = 1.2) +
  annotate("text",
           x = 0.57, y = as.numeric(paired_data$Pathogen_AMR[1]) - 0.4,
           label = "Sub-surface Runoff",
           hjust = 0, size = 4, color = "#2C3E50", fontface = "bold") +

  annotate("text",
           x = 0.52, y = as.numeric(paired_data$Pathogen_AMR[1]) - 1.1,
           label = "★", color = color_sig_surface,
           size = 5, hjust = 0, fontface = "bold") +
  annotate("text",
           x = 0.56, y = as.numeric(paired_data$Pathogen_AMR[1]) - 1.1,
           label = "p < 0.05 (Surface)",
           hjust = 0, size = 3.5, color = "#2C3E50") +

  annotate("text",
           x = 0.52, y = as.numeric(paired_data$Pathogen_AMR[1]) - 1.8,
           label = "★", color = color_sig_subsurface,
           size = 5, hjust = 0, fontface = "bold") +
  annotate("text",
           x = 0.56, y = as.numeric(paired_data$Pathogen_AMR[1]) - 1.8,
           label = "p < 0.05 (Sub-surface)",
           hjust = 0, size = 3.5, color = "#2C3E50")

ggsave(
  "Fig_Runoff_Final.pdf",
  plot = p_final,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  "Fig_Runoff_Final.png",
  plot = p_final,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "Fig_Runoff_Final_HighRes.pdf",
  plot = p_final,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  device = cairo_pdf
)

ggsave(
  "Fig_Runoff_Final_HighRes.png",
  plot = p_final,
  width = 12,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

p_standard <- ggplot(paired_data) +

  geom_vline(xintercept = 1, linetype = "dashed",
             color = color_ref, linewidth = 0.8, alpha = 0.7) +

  geom_segment(aes(x = OR_Subsurface, xend = OR_Surface,
                   y = Pathogen_AMR, yend = Pathogen_AMR),
               color = color_line, linewidth = 2.5, alpha = 0.5) +

  geom_errorbarh(aes(xmin = CI_Lower_Subsurface, xmax = CI_Upper_Subsurface,
                     y = Pathogen_AMR),
                 height = 0, color = color_subsurface,
                 linewidth = 2, alpha = 0.75) +

  geom_errorbarh(aes(xmin = CI_Lower_Surface, xmax = CI_Upper_Surface,
                     y = Pathogen_AMR),
                 height = 0, color = color_surface,
                 linewidth = 2, alpha = 0.75) +

  geom_point(aes(x = OR_Subsurface, y = Pathogen_AMR,
                 shape = Region_Type, fill = "Sub-surface"),
             color = "white", size = 6.5, stroke = 1.5) +

  geom_point(aes(x = OR_Surface, y = Pathogen_AMR,
                 shape = Region_Type, fill = "Surface"),
             color = "white", size = 6.5, stroke = 1.5) +

  geom_text(data = paired_data %>% filter(Sig_Surface),
            aes(x = OR_Surface, y = Pathogen_AMR),
            label = "★", color = color_sig_surface,
            size = 6, vjust = -1.3, fontface = "bold") +

  geom_text(data = paired_data %>% filter(Sig_Subsurface),
            aes(x = OR_Subsurface, y = Pathogen_AMR),
            label = "★", color = color_sig_subsurface,
            size = 6, vjust = -1.3, fontface = "bold") +

  scale_x_continuous(
    name = "Relative Risk (RR) with 95% Confidence Interval",
    breaks = seq(0.5, 2.5, 0.5),
    limits = c(0.45, 2.55),
    expand = c(0.01, 0.01)
  ) +

  scale_y_discrete(name = NULL) +

  scale_shape_manual(
    name = "Data Source",
    values = c("Global" = 21, "Europe" = 24)
  ) +

  scale_fill_manual(
    name = "Runoff Type",
    values = c("Surface" = color_surface, "Sub-surface" = color_subsurface),
    labels = c("Surface" = "Surface Runoff", "Sub-surface" = "Sub-surface Runoff")
  ) +

  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 5)),
    shape = guide_legend(order = 2, override.aes = list(size = 5))
  ) +

  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14,
                                color = "#2C3E50", margin = margin(t = 15, b = 5)),
    axis.text.y = element_text(size = 12.5, color = "#34495E", lineheight = 1.2),
    axis.text.x = element_text(size = 12, color = "#34495E"),
    axis.line.x = element_line(color = "#2C3E50", linewidth = 0.8),
    axis.line.y = element_line(color = "#95A5A6", linewidth = 0.5),
    axis.ticks.x = element_line(color = "#2C3E50", linewidth = 0.6),
    axis.ticks.y = element_line(color = "#95A5A6", linewidth = 0.4),
    axis.ticks.length = unit(0.2, "cm"),
    panel.grid.major.x = element_line(color = "#ECF0F1", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "#F8F9FA", linewidth = 0.4),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(face = "bold", size = 12, color = "#2C3E50"),
    legend.text = element_text(size = 11, color = "#34495E"),
    legend.key.size = unit(1.1, "cm"),
    legend.spacing.x = unit(0.6, "cm"),
    legend.margin = margin(t = 20, b = 5),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "#BDC3C7", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 20, 15, 15)
  ) +

  labs(caption = "★ p < 0.05") +
  theme(plot.caption = element_text(hjust = 0.5, size = 11,
                                    color = "#2C3E50", margin = margin(t = 10)))

ggsave(
  "Fig_Runoff_Standard_Legend.pdf",
  plot = p_standard,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)