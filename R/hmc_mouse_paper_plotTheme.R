pacman::p_load(
  broom,
  cowplot,
  foreach,
  data.table,
  dplyr,
  ggforce,
  ggplot2,
  ggpubr
)
theme_set(theme_gray()) # Disable cowplot theme

color_grey_1 <- "grey60"
color_grey_2 <- "grey75"
color_grey_3 <- "grey90"
color_other <- "#c6d3ba"
color_border <- "black"
outlier_color <- "grey10"
outlier_shape <- 16
outlier_size <- 0.2
font_size_main <- 7
font_size_gete <- 2.4

theme_publication <- function(font = font_size_main) {
  theme(
    text = element_text(family = "Helvetica", size = font, color = "black"),

    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", colour = NA),

    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    axis.line = element_line(color = "black"),
    axis.text = element_text(size = font, color = "black"), 
    axis.ticks = element_line(size = 0.3, color = "black"),
    axis.title = element_text(size = font, color = "black"),

    strip.background = element_rect(color = "white", fill = "white"),
    strip.text.x = element_text(size = font, color = "black"),
    strip.text.y = element_text(size = font, color = "black", angle = 0),

    legend.background = element_rect(
      fill = alpha("white", 0),
      color = alpha("white", 0)
    ),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.2, "cm"),    
    legend.text = element_text(size = font)
  )
}