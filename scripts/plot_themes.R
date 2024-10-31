library(ggplot2)

theme.base <- theme_minimal(base_size = 11) +
  theme(
    axis.text = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    plot.caption = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.5),
    panel.grid.minor = element_line(colour = adjustcolor("grey90", alpha.f = 0.4), linewidth = 0.5),
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.x = element_line(colour = "grey"),
    axis.line.y = element_line(colour = "grey"),
  )

theme.no_legend <- theme(legend.position = "none")

theme.no_grid <-  theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank()
)

theme.no_axis <- theme(
  axis.line.x = element_blank(),
  axis.line.y = element_blank()
)

theme.no_axis_title <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)

# Theme for timeseries
theme.ts <- theme.base +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme.no_legend
