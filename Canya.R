install.packages(c("ggplot2", "patchwork", "dplyr", "tidyr"))
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Window start positions for each variant
# A30P is centered on residue 30 → windows 11:30 to 30:49
# E46K → 27:46 to 46:65
# H50Q/G51D → 31:50–50:69 etc.

data <- list(
  A53V = c(-0.0055, 0.0008, 0.0078, 0.0093, 0.0110, 0.0135, 0.0151, 0.0146,
           0.0147, 0.0125, 0.0121, 0.0119, 0.0234, 0.0227, 0.0197, 0.0170,
           0.0147, 0.0133, 0.0130, 0.0095),
  A53T = c(0.0045, 0.0046, 0.0013, 0.0017, 0.0030, 0.0026, 0.0023, 0.0021,
           -0.0002, 0.0009, -0.0008, -0.0012, 0.0003, 0.0005, 0.0033, 0.0035,
            0.0044, 0.0022, 0.0023, 0.0025),
  H50Q = c(0.0026,-0.0038,-0.0029,-0.0032,-0.0019,-0.0006, 0.0003, 0.0003,
           -0.0019,-0.0026,-0.0026,-0.0027,-0.0017,-0.0026,-0.0026,-0.0033,
           -0.0023,-0.0035,-0.0040,-0.0021),
  A30P = c(-0.0047,-0.0011,-0.0020,-0.0018,-0.0016,-0.0029,-0.0037,-0.0048,
           -0.0046,-0.0049,-0.0045,-0.0046,-0.0050,-0.0055,-0.0055,-0.0026,
           -0.0030,-0.0038,-0.0057,-0.0053),
  E46K = c(-0.0020,-0.0046,-0.0062,-0.0062,-0.0052,-0.0061,-0.0047,-0.0059,
           -0.0041,-0.0049,-0.0002,-0.0011,-0.0024,-0.0070,-0.0089,-0.0099,
           -0.0090,-0.0119,-0.0034, 0.0017),
  G51D = c(-0.0093,-0.0047,-0.0036,-0.0043,-0.0057,-0.0054,-0.0056,-0.0071,
           -0.0076,-0.0082,-0.0087,-0.0092,-0.0090,-0.0083,-0.0111,-0.0088,
           -0.0114,-0.0103,-0.0059,-0.0050),
  A53E = c(-0.0049,-0.0022,-0.0096,-0.0126,-0.0114,-0.0145,-0.0163,-0.0161,
           -0.0181,-0.0160,-0.0172,-0.0169,-0.0207,-0.0205,-0.0201,-0.0172,
           -0.0159,-0.0149,-0.0150,-0.0098)
)

# Window start positions (20 windows per variant)
positions <- list(
  A53V = 34:53, A53T = 34:53,
  H50Q = 31:50, A30P = 11:30,
  E46K = 27:46, G51D = 32:51,
  A53E = 34:53
)

mutation_sites <- c(A53V=53, A53T=53, H50Q=50, A30P=30, E46K=46, G51D=51, A53E=53)

df <- mapply(function(variant, deltas, pos) {
  data.frame(variant = variant, position = pos, delta = deltas)
}, names(data), data, positions, SIMPLIFY = FALSE) |> bind_rows()

# Factor order: strongest positive → strongest negative
variant_order <- c("A53V","A53T","H50Q","A30P","E46K","G51D","A53E")
df$variant <- factor(df$variant, levels = variant_order)

colors <- c(
  A53V = "#3d7ec9",
  A53T = "#5ba85a",
  H50Q = "#9b6fd4",
  A30P = "#e07b39",
  E46K = "#d44f6e",
  G51D = "#3aabaa",
  A53E = "#b5892e"
)

# Vertical lines at each mutation site
vlines <- data.frame(
  variant = factor(names(mutation_sites), levels = variant_order),
  xint    = as.numeric(mutation_sites)
)

p1 <- ggplot(df, aes(x = position, y = delta, color = variant)) +
  geom_hline(yintercept = 0, color = "#c0c8d4", linewidth = 0.5) +
  geom_vline(data = vlines, aes(xintercept = xint),
             linetype = "dashed", color = "#e8a838", linewidth = 0.6, alpha = 0.8) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  scale_color_manual(values = colors) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  facet_wrap(~variant, nrow = 2, scales = "free_x") +
  labs(
    title    = "Δ CANYA Score per Window (Mutant − Wild-Type)",
    subtitle = "Amber dashed line = mutation site residue",
    x = "Window start position", y = "Δ CANYA score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#e8ecf0"),
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(color = "grey50", size = 9)
  )

summary_df <- df |>
  group_by(variant) |>
  summarise(mean_delta = mean(delta), .groups = "drop") |>
  mutate(direction = ifelse(mean_delta > 0, "Increases aggregation", "Decreases aggregation"))

p2 <- ggplot(summary_df, aes(x = variant, y = mean_delta, fill = variant)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "#666", linewidth = 0.4) +
  geom_text(aes(
    label = sprintf("%+.4f", mean_delta),
    vjust = ifelse(mean_delta >= 0, -0.5, 1.4)
  ), size = 3, fontface = "bold") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  labs(
    title = "Mean Δ CANYA Score by Variant",
    x = NULL, y = "Mean Δ score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title        = element_text(face = "bold", size = 12)
  )
