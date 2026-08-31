options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(scales)
  library(svglite)
  library(ragg)
})

all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", all_args[grep("^--file=", all_args)])
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
output_root <- Sys.getenv("ANALYSIS_OUTPUT_DIR", unset = file.path(script_dir, "outputs"))
analysis_dir <- file.path(output_root, "analysis")
figure_dir <- file.path(output_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

blue <- "#0F4D92"
blue2 <- "#3775BA"
rose <- "#B64342"
teal <- "#42949E"
neutral <- "#767676"
light <- "#D9E2F2"

theme_set(
  theme_classic(base_size = 7.5, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6.8),
      strip.text = element_text(size = 7, face = "bold"),
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 6.8, colour = neutral),
      plot.tag = element_text(size = 9, face = "bold"),
      panel.grid = element_blank()
    )
)

save_pub <- function(plot, stem, width_mm = 183, height_mm = 120, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite::svglite(file.path(figure_dir, paste0(stem, ".svg")), width = w, height = h)
  print(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(figure_dir, paste0(stem, ".pdf")), width = w, height = h, family = "sans")
  print(plot)
  dev.off()
  ragg::agg_tiff(file.path(figure_dir, paste0(stem, ".tiff")), width = w, height = h, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()
  ragg::agg_png(file.path(figure_dir, paste0(stem, ".png")), width = w, height = h, units = "in", res = 300)
  print(plot)
  dev.off()
}

# Figure 1: persistence and accumulation of the overall deficit
curve <- read_csv(file.path(analysis_dir, "figure1_overall_survival_curve_source.csv"), show_col_types = FALSE) |>
  rename(month = Interval, observed = `Observed Survival (Cum)`, expected = `Expected Survival (Cum)`,
         obs_lo = `Observed Lower CI (Cum)`, obs_hi = `Observed Upper CI (Cum)`)

p1a <- ggplot(curve, aes(x = month)) +
  geom_ribbon(aes(ymin = obs_lo, ymax = obs_hi), fill = light, alpha = 0.8, na.rm = TRUE) +
  geom_line(aes(y = expected, colour = "Expected"), linewidth = 0.75) +
  geom_line(aes(y = observed, colour = "Observed"), linewidth = 0.75) +
  scale_colour_manual(values = c("Observed" = blue, "Expected" = neutral)) +
  scale_x_continuous(breaks = seq(0, 180, 36), labels = seq(0, 15, 3)) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0.60, 1.01)) +
  labs(x = "Years since diagnosis", y = "Cumulative survival", colour = NULL,
       subtitle = "15-year cohort: diagnoses 2004-2008 (n = 91,752)") +
  theme(legend.position = c(0.18, 0.18), legend.background = element_blank(),
        plot.margin = margin(4, 6, 4, 4))

overall <- bind_rows(
  read_csv(file.path(analysis_dir, "overall_5y_rlel.csv"), show_col_types = FALSE) |> filter(tau_months == 60),
  read_csv(file.path(analysis_dir, "overall_10y_rlel.csv"), show_col_types = FALSE) |> filter(tau_months == 120),
  read_csv(file.path(analysis_dir, "overall_15y_rlel.csv"), show_col_types = FALSE) |> filter(tau_months == 180)
) |>
  mutate(horizon = factor(paste0(tau_months / 12, " y"), levels = c("5 y", "10 y", "15 y")))

p1b <- ggplot(overall, aes(x = horizon, y = rlel_months)) +
  geom_col(width = 0.58, fill = blue2) +
  geom_errorbar(aes(ymin = rlel_ci_lower, ymax = rlel_ci_upper), width = 0.12, linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.2f", rlel_months)), vjust = -0.8, size = 2.5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Restriction horizon", y = "Restricted loss in life expectancy (months)",
       subtitle = "Maximum-follow-up cohort for each horizon")

fig1 <- p1a + p1b + plot_annotation(tag_levels = "a") + plot_layout(widths = c(1.65, 1))
save_pub(fig1, "Figure1_overall_survival_deficit", 183, 105)

# Figure 2: tumor-level heterogeneity
het <- read_csv(file.path(analysis_dir, "figure2_heterogeneity_source.csv"), show_col_types = FALSE) |>
  mutate(group = factor(group, levels = rev(group[order(rlel_months)])),
         domain = factor(domain, levels = c("Behavior", "Brain/CNS recode", "C75 site")))

p2 <- ggplot(het, aes(x = rlel_months, y = group, colour = domain)) +
  geom_vline(xintercept = 0, colour = "#BDBDBD", linewidth = 0.4, linetype = 2) +
  geom_errorbar(aes(xmin = rlel_ci_lower, xmax = rlel_ci_upper),
                orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.1) +
  geom_text(aes(label = paste0("n=", comma(n))), colour = neutral, hjust = -0.05, nudge_x = 0.5, size = 2.2) +
  scale_colour_manual(values = c("Behavior" = rose, "Brain/CNS recode" = blue, "C75 site" = teal)) +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.19))) +
  labs(x = "10-year restricted loss in life expectancy (months)", y = NULL, colour = NULL) +
  theme(legend.position = "top", axis.text.y = element_text(size = 6.8))
save_pub(p2, "Figure2_heterogeneity", 183, 165)

# Figure 3: age gradient
age <- read_csv(file.path(analysis_dir, "age_10y_rlel.csv"), show_col_types = FALSE) |>
  filter(tau_months == 120) |>
  mutate(age_group = factor(`Age recode with <1 year olds and 90+_label`, levels = `Age recode with <1 year olds and 90+_label`))

p3 <- ggplot(age, aes(x = age_group, y = rlel_months, group = 1)) +
  geom_ribbon(aes(ymin = rlel_ci_lower, ymax = rlel_ci_upper, group = 1), fill = light, alpha = 0.9) +
  geom_line(colour = blue, linewidth = 0.75) +
  geom_point(aes(size = n), colour = blue, alpha = 0.9) +
  scale_size_continuous(range = c(1.5, 4.2), breaks = c(500, 5000, 20000), labels = comma) +
  labs(x = "Age at diagnosis", y = "10-year restricted loss in life expectancy (months)", size = "Cases") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1), legend.position = "top")
save_pub(p3, "Figure3_age_gradient", 183, 105)

# Figure 4: diagnosis-year trend
std <- read_csv(file.path(analysis_dir, "age_standardized_dxyear_5y.csv"), show_col_types = FALSE) |>
  transmute(year, estimate = standardized_rlel_5y_months, lo = ci_lower, hi = ci_upper,
            series = "Age-standardized")
crude <- read_csv(file.path(analysis_dir, "dxyear_5y_rlel.csv"), show_col_types = FALSE) |>
  mutate(year = suppressWarnings(as.numeric(`Year of diagnosis_label`))) |>
  filter(year >= 2004, year <= 2018, tau_months == 60) |>
  transmute(year, estimate = rlel_months, lo = rlel_ci_lower, hi = rlel_ci_upper, series = "Crude")
trend <- bind_rows(std, crude)

p4 <- ggplot(trend, aes(x = year, y = estimate, colour = series, fill = series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.10, colour = NA) +
  geom_line(linewidth = 0.75) +
  geom_point(size = 1.7) +
  geom_smooth(data = subset(trend, series == "Age-standardized"), method = "lm", se = FALSE,
              linewidth = 0.55, linetype = 2, colour = blue) +
  annotate("text", x = 2004.2, y = max(trend$hi) + 0.03,
           label = "Age-standardized slope: -0.068 months/year (HC3 95% CI -0.088 to -0.047)",
           hjust = 0, size = 2.4, colour = blue) +
  scale_colour_manual(values = c("Age-standardized" = blue, "Crude" = neutral)) +
  scale_fill_manual(values = c("Age-standardized" = blue, "Crude" = neutral)) +
  scale_x_continuous(breaks = seq(2004, 2018, 2)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.15))) +
  labs(x = "Year of diagnosis", y = "5-year restricted loss in life expectancy (months)",
       colour = NULL, fill = NULL) +
  theme(legend.position = "top")
save_pub(p4, "Figure4_diagnosis_year_trend", 183, 105)

write_lines(capture.output(sessionInfo()), file.path(figure_dir, "R_sessionInfo.txt"))
message("Figures written to: ", figure_dir)
