# =============================================================================
# Polypharmacy KM panel figure using original cohort data.
#
# Packages: install.packages(c("survival", "survminer", "cowplot", "ggplot2",
# "showtext", "sysfonts", "magick", "haven", "dplyr"))
#
# Input requirements:
# - A real cohort file (not included in this repository)
# - Columns used below (same preprocessing conventions as main analysis script)
# - Drug category levels: "0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(cowplot)
  library(showtext)
  library(magick)
  library(sysfonts)
  library(haven)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# Real data input (set path externally; file is not part of this repository)
# -----------------------------------------------------------------------------
source_data_path <- Sys.getenv("POLYPHARMACY_SOURCE_DATA", unset = "")
if (source_data_path == "") {
  stop(
    paste(
      "Set env var POLYPHARMACY_SOURCE_DATA to the original cohort file path",
      "(e.g. a .sav export) before running this script."
    )
  )
}

cohort_df <- read_sav(source_data_path) %>% as.data.frame()
fixed_end_date <- as.Date("2021-12-31")

required_columns <- c(
  "DateISorLoppuOrDeath",
  "DateLoppuOrDeath",
  "CohortEntryDate",
  "ISaftercohortall",
  "Kuollut",
  "FirsteverICHAfterCohort",
  "RecurrentICHAfterCohort",
  "kuolpvmSPSSdate",
  "CombinedFirstEverRecurrentAfterCohort",
  "Unique_Drugs_ATC_Category"
)
missing_columns <- setdiff(required_columns, names(cohort_df))
if (length(missing_columns) > 0) {
  stop(
    sprintf(
      "Input dataset is missing required columns: %s",
      paste(missing_columns, collapse = ", ")
    )
  )
}

level_order <- c("0\u20132 drugs", "3\u20134 drugs", "5\u20139 drugs", "10+ drugs")
cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_Category,
  levels = level_order
)

# -----------------------------------------------------------------------------
# Real-data survival objects
# -----------------------------------------------------------------------------
cohort_df$duration_stroke <- as.numeric(
  difftime(cohort_df$DateISorLoppuOrDeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25
tmp_stroke <- cohort_df %>% dplyr::filter(!is.na(duration_stroke) & duration_stroke >= 0)
surv_stroke <- with(tmp_stroke, Surv(duration_stroke, ISaftercohortall == 1))
km_stroke <- survfit(surv_stroke ~ Unique_Drugs_ATC_Category, data = tmp_stroke)

cohort_df$duration_mortality <- as.numeric(
  difftime(cohort_df$DateLoppuOrDeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25
tmp_mortality <- cohort_df %>% dplyr::filter(!is.na(duration_mortality) & duration_mortality >= 0)
surv_mortality <- with(tmp_mortality, Surv(duration_mortality, Kuollut == 1))
km_mortality <- survfit(surv_mortality ~ Unique_Drugs_ATC_Category, data = tmp_mortality)

cohort_df$ICHfirstdayaftercohortfirstorrecurrent <- pmin(
  cohort_df$FirsteverICHAfterCohort,
  cohort_df$RecurrentICHAfterCohort,
  na.rm = TRUE
)
cohort_df$ENDdateICH_loppuICHordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$ICHfirstdayaftercohortfirstorrecurrent,
  fixed_end_date,
  na.rm = TRUE
)
cohort_df$ICHallaftercohort <- ifelse(
  !is.na(cohort_df$ICHfirstdayaftercohortfirstorrecurrent), 1, 0
)
cohort_df$duration_ich <- as.numeric(
  difftime(cohort_df$ENDdateICH_loppuICHordeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25
tmp_ich <- cohort_df %>% dplyr::filter(!is.na(duration_ich) & duration_ich >= 0)
surv_ich <- with(tmp_ich, Surv(duration_ich, ICHallaftercohort == 1))
km_ich <- survfit(surv_ich ~ Unique_Drugs_ATC_Category, data = tmp_ich)

cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$CombinedFirstEverRecurrentAfterCohort,
  fixed_end_date,
  na.rm = TRUE
)
cohort_df$Anybleedallaftercohort <- ifelse(
  !is.na(cohort_df$CombinedFirstEverRecurrentAfterCohort), 1, 0
)
cohort_df$duration_bleed <- as.numeric(
  difftime(cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25
tmp_bleed <- cohort_df %>% dplyr::filter(!is.na(duration_bleed) & duration_bleed >= 0)
surv_bleed <- with(tmp_bleed, Surv(duration_bleed, Anybleedallaftercohort == 1))
km_bleed <- survfit(surv_bleed ~ Unique_Drugs_ATC_Category, data = tmp_bleed)

# --- Font setup ---
font_add_google("Rosario", "rosario")
showtext_auto()

# =====================================================
# Common Plot Settings
# =====================================================
# Use one color per drug category
custom_colors <- c("#588692", "#E6A960", "#B96D40", "#8C3F5D")

plot_fun <- function(fit, data, title_label, ylim_value) {
  p <- ggsurvplot(
    fit = fit,
    data = data,
    fun = "event",
    conf.int = TRUE,
    risk.table = FALSE,
    censor = FALSE,
    palette = c("#588692", "#E6A960", "#B96D40", "#8C3F5D"),  # optional custom palette
    pval = FALSE,
    xlab = "Follow-up time (years)",
    ylab = "Cumulative incidence",
    break.x.by = 1,
    legend.title = "Number of drugs",  # ✅ Clean title
    legend.labs = c("0-2", "3-4", "5-9", "10+"),  # ✅ No variable name prefix
    legend = c(0.5, 1),
    ggtheme = theme_classic(base_size = 14, base_family = "rosario") +
      theme(
        legend.title = element_text(size = 15, face = "plain"),
        legend.text = element_text(size = 14, face = "plain"),
        axis.title = element_text(size = 17, face = "bold"),
        axis.text = element_text(size = 16),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(margin = margin(t = 18)),
        axis.title.y = element_text(margin = margin(r = 18)),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = rgb(1, 1, 1, alpha = 0.85), colour = NA)
      )
  )

  # Zoom to 5 years and adjust y-limits
  p$plot <- p$plot +
    coord_cartesian(xlim = c(0, 5), ylim = c(0, ylim_value), clip = "off") +
    ggtitle(title_label) +
    theme(
      plot.title = element_text(face = "bold", size = 17),
      legend.position = "inside",
      legend.position.inside = c(0.5, 1),
      legend.justification = c(0.5, 1)
    )

  return(p)
}

# =====================================================
# Create plots A–D (no risk tables)
# =====================================================
plotA <- plot_fun(km_stroke, tmp_stroke, "A. Ischemic stroke", 0.15)
plotB <- plot_fun(km_mortality, tmp_mortality, "B. Mortality", 0.60)
plotC <- plot_fun(km_ich, tmp_ich, "C. Intracranial hemorrhage", 0.10)
plotD <- plot_fun(km_bleed, tmp_bleed, "D. Major bleeding", 0.30)

# =====================================================
# Combine into panel (no risk tables)
# =====================================================
null_plot <- ggplot() + theme_void()

panel_plot <- plot_grid(
  plotA$plot, null_plot,  plotB$plot,
  null_plot, null_plot,   null_plot,
  plotC$plot, null_plot,  plotD$plot,
  ncol = 3,
  rel_widths = c(1, 0.07, 1),
  rel_heights = c(1, 0.07, 1),
  labels = NULL
)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------

ggsave("km_panel_polypharmacy_polished.pdf",
       plot = panel_plot, width = 16, height = 10)

pdf_image <- image_read_pdf("km_panel_polypharmacy_polished.pdf", density = 300)
image_write(pdf_image, path = "km_panel_polypharmacy_polished.png", format = "png", density = 300)
image_write(pdf_image, path = "km_panel_polypharmacy_polished.tiff", format = "tiff", density = 300, compression = "LZW")

message("Wrote km_panel_polypharmacy_polished.pdf, .png, .tiff (styling matches main analysis script).")
