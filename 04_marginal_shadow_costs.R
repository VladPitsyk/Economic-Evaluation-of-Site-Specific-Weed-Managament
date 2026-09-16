# =============================================================================
# Marginal shadow cost of one additional weed plant.
#
#
# Description:
#     Produces the four-panel figure of marginal shadow costs against local
#     weed density. The marginal shadow cost is the increase in the discounted
#     ten-year weed liability caused by adding one weed plant to a grid cell
#     before the Year-0 herbicide application, holding that cell's spray
#     decision fixed. It is computed per cell by 02_bioeconomic_simulation.py
#     and written to the grid layer as Marg_SC_D (dicots) and Marg_SC_M
#     (monocots).
#
#     Panels are split by weed class and by whether the cell was sprayed. The
#     point cloud shows individual grid cells pooled across fields; the solid
#     line is an inverse-polynomial fit to the binned maximum, which traces the
#     declining upper boundary of the cloud rather than its centre. The upper
#     envelope is used because the quantity of interest is how large the
#     marginal cost can be at a given density, not its average.
#
# Inputs:
#     A folder of CSV exports of the grid layer, one per field, each containing
#     at least the columns:
#         Spray_Dec  "Spray" or "No", the Year-0 decision for that cell
#         P_Dicot    detected dicot count in the cell
#         P_Monocot  detected monocot count in the cell
#         TrueDenD   estimated true pre-treatment dicot density (plants m-2)
#         TrueDenM   estimated true pre-treatment monocot density (plants m-2)
#         Marg_SC_D  marginal shadow cost of one additional dicot (EUR, 10-yr NPV)
#         Marg_SC_M  marginal shadow cost of one additional monocot (EUR, 10-yr NPV)
#
#     IMPORTANT: the grid layer carries no column identifying which decision
#     strategy produced it, so INPUT_FOLDER must contain exports from a single
#     strategy at a single resolution only. The figure in the paper uses the
#     economic threshold (ET) strategy at 1 x 1 m across all 19 fields. Mixing
#     strategies would pool cells treated at different doses, and the dual-dose
#     strategy has no untreated cells at all, so the lower panels would be
#     empty for it.
#
# Outputs:
#     - Marginal_shadow_costs.png, a 2 x 2 panel figure at 300 dpi.
#
# Requirements:
#     R >= 4.2 with ggplot2, dplyr, gridExtra, purrr, readr.
#
# Usage:
#     Set INPUT_FOLDER and OUTPUT_FILE below, then source the script.
# =============================================================================

library(ggplot2)
library(dplyr)
library(gridExtra)
library(purrr)
library(readr)


# =============================================================================
# USER CONFIGURATION
# =============================================================================

# Folder holding the per-field CSV exports of the grid layer. Every .csv in
# this folder is read and combined, so it must contain exports from one
# decision strategy and one resolution only (see the note in the header).
INPUT_FOLDER <- "C:/path/to/marginal_shadow_cost_exports"

# Output figure.
OUTPUT_FILE <- "C:/path/to/output/Marginal_shadow_costs.png"

# Number of fields expected in INPUT_FOLDER. Set to NA to skip the check.
EXPECTED_N_FILES <- 19

# Axis limits. The x range is restricted to the density interval containing
# most observations, where the marginal cost changes most steeply; the sprayed
# panels cover a wider density range than the unsprayed ones because sprayed
# cells are the denser ones by construction.
max_x_top      <- 60   # sprayed cells, both classes
max_x_bottom_d <- 10   # unsprayed cells, dicots
max_x_bottom_m <- 5    # unsprayed cells, monocots


# =============================================================================
# 1. IMPORT
# =============================================================================

file_list <- list.files(path = INPUT_FOLDER, pattern = "\\.csv$", full.names = TRUE)

if (length(file_list) == 0)
  stop("No CSV files found in INPUT_FOLDER: ", INPUT_FOLDER)

if (!is.na(EXPECTED_N_FILES) && length(file_list) != EXPECTED_N_FILES)
  warning(sprintf(paste("Found %d CSV files but expected %d. Check that the folder",
                        "holds one export per field, from a single strategy and",
                        "resolution."),
                  length(file_list), EXPECTED_N_FILES))

message("Reading ", length(file_list), " files...")
df_master <- map_dfr(file_list, read_csv, show_col_types = FALSE)

required <- c("Spray_Dec", "P_Dicot", "P_Monocot",
              "TrueDenD", "TrueDenM", "Marg_SC_D", "Marg_SC_M")
missing_cols <- setdiff(required, names(df_master))
if (length(missing_cols) > 0)
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

message("Loaded ", nrow(df_master), " grid cells.")


# =============================================================================
# 2. PREPARE THE POINT CLOUDS
# =============================================================================

# Densities are the estimated TRUE pre-treatment densities, i.e. after the
# recall correction and the spatial allocation of undetected weeds, so the
# x axis reflects the weeds actually present rather than those detected.
df_clean <- df_master %>%
  mutate(Dens_D = TrueDenD,
         Dens_M = TrueDenM)

# Cells with a zero marginal cost carry no information about the shape of the
# relationship and are excluded, as are cells with no weeds of the class being
# plotted. The sprayed panels additionally require a detection in that class,
# since a cell is sprayed on the basis of what was detected.
df_ds <- df_clean %>%
  filter(Spray_Dec == "Spray", P_Dicot   > 0,
         Dens_D > 0, Dens_D <= max_x_top,      Marg_SC_D > 0)
df_ms <- df_clean %>%
  filter(Spray_Dec == "Spray", P_Monocot > 0,
         Dens_M > 0, Dens_M <= max_x_top,      Marg_SC_M > 0)
df_dn <- df_clean %>%
  filter(Spray_Dec == "No",
         Dens_D > 0, Dens_D <= max_x_bottom_d, Marg_SC_D > 0)
df_mn <- df_clean %>%
  filter(Spray_Dec == "No",
         Dens_M > 0, Dens_M <= max_x_bottom_m, Marg_SC_M > 0)

message(sprintf("Cells plotted: dicot/sprayed %d, monocot/sprayed %d, dicot/unsprayed %d, monocot/unsprayed %d",
                nrow(df_ds), nrow(df_ms), nrow(df_dn), nrow(df_mn)))


# =============================================================================
# 3. UPPER ENVELOPES
# =============================================================================

# The maximum marginal cost within each density bin. Bin widths are narrower
# for the unsprayed panels because those cover a shorter density range.
upper_envelope <- function(df, dens_col, cost_col, bin_per_unit) {
  df %>%
    mutate(x_bin = round(.data[[dens_col]] * bin_per_unit) / bin_per_unit) %>%
    group_by(x_bin) %>%
    summarise(y_max = max(.data[[cost_col]], na.rm = TRUE), .groups = "drop")
}

env_ds <- upper_envelope(df_ds, "Dens_D", "Marg_SC_D",  5)
env_ms <- upper_envelope(df_ms, "Dens_M", "Marg_SC_M",  5)
env_dn <- upper_envelope(df_dn, "Dens_D", "Marg_SC_D", 10)
env_mn <- upper_envelope(df_mn, "Dens_M", "Marg_SC_M", 20)


# =============================================================================
# 4. PANELS
# =============================================================================

pub_theme <- theme_minimal() +
  theme(
    plot.title   = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 16),
    axis.title.x = element_text(face = "bold", size = 16),
    axis.text    = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
  )

x_label <- bquote(bold("Estimated true pre-treatment density (plants/" ~ m^2 ~ ")"))

# The envelope is fitted as a linear model in 1/(x + c) and 1/(x + c)^2, which
# gives a monotonically declining curve that flattens at high density without
# imposing a functional form on the underlying process. The offset c prevents
# the fit from diverging as density approaches zero and is set to roughly one
# bin width, so it is smaller for the panels with the narrower density range.
make_panel <- function(df, env, dens_col, cost_col, title, x_max, y_max,
                       offset, x_breaks, formula_rhs) {
  ggplot(df, aes(x = .data[[dens_col]], y = .data[[cost_col]])) +
    geom_point(alpha = 0.08, size = 1.2, color = "#555555") +
    geom_smooth(data = env, aes(x = x_bin, y = y_max),
                method = "lm", formula = formula_rhs,
                se = FALSE, color = "black", linewidth = 1.2) +
    coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
    scale_x_continuous(n.breaks = x_breaks) +
    scale_y_continuous(n.breaks = 8) +
    labs(title = title, x = x_label, y = "Marginal Shadow Cost (\u20ac/weed)") +
    pub_theme
}

p_dicot_spray <- make_panel(
  df_ds, env_ds, "Dens_D", "Marg_SC_D", "Dicots in sprayed cells",
  max_x_top, 0.032, 0.5, 8,
  y ~ I(1 / (x + 0.5)) + I(1 / (x + 0.5)^2))

p_monocot_spray <- make_panel(
  df_ms, env_ms, "Dens_M", "Marg_SC_M", "Monocots in sprayed cells",
  max_x_top, 0.012, 0.5, 8,
  y ~ I(1 / (x + 0.5)) + I(1 / (x + 0.5)^2))

p_dicot_no <- make_panel(
  df_dn, env_dn, "Dens_D", "Marg_SC_D", "Dicots in non-sprayed cells",
  max_x_bottom_d, 0.075, 0.1, 6,
  y ~ I(1 / (x + 0.1)) + I(1 / (x + 0.1)^2))

# The monocot/unsprayed panel is fitted through the origin (0 +) because the
# few observations at very low density would otherwise dominate the intercept.
p_monocot_no <- make_panel(
  df_mn, env_mn, "Dens_M", "Marg_SC_M", "Monocots in non-sprayed cells",
  max_x_bottom_m, 0.075, 0.03, 6,
  y ~ 0 + I(1 / (x + 0.03)) + I(1 / (x + 0.03)^2))


# =============================================================================
# 5. ASSEMBLE AND SAVE
# =============================================================================

message("Assembling the 2 x 2 panel figure...")

final_grid <- arrangeGrob(
  p_dicot_spray, p_monocot_spray,
  p_dicot_no,    p_monocot_no,
  ncol = 2
)

grid::grid.newpage()
grid::grid.draw(final_grid)

ggsave(
  filename = OUTPUT_FILE,
  plot     = final_grid,
  width    = 12,
  height   = 10,
  dpi      = 300,
  bg       = "white"
)

message("Figure written to: ", OUTPUT_FILE)
