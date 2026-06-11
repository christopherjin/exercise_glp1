# combine_figure_panels.R
# ------------------------------------------------------------------------------
# Assembles the individual panel files in figures/final_figs/ into the composite
# multi-panel manuscript figures, following the panel labeling annotated in
# glp1a_exploration.Rmd:
#
#   Figure 1 = Fig_1A (schematic) + Fig_1B (detection summary table)
#              + Fig_1C (rat/human detection grid)
#   Figure 2 = Fig_2A (per-tissue z-score heatmap) + Fig_2B (triangle direction)
#              + Fig_2C (Cckbr HIPPOC epigenetic/transcript)
#   Figure 3 = Fig_3 (GIP/GLP1 immunoassay heatmaps)
#
# Each source PDF is rasterized at high density and each grid cell is sized to
# the panel's native aspect ratio, so panels are reproduced sharply with no
# stretching and no overlap. Outputs land in figures/composite_figures/ as PDFs.
# ------------------------------------------------------------------------------

suppressMessages({
  library(magick)
  library(cowplot)
  library(here)
})

panel_dir = here("figures", "final_figs")
out_dir = here("figures", "composite_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

target_long_px = 10800   # ~12 in at 900 dpi: pixels along a full-width panel
max_density = 1200       # ceiling on PDF rasterization density (memory guard)
label_size = 22

# Read a panel (PDF or raster), flatten onto white, and rasterize PDFs at a
# density chosen so the longest side lands near target_long_px. This keeps small
# panels crisp at 900 dpi without rendering very wide heatmaps at absurd sizes.
read_panel = function(file) {
  path = file.path(panel_dir, file)
  if (!file.exists(path)) stop("Missing panel: ", path)
  is_pdf = grepl("\\.pdf$", file, ignore.case = TRUE)
  if (is_pdf) {
    inf72 = image_info(image_read(path, density = 72))[1, ]   # page size: px/72 = inches
    long_in = max(inf72$width, inf72$height) / 72
    density = min(max_density, ceiling(target_long_px / long_in))
    img = image_read(path, density = density)
  } else {
    img = image_read(path)
  }
  image_flatten(image_background(img, "white"))
}

# Aspect ratio (width / height) of a magick image.
panel_ar = function(img) {
  inf = image_info(img)[1, ]
  inf$width / inf$height
}

# Wrap a magick image as a ggdraw object for use inside plot_grid.
as_panel = function(img) ggdraw() + draw_image(img)

# ---- Figure 1 ----------------------------------------------------------------
# Top row: A (schematic) beside B (summary table), matched on height.
# Bottom row: C spanning the full width.
p1a = read_panel("Fig_1A.png")
p1b = read_panel("Fig_1B_detection_summary_table.png")
p1c = read_panel("Fig_1C_combined_split_detection_table.pdf")

ar_1a = panel_ar(p1a)
ar_1b = panel_ar(p1b)
ar_1c = panel_ar(p1c)

canvas_w1 = 12
top_row_h = canvas_w1 / (ar_1a + ar_1b)   # shared height that fits A and B
bottom_row_h = canvas_w1 / ar_1c          # full-width C
canvas_h1 = top_row_h + bottom_row_h

fig1_top = plot_grid(
  as_panel(p1a), as_panel(p1b),
  nrow = 1,
  rel_widths = c(ar_1a, ar_1b),
  labels = c("A", "B"),
  label_size = label_size
)
fig1 = plot_grid(
  fig1_top, as_panel(p1c),
  ncol = 1,
  rel_heights = c(top_row_h, bottom_row_h),
  labels = c("", "C"),
  label_size = label_size
)

save_plot(file.path(out_dir, "Figure_1.pdf"), fig1,
          base_width = canvas_w1, base_height = canvas_h1)
message("Wrote Figure_1 (", round(canvas_w1, 2), " x ", round(canvas_h1, 2), " in)")

# ---- Figure 2 ----------------------------------------------------------------
# Top: A heatmap full width. Bottom: B and C side by side (equal aspect ratios).
p2a = read_panel("Fig_2A_all_tissue_combined_rna_zscore.pdf")
p2b = read_panel("Fig_2B_rna_triangle_direction.pdf")
p2c = read_panel("Fig_2C_Cckbr_HIPPOC.pdf")

ar_2a = panel_ar(p2a)
ar_2b = panel_ar(p2b)
ar_2c = panel_ar(p2c)

canvas_w2 = 12
row1_h2 = canvas_w2 / ar_2a               # full-width A
row2_h2 = canvas_w2 / (ar_2b + ar_2c)     # B + C share the width
canvas_h2 = row1_h2 + row2_h2

fig2_bottom = plot_grid(
  as_panel(p2b), as_panel(p2c),
  nrow = 1,
  rel_widths = c(ar_2b, ar_2c),
  labels = c("B", "C"),
  label_size = label_size
)
fig2 = plot_grid(
  as_panel(p2a), fig2_bottom,
  ncol = 1,
  rel_heights = c(row1_h2, row2_h2),
  labels = c("A", ""),
  label_size = label_size
)

save_plot(file.path(out_dir, "Figure_2.pdf"), fig2,
          base_width = canvas_w2, base_height = canvas_h2)
message("Wrote Figure_2 (", round(canvas_w2, 2), " x ", round(canvas_h2, 2), " in)")

# ---- Figure 3 ----------------------------------------------------------------
# Single combined immunoassay panel; rendered at full resolution for the set.
p3 = read_panel("Fig_3_gip_glp_immunoassay.pdf")
ar_3 = panel_ar(p3)
canvas_w3 = 9
canvas_h3 = canvas_w3 / ar_3
fig3 = as_panel(p3)

save_plot(file.path(out_dir, "Figure_3.pdf"), fig3,
          base_width = canvas_w3, base_height = canvas_h3)
message("Wrote Figure_3 (", round(canvas_w3, 2), " x ", round(canvas_h3, 2), " in)")

message("Composite figures written to ", out_dir)
