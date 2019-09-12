source("libSVPaper.R")
source("anchorSupportSummary.R")
source("SVPaper_manta_gridss_consistency.R")
source("SVPaperLine.R")

plot_gridss_top = grid.arrange(
  plot_assembly_phasing,
  #plot_sgllongrepeat,
  plot_manta_vs_gridss_tp_colo829,
  plot_manta_vs_gridss_fp_colo829,
  plot_manta_vs_gridss_fp_na12878,
  probe_results_vs_manta_over_50bp,
  plot_probe_results_vs_stelka_under_50bp,
  widths = c(5, 1, 1),
  heights = c(1.2,1,1,0.2,0.5,1.5), # TODOO adjust heights so legend can go at bottom
  layout_matrix = rbind(
    c(1, 2, 2),
    c(1, 3, 3),
    c(1, 4, 4),
    c(1, NA, NA),
    c(1, 5, 6),
    c(NA, 5, 6)))
plot_gridss_top
figsave("figure2_gridss_top", plot_gridss_top, width=10, height=8)

