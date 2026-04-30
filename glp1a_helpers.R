conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("group_by", "dplyr")
conflict_prefer("mutate", "dplyr")
conflict_prefer("summarise", "dplyr")
conflict_prefer("arrange", "dplyr")
conflict_prefer("intersect", "dplyr")
conflict_prefer("setdiff", "dplyr")

find_interested_features = function(tissue_types_1b){
  all_da_results_subset_filter = data.frame()
  for(tissue in tissue_types_1b){
    all_da_results_subset = MotrpacRatTraining6mo::combine_da_results(tissues = tissue, assays = c("ATAC", "METHYL", "TRNSCRPT"), include_epigen = T)
    # all_da_results_subset = MotrpacRatTraining6mo::combine_da_results(tissues = tissue, assays = c("ATAC", "METHYL", "TRNSCRPT"), include_epigen = F)
    all_da_results_subset_filter_single = merge(all_da_results_subset, all_feat_in, by = 'feature_ID') %>%
      dplyr::select(feature_ID, gene_symbol, assay, tissue, sex, comparison_group, p_value, adj_p_value, logFC, zscore) %>%
      dplyr::group_by(gene_symbol)
    all_da_results_subset_filter = rbind(all_da_results_subset_filter, all_da_results_subset_filter_single)
    gc()
  }
  all_da_results_subset_filter = all_da_results_subset_filter[!duplicated(all_da_results_subset_filter), ]
  return(all_da_results_subset_filter)
}

get_expression_stats = function(desired_feature,
                                tissue_list,
                                assay_type = c("TRNSCRPT", "IMMUNO"),
                                group_filter = "control",
                                id_column = "feature_ID",
                                sample_id_column = c("viallabel", "pid")) {

  assay_type = match.arg(assay_type)
  sample_column = ifelse(assay_type == "TRNSCRPT", sample_id_column[1], sample_id_column[2])

  all_counts = list()

  for (tissue in tissue_list) {
    ome_actual_data = load_sample_data(tissue = tissue, assay = assay_type, normalized = (assay_type != "TRNSCRPT"))

    if (nrow(ome_actual_data) == 0) next

    rownames(ome_actual_data) = NULL
    ome_actual_data = ome_actual_data %>%
      distinct(feature_ID, .keep_all = TRUE) %>%
      tibble::column_to_rownames(id_column) %>%
      select(!matches("feature|tissue|assay|dataset"))

    relevant_samples = MotrpacRatTraining6moData::PHENO %>%
      filter(.data[[sample_column]] %in% colnames(ome_actual_data)) %>%
      filter(group == group_filter) %>%
      distinct(.data[[sample_column]], .keep_all = TRUE)

    selected_data = ome_actual_data %>%
      select(as.character(relevant_samples[[sample_column]]))

    if (assay_type == "TRNSCRPT") {
      selected_data = edgeR::cpm(selected_data)
    }

    avg = rowMeans(selected_data)
    sdv = apply(selected_data, 1, sd)

    if (!desired_feature %in% names(avg)) next

    all_counts[[tissue]] = data.frame(
      Tissue = tissue,
      Average = avg[[desired_feature]],
      SD = sdv[[desired_feature]]
    )
  }

  counts_df = do.call(rbind, all_counts)
  rownames(counts_df) = NULL
  return(counts_df)
}

#slightly modified function to get stats by sex. mostly copied over from other function
get_expression_stats_by_sex = function(desired_feature,
                                       tissue_list,
                                       assay_type = c("TRNSCRPT", "IMMUNO"),
                                       group_filter = "control",
                                       id_column = "feature_ID",
                                       sample_id_column = c("viallabel", "pid")) {

  assay_type = match.arg(assay_type)
  sample_column = ifelse(assay_type == "TRNSCRPT", sample_id_column[1], sample_id_column[2])

  all_counts = list()

  for (tissue in tissue_list) {
    ome_actual_data = load_sample_data(tissue = tissue, assay = assay_type, normalized = (assay_type != "TRNSCRPT"))

    if (nrow(ome_actual_data) == 0) next

    rownames(ome_actual_data) = NULL
    ome_actual_data = ome_actual_data %>%
      distinct(feature_ID, .keep_all = TRUE) %>%
      tibble::column_to_rownames(id_column) %>%
      select(!matches("feature|tissue|assay|dataset"))

    relevant_samples = MotrpacRatTraining6moData::PHENO %>%
      filter(.data[[sample_column]] %in% colnames(ome_actual_data)) %>%
      filter(group == group_filter) %>%
      distinct(.data[[sample_column]], .keep_all = TRUE)

    selected_data = ome_actual_data %>%
      select(as.character(relevant_samples[[sample_column]]))

    if (assay_type == "TRNSCRPT") {
      selected_data = edgeR::cpm(selected_data)
    }

    sex_vec = relevant_samples$sex
    names(sex_vec) = relevant_samples[[sample_column]]

    avg = sapply(split(names(sex_vec), sex_vec),
                 function(s) rowMeans(selected_data[, s, drop = FALSE]))
    sdv = sapply(split(names(sex_vec), sex_vec),
                 function(s) apply(selected_data[, s, drop = FALSE], 1, sd))

    if (!desired_feature %in% rownames(avg)) next

    all_counts[[tissue]] =
      data.frame(
        Tissue = tissue,
        Sex = colnames(avg),
        Average = avg[desired_feature, ],
        SD = sdv[desired_feature, ]
      )
  }

  counts_df = do.call(rbind, all_counts)
  rownames(counts_df) = NULL
  return(counts_df)
}

#needs breaks because the CPM dist for GLP1R is so different
#immuno breaks not needed for the main 2 but useful for if i need to make others
plot_rna_immuno_combined = function(rna_data,
                                    immuno_data,
                                    gene_name,
                                    include_first_legend = TRUE,
                                    rna_y_breaks = NULL,
                                    immuno_y_breaks = NULL,
                                    immuno_reverse = TRUE) {
  merged_data = full_join(rna_data, immuno_data, by = "Tissue")
  merged_data$Tissue = factor(merged_data$Tissue, levels = names(MotrpacBicQC::tissue_cols))
  #tissues run in transcript but not immuno
  merged_data$missing_immuno = is.na(merged_data$Average.y)
  #so here - the luminex assays were run in 17 tissues, not run in 2.
  #in some of the tissues, the feature was not detected -> is na
  merged_data[merged_data$Tissue %in% c("HYPOTH", "VENACV"),]$missing_immuno = "Not run"
  merged_data$missing_transcript = merged_data$Average.x < 0.5

  if(gene_name == "GLP1R") y_adjust = 10
  if(gene_name == "GIPR") y_adjust = 4

  rna_plot = ggplot(merged_data, aes(x = Tissue, y = Average.x, fill = Tissue)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_errorbar(aes(ymin = Average.x - SD.x, ymax = Average.x + SD.x), width = 0.3) +
    geom_text(
      data = merged_data %>% filter(missing_transcript),
      aes(x = Tissue, y = y_adjust, label = "< 0.2 CPM."),
      inherit.aes = FALSE,
      size = 2.5,
      angle = 90,
      vjust = -0.3
    ) +
    scale_fill_manual(values = MotrpacBicQC::tissue_cols) +
    ggtitle(paste("RNA-seq Expression -", gene_name)) +
    ylab("Counts per million (CPM)") +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = if (include_first_legend) "right" else "none"
    )

  if (!is.null(rna_y_breaks)) {
    for (brk in rna_y_breaks) {
      rna_plot = rna_plot +
        ggbreak::scale_y_zoom(brk, zoom = 0.5)
    }
  }

  #----make sure its putting the metabolite level not the gene name-------
  signalling_molecule_name = sub("R$", "", gene_name)
  #---this just pops off the R at the end to remove the receptor----------

  immuno_plot = ggplot(merged_data, aes(x = Tissue, y = Average.y, fill = Tissue)) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_errorbar(aes(ymin = Average.y - SD.y, ymax = Average.y + SD.y), width = 0.3) +
    geom_text(
      data = merged_data %>% filter(missing_immuno == "TRUE"),
      aes(x = Tissue, y = 5, label = "N.D."),
      inherit.aes = FALSE,
      size = 2.5,
      angle = 0,
      vjust = -0.3
    ) +
    geom_text(
      data = merged_data %>% filter(missing_immuno == "Not run"),
      aes(x = Tissue, y = 5, label = "N.R."),
      inherit.aes = FALSE,
      size = 2.5,
      angle = 0,
      vjust = -0.3
    ) +
    scale_fill_manual(values = MotrpacBicQC::tissue_cols) +
    ylab("Immunoassay\n(Values relative to plate standards)") +
    ggtitle(paste("Immunoassay Expression -", signalling_molecule_name)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    )

  if (!is.null(immuno_y_breaks)) {
    for (brk in immuno_y_breaks) {
      immuno_plot = immuno_plot +
        ggbreak::scale_y_break(brk, scales = "fixed", space = 0.03)
    }
    #if we want to make the directionality go down, like if its kinda split
    #in the middle. I thought it would look better than it does.
  }
  if (immuno_reverse) {
    immuno_plot = immuno_plot + scale_y_reverse()
  }
  return((rna_plot / immuno_plot) + plot_layout(heights = c(1.5, 1), guides = "collect"))
}



plot_rna_immuno_jitterbysex = function(rna_data,
                                       immuno_data = NULL,
                                       gene_name,
                                       include_first_legend = TRUE,
                                       rna_y_breaks = NULL,
                                       immuno_y_breaks = NULL,
                                       immuno_reverse = TRUE) {
  has_immuno = !is.null(immuno_data)

  if (has_immuno) {
    merged_data = full_join(rna_data, immuno_data, by = c("Tissue", "Sex"))
    merged_data$missing_immuno = is.na(merged_data$Average.y)
    merged_data[merged_data$Tissue %in% c("HYPOTH", "VENACV"),]$missing_immuno = "Not run"
  } else {
    merged_data = rna_data %>% rename(Average.x = Average, SD.x = SD)
    merged_data$Average.y = NA_real_
    merged_data$SD.y = NA_real_
  }

  merged_data$Tissue = factor(merged_data$Tissue, levels = names(MotrpacBicQC::tissue_cols))
  merged_data$missing_transcript = merged_data$Average.x < 0.5

  y_adjust = 0
  if (gene_name == "GLP1R") y_adjust = 5
  if (gene_name == "GIPR") y_adjust = 4

  # female = solid fill, male = crosshatch overlay to distinguish sex within tissue color
  sex_patterns = c("female" = "none", "male" = "crosshatch")

  rna_plot = ggplot(merged_data, aes(x = Tissue, y = Average.x, fill = Tissue, group = Sex)) +
    ggpattern::geom_bar_pattern(
      aes(pattern = Sex),
      stat = "identity",
      width = 0.7,
      position = position_dodge(width = 0.8),
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_density = 0.15,
      pattern_spacing = 0.025,
      pattern_angle = 45,
      colour = "grey30",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(ymin = Average.x - SD.x, ymax = Average.x + SD.x),
      width = 0.25,
      position = position_dodge(width = 0.8)
    ) +
    scale_fill_manual(values = MotrpacBicQC::tissue_cols) +
    ggpattern::scale_pattern_manual(
      values = sex_patterns,
      labels = c("female" = "Female (solid)", "male" = "Male (crosshatch)"),
      name = "Sex"
    ) +
    ggtitle(paste("RNA-seq Expression -", gene_name)) +
    ylab("Counts per million (CPM)") +
    theme_bw() +
    theme(
      axis.text.x = if (has_immuno) element_blank() else element_text(size = 9, angle = 45, hjust = 1, color = "black"),
      axis.title.x = element_blank(),
      axis.ticks.x = if (has_immuno) element_blank() else element_line(),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = if (include_first_legend) "right" else "none"
    )

  if (!is.null(rna_y_breaks))
    # rna_plot = rna_plot + ggbreak::scale_y_break(rna_y_breaks, scales = "fixed", space = 0.1)
    rna_plot = rna_plot + ggforce::facet_zoom(y = Average.x < rna_y_breaks[1], zoom.size = 2)
  # ggforce::facet_zoom(y = Average.x < rna_y_breaks[1], zoom.size = 1)

  if (!has_immuno) return(rna_plot)

  #----make sure its putting the metabolite level not the gene name-------
  signalling_molecule_name = sub("R$", "", gene_name)
  #---this just pops off the R at the end to remove the receptor----------

  immuno_plot = ggplot(merged_data, aes(x = Tissue, y = Average.y, fill = Tissue, group = Sex)) +
    ggpattern::geom_bar_pattern(
      aes(pattern = Sex),
      stat = "identity",
      width = 0.7,
      position = position_dodge(0.8),
      pattern_fill = "black",
      pattern_colour = "black",
      pattern_density = 0.15,
      pattern_spacing = 0.025,
      pattern_angle = 45,
      colour = "grey30",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(ymin = Average.y - SD.y, ymax = Average.y + SD.y),
      width = 0.25, position = position_dodge(0.8)
    ) +
    geom_text(
      data = merged_data %>% filter(missing_immuno == "TRUE"),
      aes(x = Tissue, y = 5, label = "N.D.", group = Sex),
      # position = position_dodge(0.8),
      inherit.aes = FALSE, size = 2.5, vjust = -0.3
    ) +
    geom_text(
      data = merged_data %>% filter(missing_immuno == "Not run"),
      aes(x = Tissue, y = 5, label = "N.R.", group = Sex),
      # position = position_dodge(0.8),
      inherit.aes = FALSE, size = 2.5, vjust = -0.3
    ) +
    scale_fill_manual(values = MotrpacBicQC::tissue_cols) +
    ggpattern::scale_pattern_manual(
      values = sex_patterns,
      labels = c("female" = "Female (solid)", "male" = "Male (crosshatch)"),
      name = "Sex"
    ) +
    ylab("Immunoassay\n(Luminex Normalized Protein Expression)") +
    ggtitle(paste("Immunoassay Expression -", signalling_molecule_name)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    )

  if (!is.null(immuno_y_breaks))
    immuno_plot = immuno_plot + ggbreak::scale_y_break(immuno_y_breaks, scales = "fixed", space = 0.03)

  #if we want to make the directionality go down, like if its kinda split
  #in the middle. I thought it would look better than it does.
  if (immuno_reverse) {
    immuno_plot = immuno_plot + scale_y_reverse()
  }
  return((rna_plot / immuno_plot) + plot_layout(heights = c(1.5, 1), guides = "collect"))
}


make_a_plot_for_tissue_heatmaps = function(select_tissue_types,
                                           interested_genes,
                                           value_col = "zscore",
                                           cutoff = 0.1,
                                           show_heatmap_legend = TRUE,
                                           legend_tissue_types = select_tissue_types) {
  combined_ome_metadata = list(); combined_ome_fdr = list(); current_ome = "TRNSCRPT"

  for (tissue_value in select_tissue_types){
    da_res_check = all_da_results_subset_filter %>%
      dplyr::filter(tissue == tissue_value) %>%
      dplyr::filter(assay == current_ome) %>%
      dplyr::mutate(adj_p_value = p.adjust(p_value, method = "BH"), .by = c("tissue", "assay", "comparison_group"))

    if(nrow(da_res_check) == 0){next} #skip if no values

    all_genes = unique(all_da_results_subset_filter$gene_symbol)
    genes_found_in_tissue = unique(da_res_check$gene_symbol)
    genes_not_found_in_tissue = setdiff(all_genes, genes_found_in_tissue)

    expected_groups = c("1w", "2w", "4w", "8w")
    da_res_tissue = da_res_check %>%
      rename(group = comparison_group)  %>%
      group_by(sex, tissue) %>%  # Group only by sex and tissue
      tidyr::complete(gene_symbol = all_genes, group = expected_groups, fill = setNames(list(NaN), value_col)) %>%
      ungroup() %>%
      select(gene_symbol, assay, tissue, group, zscore, logFC, adj_p_value, sex) %>%
      mutate(assay = "TRNSCRPT")

    make_ome_metadata_matrix = function(value_col, group_levels = c("1w", "2w", "4w", "8w")) {

      MotrpacRatTraining6moData::PHENO %>%
        select(sex, group, tissue) %>%
        distinct() %>%
        right_join(da_res_tissue, by = c("sex", "group", "tissue")) %>%
        mutate(group = factor(group, levels = group_levels)) %>%
        arrange(sex, group) %>%
        distinct(sex, group, assay, tissue, gene_symbol, .keep_all = TRUE) %>%
        rename(value = !!value_col) %>%
        group_by(sex, group, assay, tissue) %>%
        pivot_wider(
          id_cols = c(sex, group, tissue),
          values_from = value,
          names_from = gene_symbol
        ) %>%
        as.data.frame() %>%
        `rownames<-`(with(., interaction(sex, group, tissue))) %>%
        select(-c(sex, group, tissue)) %>%
        t() %>%
        as.matrix()
    }

    ome_metadata_values = make_ome_metadata_matrix(value_col)
    ome_metadata_fdr = make_ome_metadata_matrix("adj_p_value")

    combined_ome_metadata[[tissue_value]] <- ome_metadata_values
    combined_ome_fdr[[tissue_value]] <- ome_metadata_fdr
  }
  common_genes <- Reduce(intersect, lapply(combined_ome_metadata, rownames))

  filtered_ome_metadata <- lapply(combined_ome_metadata, function(mat) {
    mat[common_genes, , drop = FALSE]
  })
  final_ome_metadata <- do.call(cbind, filtered_ome_metadata)

  #------again the fdr is only for the asterisks-------
  filtered_fdr <- lapply(combined_ome_fdr, function(mat) {
    mat[common_genes, , drop = FALSE]
  })
  filtered_fdr <- do.call(cbind, filtered_fdr)
  #-----------------------------------------------

  annotation_df = colnames(final_ome_metadata) %>%
    strsplit(split = "\\.") %>%
    do.call(rbind, .) %>%
    as.data.frame(stringsAsFactors = FALSE)

  #we want to order the legend in the same way as the columns
  colnames(annotation_df) = c("sex", "group", "tissue")
  color_order = annotation_df %>% pull(tissue) %>% unique()

  annotation_df = annotation_df %>%
    mutate(Sex = sex,
           Timepoint = group,
           Tissue = factor(tissue,
                           levels = color_order))  # Keep these for plotting

  combined_annotation = ComplexHeatmap::HeatmapAnnotation(
    Tissue = ComplexHeatmap::anno_block(
      labels = color_order,
      labels_gp = gpar(fontsize = 6 * scale, fontface = "bold"),
      gp = gpar(fill = NA, col = NA),
      height = unit(8, "pt") * scale
    ),
    Sex = annotation_df$Sex,
    Timepoint = annotation_df$Timepoint,
    border = TRUE,
    gp = gpar(col = "black"),
    gap = 0,
    which = "column",
    show_legend = FALSE,
    simple_anno_size = unit(6, "pt") * scale,
    col = list(
      Sex = c("female" = "#ff6eff", "male" = "#5555ff"),
      Timepoint = c('1w' = '#F7FCB9',
                    '2w' = '#ADDD8E',
                    '4w' = '#238443',
                    '8w' = '#002612')
    ),
    annotation_name_gp = gpar(fontsize = 7 * scale),
    annotation_legend_param = list(Timepoint = list(
      at = c("1w", "2w", "4w", "8w")),
      border = "black",
      labels_gp = gpar(fontsize = 6.5 * scale),
      title_gp = gpar(fontsize = 7 * scale, fontface = "bold")
    )
  )
  if (value_col == "zscore") {
    col_breaks = c(-1.6, 0, 1.8)
    legend_title = "Z-Score \nper Group/Sex"
    legend_at = c(-3, -1:1, 3)
  } else {
    col_breaks = c(-1, 0, 1)
    legend_title = "logFC VS\nsex-matched control"
    legend_at = c(-5, -3, 0, 3, 5)
  }

  col_fun = circlize::colorRamp2(breaks = col_breaks, colors = c("#3366ff", "white", "darkred"))

  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    show_heatmap_legend = FALSE,
    matrix = final_ome_metadata,
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    clustering_distance_rows = "pearson",
    top_annotation = combined_annotation,
    border = "black",
    row_names_gp = gpar(fontsize = 5 * scale),
    height = nrow(final_ome_metadata) * unit(5.5, "pt") * scale,
    width = ncol(final_ome_metadata) * unit(5.5, "pt") * scale,
    column_split = rep(seq_along(color_order), each = 8),
    column_title = NULL,
    cell_fun = function(j, i, x, y, width, height, fill) {
      #Heatmap grid
      grid.rect(x = x, y = y, width, height,
                gp = gpar(col = "#555555"))
      if ((j - 1) %% 8 == 3) {
        grid.lines(c(x + width/2, x + width/2),
                   c(y - height/2, y + height/2),
                   gp = gpar(col = "black", lwd = 5))
      }
      #update to use FDR instead of easy Zscore filt
      if(!is.na(filtered_fdr[i, j]) && abs(filtered_fdr[i, j]) < cutoff) {
        gb = textGrob("*")
        gb_w = convertWidth(grobWidth(gb), "mm")
        gb_h = convertHeight(grobHeight(gb), "mm")
        grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
      }
    }
  )

  if (!show_heatmap_legend) {
    return(list(ht = ht, lgd = NULL))
  }

  lgd_sex = ComplexHeatmap::Legend(
    labels = c("female", "male"),
    legend_gp = gpar(fill = c("female" = "#ff6eff", "male" = "#5555ff")),
    title = "Sex",
    title_gp = gpar(fontsize = 7 * scale, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5 * scale),
    border = "black"
  )
  lgd_timepoint = ComplexHeatmap::Legend(
    labels = c("1w", "2w", "4w", "8w"),
    legend_gp = gpar(fill = c('#F7FCB9', '#ADDD8E', '#238443', '#002612')),
    title = "Timepoint",
    title_gp = gpar(fontsize = 7 * scale, fontface = "bold"),
    labels_gp = gpar(fontsize = 6.5 * scale),
    border = "black"
  )
  lgd_color = ComplexHeatmap::Legend(
    col_fun = col_fun,
    at = legend_at,
    title = legend_title,
    title_gp = gpar(fontsize = 7 * scale, fontface = "bold"),
    labels_gp = gpar(fontsize = 6 * scale),
    legend_height = 5 * scale * unit(8, "pt"),
    border = "black"
  )

  lgd = ComplexHeatmap::packLegend(lgd_sex, lgd_timepoint, lgd_color,
                                   direction = "vertical", gap = unit(3, "mm"))

  return(list(ht = ht, lgd = lgd))
}

make_immuno_heatmap = function(desired_feature, color_range = NULL, show_legend = TRUE, show_anno_names = TRUE){
  all_immuno_da_results = MotrpacRatTraining6mo::combine_da_results(tissues = immuno_tissue_types,
                                                                    assays = "IMMUNO",
                                                                    include_epigen = F) %>%
    filter(feature_ID %in% c("GIP", "GLP1", "GLUCAGON")) %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"), .by = c("tissue", "assay", "comparison_group")) %>%
    filter(feature_ID == desired_feature) %>%
    select(tissue, sex, comparison_group, logFC) %>%
    tidyr::pivot_wider(id_cols = tissue,
                       names_from = c(sex, comparison_group),
                       values_from = logFC) %>%
    column_to_rownames("tissue") %>%
    as.matrix()

  x_range <- if (!is.null(color_range)) color_range else range(all_immuno_da_results, na.rm = TRUE)
  color_breaks <- c(x_range[1], 0, x_range[2])
  legend_breaks <- 0.1 * c(floor(color_breaks[1] / 0.1), 0,
                            ceiling(color_breaks[3] / 0.1))

  immuno_pvalue_annotations = MotrpacRatTraining6mo::combine_da_results(tissues = immuno_tissue_types,
                                                                        assays = "IMMUNO",
                                                                        include_epigen = F) %>%
    filter(feature_ID %in% c("GIP", "GLP1", "GLUCAGON")) %>%
    mutate(adj_p_value = p.adjust(p_value, method = "BH"), .by = c("tissue", "assay", "comparison_group")) %>%
    filter(feature_ID == desired_feature) %>%
    select(tissue, sex, comparison_group, adj_p_value) %>%
    tidyr::pivot_wider(id_cols = tissue,
                       names_from = c(sex, comparison_group),
                       values_from = adj_p_value) %>%
    column_to_rownames("tissue") %>%
    as.data.frame()
  annotation_matrix <- ifelse(immuno_pvalue_annotations < 0.1, "*", "")

  col_meta = data.frame(
    col_name = colnames(all_immuno_da_results),
    sex = sub("_.*", "", colnames(all_immuno_da_results)),
    timepoint = sub(".*_", "", colnames(all_immuno_da_results))
  )

  top_anno = ComplexHeatmap::HeatmapAnnotation(
    Sex = col_meta$sex,
    Timepoint = col_meta$timepoint,
    col = list(
      Sex = c("female" = "#ff6eff", "male" = "#5555ff"),
      Timepoint = c("1w" = "#F7FCB9", "2w" = "#ADDD8E", "4w" = "#238443", "8w" = "#002612")
    ),
    show_legend = FALSE,
    show_annotation_name = show_anno_names,
    border = TRUE,
    gp = gpar(col = "black"),
    gap = 0,
    simple_anno_size = unit(9, "pt"),
    annotation_name_gp = gpar(fontsize = 7)
  )

  p = Heatmap(
    matrix = all_immuno_da_results,
    col = circlize::colorRamp2(breaks = color_breaks, colors = c("#3366ff", "white", "darkred")),
    na_col = "grey85",
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_title = paste0(desired_feature, " peptide"),
    column_title_gp = gpar(fontsize = 9, fontface = "bold"),
    show_column_names = TRUE,
    column_labels = col_meta$timepoint,
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_split = col_meta$sex,
    column_gap = unit(2, "mm"),
    top_annotation = top_anno,
    border = "black",
    rect_gp = gpar(col = "white", lwd = 0.5),
    row_names_gp = gpar(fontsize = 7),
    row_names_side = "left",
    width = ncol(all_immuno_da_results) * unit(24, "pt"),
    height = nrow(all_immuno_da_results) * unit(24, "pt"),
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(
      title = "logFC vs.\nControl",
      at = legend_breaks,
      labels = legend_breaks,
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      border = "black"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (annotation_matrix[i, j] == "*") {
        grid.text("*", x = x, y = y + height * 0.1,
                  gp = gpar(fontsize = 12, col = "black"))
      }
    }
  )

  annotation_legend = Legend(
    labels = "FDR < 0.1",
    title = "Significance",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    pch = "*",
    type = "points"
  )
  figure_return = list()
  figure_return$plot = p
  figure_return$legend = annotation_legend
  return(figure_return)
}

#as a filter, we implement the DA based on the
get_motrpac_expression = function(genes, tissues, motrpac_da) {
  transcript_filt = motrpac_da %>%
    filter(assay == "TRNSCRPT") %>%
    distinct(gene_symbol, tissue)

  trns_feat_to_gene = MotrpacRatTraining6mo::load_feature_annotation(assay = "TRNSCRPT")
  gene_lookup = trns_feat_to_gene %>%
    dplyr::filter(tolower(gene_name) %in% tolower(genes)) %>%
    dplyr::distinct(gene_id, gene_name)

  expr_table = lapply(seq_len(nrow(gene_lookup)), function(i) {
    df = get_expression_stats(gene_lookup$gene_id[i], tissues, "TRNSCRPT")
    if (nrow(df) == 0) return(NULL)
    df$Gene = gene_lookup$gene_name[i]
    df
  }) %>%
    dplyr::bind_rows()  %>%
    dplyr::mutate(detected = paste(toupper(Gene), as.character(Tissue)) %in%
                   paste(toupper(transcript_filt$gene_symbol), transcript_filt$tissue))

  expr_table$Tissue = factor(expr_table$Tissue, levels = names(MotrpacBicQC::tissue_cols))
  expr_table$Gene = factor(expr_table$Gene)

  tidyr::expand_grid(tissue = tissues, gene = toupper(genes)) %>%
    dplyr::left_join(
      expr_table %>%
        dplyr::mutate(tissue = as.character(Tissue), gene = toupper(as.character(Gene))) %>%
        dplyr::select(tissue, gene, detected) %>%
        dplyr::distinct(),
      by = c("tissue", "gene")
    ) %>%
    dplyr::mutate(source = "rat")
}

get_gtex_expression = function(genes, rat_to_gtex, tissues = names(rat_to_gtex)) {
  gtex_gene_ids = lapply(genes, function(g) {
    res = get_gene_search(geneId = g)
    if (nrow(res) == 0) return(NULL)
    res %>%
      dplyr::filter(geneSymbol == g) %>%
      dplyr::select(geneSymbol, gencodeId)
  }) %>%
    dplyr::bind_rows()

  gtex_expr = get_median_gene_expression(
    gencodeIds = gtex_gene_ids$gencodeId,
    tissueSiteDetailIds = unique(rat_to_gtex)
  ) %>%
    dplyr::select(geneSymbol, tissueSiteDetailId, median_tpm = median)

  rat_tissue_map_df = tibble::tibble(
    rat_tissue = names(rat_to_gtex),
    tissueSiteDetailId = unname(rat_to_gtex)
  )

  gtex_with_rat = gtex_expr %>%
    dplyr::left_join(rat_tissue_map_df, by = "tissueSiteDetailId",
                     relationship = "many-to-many") %>%
    dplyr::mutate(detected = median_tpm > 0.1)

  tidyr::expand_grid(tissue = tissues, gene = genes) %>%
    dplyr::left_join(
      gtex_with_rat %>%
        dplyr::mutate(tissue = as.character(rat_tissue), gene = as.character(geneSymbol)) %>%
        dplyr::distinct(tissue, gene, detected),
      by = c("tissue", "gene")
    ) %>%
    dplyr::mutate(source = "human")
}

make_triangle_poly = function(df, tissue_levels, gene_levels) {
  lapply(seq_len(nrow(df)), function(i) {
    x_num = which(tissue_levels == df$tissue[i])
    y_num = which(gene_levels == df$gene[i])
    fill_val = if (is.na(df$detected[i])) "no_data" else if (df$detected[i]) "detected" else "not_detected"
    if (df$source[i] == "rat") {
      px = c(x_num - 0.5, x_num + 0.5, x_num + 0.5)
      py = c(y_num + 0.5, y_num + 0.5, y_num - 0.5)
    } else {
      px = c(x_num - 0.5, x_num - 0.5, x_num + 0.5)
      py = c(y_num + 0.5, y_num - 0.5, y_num - 0.5)
    }
    data.frame(px = px, py = py, fill_val = fill_val,
               group = paste(df$tissue[i], df$gene[i], df$source[i], sep = "_"))
  }) %>% dplyr::bind_rows()
}

summarize_detection = function(rat_data, gtex_data, tissues = NULL) {
  combined = dplyr::inner_join(
    rat_data %>% dplyr::select(tissue, gene, rat_detected = detected),
    gtex_data %>% dplyr::select(tissue, gene, gtex_detected = detected),
    by = c("tissue", "gene")
  )

  fmt = function(x) paste0(toupper(substr(x, 1, 1)), tolower(substring(x, 2)))

  tissue_iter = if (!is.null(tissues)) tissues else unique(combined$tissue)

  summaries = lapply(tissue_iter, function(t) {
    td = combined %>% dplyr::filter(tissue == t)

    rat_genes = td %>%
      dplyr::filter(!is.na(rat_detected) & rat_detected) %>%
      dplyr::pull(gene)

    disagree_td = td %>%
      dplyr::filter(!is.na(rat_detected) & !is.na(gtex_detected) & rat_detected != gtex_detected)

    rat_only = disagree_td %>% dplyr::filter(rat_detected & !gtex_detected) %>% dplyr::pull(gene)
    human_only = disagree_td %>% dplyr::filter(!rat_detected & gtex_detected) %>% dplyr::pull(gene)

    rat_str = if (length(rat_genes) == 0) {
      "no transcript expression detected in rat"
    } else {
      paste0("Transcriptomic expression for ", paste(fmt(rat_genes), collapse = ", "), " was detected in rats")
    }

    disagree_parts = c(
      if (length(rat_only) > 0) paste0(paste(fmt(rat_only), collapse = ", "), " detected in rats only"),
      if (length(human_only) > 0) paste0(paste(fmt(human_only), collapse = ", "), " detected in humans only")
    )
    disagree_str = if (length(disagree_parts) == 0) {
      "GTEx human data agrees"
    } else {
      paste0("GTEx human data disagrees: ", paste(disagree_parts, collapse = "; "))
    }

    paste0(t, ": ", rat_str, ". ", disagree_str, ".")
  })

  paste(summaries, collapse = "\n")
}

compute_detection_summary = function(rat_data, gtex_data, tissues = NULL) {
  combined = dplyr::inner_join(
    rat_data %>% dplyr::select(tissue, gene, rat_detected = detected),
    gtex_data %>% dplyr::select(tissue, gene, gtex_detected = detected),
    by = c("tissue", "gene")
  )

  if (!is.null(tissues)) {
    combined = combined %>% dplyr::filter(tissue %in% tissues) %>%
      dplyr::mutate(tissue = factor(tissue, levels = tissues))
  }

  by_tissue = combined %>%
    dplyr::group_by(tissue) %>%
    dplyr::summarise(
      n_genes = dplyr::n(),
      n_comparable = sum(!is.na(gtex_detected)),
      rat_n_detected = sum(rat_detected, na.rm = TRUE),
      rat_pct_detected = round(100 * rat_n_detected / n_genes, 1),
      human_n_detected = sum(gtex_detected, na.rm = TRUE),
      human_pct_detected = round(100 * human_n_detected / n_comparable, 1),
      n_agree = sum(rat_detected == gtex_detected, na.rm = TRUE),
      pct_agree = round(100 * n_agree / n_comparable, 1),
      .groups = "drop"
    )

  overall = combined %>%
    dplyr::summarise(
      tissue = factor("Overall"),
      n_genes = dplyr::n(),
      n_comparable = sum(!is.na(gtex_detected)),
      rat_n_detected = sum(rat_detected, na.rm = TRUE),
      rat_pct_detected = round(100 * rat_n_detected / n_genes, 1),
      human_n_detected = sum(gtex_detected, na.rm = TRUE),
      human_pct_detected = round(100 * human_n_detected / n_comparable, 1),
      n_agree = sum(rat_detected == gtex_detected, na.rm = TRUE),
      pct_agree = round(100 * n_agree / n_comparable, 1)
    )

  dplyr::bind_rows(overall, by_tissue) %>%
    dplyr::select(tissue, rat_n_detected, n_genes, rat_pct_detected,
                  human_n_detected, n_comparable, human_pct_detected,
                  n_agree, pct_agree)
}

