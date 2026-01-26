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
  merged_data$missing_transcript = merged_data$Average.x < 0.2

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
        ggbreak::scale_y_break(brk, scales = "fixed", space = 0.1)
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
                                       immuno_data,
                                       gene_name,
                                       include_first_legend = TRUE,
                                       rna_y_breaks = NULL,
                                       immuno_y_breaks = NULL,
                                       immuno_reverse = TRUE) {
  merged_data = full_join(rna_data, immuno_data, by = c("Tissue", "Sex"))
  merged_data$Tissue = factor(merged_data$Tissue, levels = names(MotrpacBicQC::tissue_cols))
  #tissues run in transcript but not immuno
  merged_data$missing_immuno = is.na(merged_data$Average.y)
  #so here - the luminex assays were run in 17 tissues, not run in 2.
  #in some of the tissues, the feature was not detected -> is na
  merged_data[merged_data$Tissue %in% c("HYPOTH", "VENACV"),]$missing_immuno = "Not run"
  merged_data$missing_transcript = merged_data$Average.x < 0.2

  if(gene_name == "GLP1R") y_adjust = 10
  if(gene_name == "GIPR") y_adjust = 4

  rna_plot = ggplot(merged_data, aes(x = Tissue, y = Average.x, fill = Tissue, group = Sex)) +
    geom_bar(
      stat = "identity",
      width = 0.7,
      position = position_dodge(width = 0.8)
    ) +
    geom_errorbar(
      aes(ymin = Average.x - SD.x, ymax = Average.x + SD.x),
      width = 0.25,
      position = position_dodge(width = 0.8)
    ) +
    geom_text(
      data = merged_data %>% filter(missing_transcript),
      aes(x = Tissue, y = y_adjust, label = "< 0.2 CPM.", group = Sex),
      # position = position_dodge(width = 0.8), #don't have 2 sets of description for this for both sex
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

  if (!is.null(rna_y_breaks))
    # rna_plot = rna_plot + ggbreak::scale_y_break(rna_y_breaks, scales = "fixed", space = 0.1)
    rna_plot = rna_plot + ggforce::facet_zoom(y = Average.x < rna_y_breaks[1], zoom.size = 2)
  # ggforce::facet_zoom(y = Average.x < rna_y_breaks[1], zoom.size = 1)

  #----make sure its putting the metabolite level not the gene name-------
  signalling_molecule_name = sub("R$", "", gene_name)
  #---this just pops off the R at the end to remove the receptor----------

  immuno_plot = ggplot(merged_data, aes(x = Tissue, y = Average.y, fill = Tissue, group = Sex)) +
    geom_bar(stat = "identity", width = 0.7, position = position_dodge(0.8)) +
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
    ylab("Immunoassay\n(Values relative to plate standards)") +
    ggtitle(paste("Immunoassay Expression -", signalling_molecule_name)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
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
                                           interested_genes){
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
      tidyr::complete(gene_symbol = all_genes, group = expected_groups, fill = list(zscore = NaN)) %>%
      ungroup() %>%
      select(gene_symbol, assay, tissue, group, zscore, adj_p_value, sex) %>%
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

    ome_metadata_zscore = make_ome_metadata_matrix("zscore")
    ome_metadata_fdr = make_ome_metadata_matrix("adj_p_value")

    combined_ome_metadata[[tissue_value]] <- ome_metadata_zscore
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
    df = annotation_df %>% select(Tissue, Sex, Timepoint),
    border = TRUE,
    gp = gpar(col = "black"),
    gap = 0,
    which = "column",
    height = unit(6 * 2, "pt")*scale,
    col = list(
      Tissue = MotrpacBicQC::tissue_cols,
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
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    matrix = final_ome_metadata,
    col = circlize::colorRamp2(
      breaks = c(-1.6, 0, 1.8),
      colors = c("#3366ff", "white", "darkred")
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    show_column_names = FALSE,
    clustering_distance_rows = "pearson",
    top_annotation = combined_annotation,
    border = "black",
    row_names_gp = gpar(fontsize = 5 * scale),
    height = nrow(final_ome_metadata) * unit(5.5, "pt") * scale,
    width = ncol(final_ome_metadata) * unit(5.5, "pt") * scale,
    column_split = rep(1:(2*length(select_tissue_types)), each = 4),
    column_title = NULL,
    heatmap_legend_param = list(
      title = "Z-Score \nper Group/Sex",
      at = c(-3, -1:1, 3),
      title_gp = gpar(fontsize = 7 * scale, fontface = "bold"),
      labels_gp = gpar(fontsize = 6 * scale),
      legend_height = 5 * scale * unit(8, "pt"),
      border = "black"
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      #Heatmap grid
      grid.rect(x = x, y = y, width, height,
                gp = gpar(col = "#555555"))
      #update to use FDR instead of easy Zscore filt
      if(!is.na(filtered_fdr[i, j]) && abs(filtered_fdr[i, j]) < 0.1) {
        gb = textGrob("*")
        gb_w = convertWidth(grobWidth(gb), "mm")
        gb_h = convertHeight(grobHeight(gb), "mm")
        grid.text("*", x, y - gb_h*0.5 + gb_w*0.4)
      }

    }
  )

}

make_immuno_heatmap = function(desired_feature){
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

  column_labels <- colnames(all_immuno_da_results)
  x_range <- range(all_immuno_da_results, na.rm = TRUE) #range of logFC
  color_breaks <- c(x_range[1], 0, x_range[2])
  # Extend range of values out to the nearest 0.1
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
  #use a 0.1 FDR cutoff.
  annotation_matrix <- ifelse(immuno_pvalue_annotations < 0.1, "*", "")

  p = Heatmap(
    matrix = all_immuno_da_results,
    col = circlize::colorRamp2(
      breaks = color_breaks,
      colors = c("#5555ff", "white", "darkred")
    ),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    column_title = desired_feature,
    column_labels = column_labels,
    heatmap_legend_param = list(
      title = paste("logFC", desired_feature),
      at = legend_breaks,
      labels = legend_breaks
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      # Add stars for cells with significant p-values
      if (annotation_matrix[i, j] == "*") {
        grid.text(
          "*", x = x, y = y,
          gp = gpar(fontsize = 18, col = "black")
        )
      }
    }
  )

  annotation_legend = Legend(
    labels = c("p < 0.05"),
    title = "Significance",
    legend_gp = gpar(fontsize = 10),
    pch = "*",  # Use the star symbol
    type = "points"
  )
  figure_return = list()
  figure_return$plot = p
  figure_return$legend = annotation_legend
  return(figure_return)
}

