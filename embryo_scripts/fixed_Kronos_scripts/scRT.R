scRTplot = function(pseudoBulkRT,
                    S_scCN,
                    Coordinates,
                    Plot = 'scRT',
                    rasterized_heatmap = F,
                    heatmap_colors = NULL,
                    sample_colors = NULL) {
  #load required operators
  `%>%` = tidyr::`%>%`

  #set theme
  ggplot2::theme_set(new = ggplot2::theme_bw())


  #check Plot

  if (!Plot %in% c('scRT', 'S_scCN', 'Norm. S_scCN')) {
    stop('Plot can be either scRT, S_scCN or Norm. S_scCN')
  }


  #calculate extremes
  extremes = S_scCN %>%
    dplyr::ungroup() %>%
    dplyr::summarise(CN_bg = round(stats::quantile(CN_bg, c(0.01, 0.99)), 1),
                     CN = round(stats::quantile(CN, c(0.01, 0.99)), 1))
  #filter data
  pseudoBulkRT = pseudoBulkRT %>%
    dplyr::filter(chr == Coordinates$chr,
                  start >= Coordinates$start,
                  end <= Coordinates$end)
  S_scCN = S_scCN %>%
    dplyr::filter(chr == Coordinates$chr,
                  start >= Coordinates$start,
                  end <= Coordinates$end)

  #depending on selection change aesthetics of the plot
  if (Plot == 'scRT') {
    if (is.null(heatmap_colors)) {
      heatmap_colors = c("Replicated" = '#a7001b',
                         "Unreplicated" = '#005095')
    } else if (any(is.na(names(heatmap_colors)))) {
      names(heatmap_colors) = c("Replicated", "Unreplicated")
    }

    S_scCN = S_scCN %>%
      dplyr::rename(Value = Rep) %>%
      dplyr::mutate(Value = ifelse(Value, "Replicated", "Unreplicated"))

    ggEXTRA = ggplot2::ggplot() +
      ggplot2::scale_fill_manual(values = heatmap_colors) +
      ggplot2::labs(fill = 'State', color = 'Sample')

  } else if (Plot == 'S_scCN') {
    if (is.null(heatmap_colors)) {
      heatmap_colors = c('#c56700',
                         '#6d0042')
    }

    S_scCN = S_scCN %>%
      dplyr::rename(Value = CN)

    ggEXTRA = ggplot2::ggplot() +
      ggplot2::scale_fill_gradient(
        low = heatmap_colors[1],
        high = heatmap_colors[2],
        limits = c(extremes$CN[1], extremes$CN[2]),
        oob = scales::squish
      ) +
      ggplot2::labs(fill = 'CNV ', color = 'Sample')

  } else if (Plot == 'Norm. scCN') {
    if (is.null(heatmap_colors)) {
      heatmap_colors = c('#00c5a4',
                         '#e58225')
    }

    S_scCN = S_scCN %>%
      dplyr::rename(Value = CN_bg)

    ggEXTRA = ggplot2::ggplot() +
      ggplot2::scale_fill_gradient(
        low = heatmap_colors[1],
        high = heatmap_colors[2],
        limits = c(extremes$CN_bg[1], extremes$CN_bg[2]),
        oob = scales::squish
      ) +
      ggplot2::labs(fill = expression(over(S[CNV[i]],
                                           bar(G1 / G2[CNV]))),
                    color = 'Sample')

  }

  Chr = Coordinates$chr
  Maxi = max(S_scCN$newIndex)
  n = stringr::str_count(Maxi)
  basenames_colors = unique(pseudoBulkRT$basename)

  #set colors for basename
  if (is.null(sample_colors)) {
    sample_colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(basenames_colors))
    names(sample_colors) = basenames_colors
  } else if (length(sample_colors) < length(basenames_colors) &
             is.null(names(sample_colors))) {
    sample_colors = c(sample_colors,
                      colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(basenames_colors) -
                                                                               length(sample_colors)))

    #name the vector
    names(sample_colors) = basenames_colors

  } else if (length(sample_colors) < length(basenames_colors) &
             !any(is.null(names(sample_colors)))) {
    #if the given colors are less than the number of types and the vector is not named
    default_colors = colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(basenames_colors) -
                                                                              length(sample_colors))

    default_colors = default_colors[!names(default_colors) %in% names(sample_colors)]

    #merge the two
    sample_colors = c(sample_colors, default_colors)
  }

    plot =
    ggEXTRA +
    ggplot2::geom_path(
      data = pseudoBulkRT %>%
        dplyr::mutate(mid = (start + end) / 2,
                      RT = RT * Maxi / 6 + Maxi / 20) %>%
        tidyr::gather(Plot, pos, start, end) %>%
        dplyr::arrange(mid, pos),
      ggplot2::aes(pos, RT, color = basename)
    ) +
    ggplot2::annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = Maxi / 100,
      ymax = Maxi / 20 - Maxi / 100,
      fill = 'white'
    ) +
    ggplot2::annotate(
      "rect",
      xmin = Coordinates$start,
      xmax = Coordinates$end,
      ymin = - Maxi - 1,
      ymax = 0,
      fill = heatmap_colors[2])
    
  if (rasterized_heatmap) {
    plot = plot + ggrastr::geom_tile_rast(data = S_scCN ,
                                          ggplot2::aes(
                                            x = start,
                                            y = -newIndex,
                                            fill = Value
                                          ))


  } else{
    plot = plot +
      ggplot2::geom_rect(
        data = S_scCN ,
        ggplot2::aes(
          xmin = start,
          xmax = end,
          ymin = -newIndex,
          ymax = -newIndex - 1,
          fill = Value
        )
      )
  }

  plot = plot +
    ggplot2::scale_y_continuous(
      breaks = c(Maxi / 6 + Maxi / 20, Maxi / 3 + Maxi / 20, Maxi / 20),
      labels = c('Early - 1', 'Mind - 0.5', 'Late - 0'),
      name = 'RT',
      sec.axis = ggplot2::sec_axis(
        ~ .,
        breaks = c(-seq(1, round(Maxi / 10 ^ (
          n - 1
        )), 1) * 10 ^ (n - 1)) - 0.5,
        labels = as.character(c(seq(
          1, round(Maxi / 10 ^ (n - 1)), 1
        ) * 10 ^ (n - 1))),
        name = 'Single Cell tracks ordered by S-phase progression'
      )
    ) +
    ggplot2::scale_x_continuous(
      labels = function(x)
        paste(x / 10 ^ 6, 'Mb', sep = ' ')
    ) + ggplot2::theme(
      legend.position = 'right',
      legend.direction = "vertical",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.title.y.right = ggplot2::element_text(hjust = 0.6, vjust =
                                                   2),
      axis.title.y.left  = ggplot2::element_text(hjust = 0.92, vjust =
                                                   2)

    ) + ggplot2::xlab(Chr) +
    ggplot2::scale_color_manual(values = sample_colors) +
    ggplot2::facet_grid( ~ group)
  return(plot)

}
