#' read pca result from a prefix
#' 
#' assuming all files is available with the same prefix
#' NOTES: this is also the constructor of the pca result object
#' 
#' @param prefix prefix of pca result
#' @param loading if TRUE, read loading of principle components, default FALSE
#' @return list of df
#' @export
read_pca_result <- function(prefix, loading = FALSE) {
  # read eigenval of principle components
  eigenvals = read_lines(paste0(prefix, ".eigenval")) %>% as.numeric()
  
  # read the sum of covariance
  sum = read_lines(paste0(prefix, ".sum")) |> as.numeric()
  
  # compute the proportion of variance explained by each PC
  pves = eigenvals / sum
  
  # read the loading of principle components
  if (loading) {
    # possible files for loading
    files = paste0(prefix, c(".eigenvec.var.gz", ".eigenvec.var.zst", ".eigenvec.var"))
    
    # check if file exist
    files = files[file.exists(files)]
    if (length(files) > 0) {
      # read the file, priotitize the compressed file
      file = files[1]
      # get number of fields of the `eigenvec.var` file from the 1st row of file
      nfields = read_lines(file, n_max = 1) %>% str_split("\\s+") %>% .[[1]] %>% length()
      df.load = fread(file, header = FALSE,
                      col.names = c("chr", "id", "minor", "major", paste0("PC", 1:(nfields - 4)))
      ) |> as_tibble()
    } else {
      df.load = NULL
    }
  } else {
    df.load = NULL
  }
  
  # read the pca score of principle components
  file = paste0(prefix, ".eigenvec")
  # check if the eigenvec file have header
  line1 = readLines(paste0(prefix, ".eigenvec"), n = 1)
  have_header = read_lines(file, n_max = 1) %>% str_detect("^#|^FID")
  delim = str_extract(line1, "\\s")
  if (have_header) {
    # read the pca score of principle components # plink1 output pca are space delimited
    df.score = read_delim(file, delim = delim, show_col_types = FALSE) |>
      # remove `#` from colnames 
      rename_with(~str_remove(.x, "^#"))
  } else {
    # get number of fields of the `eigenvec` file from the 1st row of file
    nfields = read_lines(file, n_max = 1) %>% str_split("\\s+") %>% .[[1]] %>% length()
    
    # compare to the number of vector in eigenvals
    if (nfields - length(eigenvals) > 1) {
      .colnames = c("fid", "gaid", paste0("PC", 1:length(eigenvals)))
    } else {
      .colnames = c("gaid", paste0("PC", 1:length(eigenvals)))
    }
    
    # read the pca score of principle components
    df.score = fread(file, header = FALSE,
                     col.names = .colnames
    ) |> as_tibble()
  }
  
  # return list of dataframes with comment status
  ls.out = list("pves" = pves,
                "load" = df.load,
                "score" = df.score)
  comment(ls.out) = "pca"
  return(ls.out)
}

#' Evaluate highlight condition
#' 
#' wrap this feature to reused in multiple functions. apparently, it's useful for making parsed valued avail in sub-env
#' @param df data frame to evaluate the highlight condition
#' @param highlight a string equivalence of the R expression that can be used to filter the data frame for highlight
#' @return df for plotting
#' @export
eval_hli_rule <- function(df, highlight) {
  # highlight = params$highlight # ls.pca = results
  # parse highlight condition
  hli.value <- str_extract(highlight, '(?<=\\").*(?=\\")')
  hli.key <- str_replace(highlight, ' .*', '')
  
  # check if the hli.value contain any special character
  if (str_detect(hli.value, "[[:alnum:]]")) {
    # if pattern contain special character # turn it into a working R expression for matching regular expression
    # ERROR: cannot change the value of a locked binding for 'highlight'
    highlight <- stringr::str_glue("grepl('{hli.value}', {hli.key})")
  }
  
  # extract df of pca score from the list
  df = df %>%
    # mutate highlight var info based on highlight condition
    mutate(hli = case_when(eval(parse(text = highlight)) ~ hli.value,
                           TRUE ~ "Others"))
    # # turn NA into "NA" # for label function to work properly
    # mutate(region = forcats::fct_explicit_na(region, "NA"))
  
  # return df for plotting
  return(df)
}

#' Plot PCA score in 2d space, annotate number of sample in legend and highlight sample of interest
#' 
#' add highlight points for another data frame based on sample_provider == "PRISM"
#' add another hli column to save the highlight info so that the original info can be used for labeling
#' refactor to work with mds also
#' 
#' @param results list of object from the function `read_pca_result`, with score updated with metadata
#' @param vois Variables of interest. Names of 2 vars to plot on x and y axis
#' @param highlight a string equivalence of the R expression that can be used to filter the data frame for highlight
#' @param hli.label logical, add geom_label_repel to label the highlighted points
#' @param color_by_var variable to color by, default to region
#' @param ... other arguments to pass to custom plot
#' @return a ggplot object
#' @examples
#' # example code
#' plot_pca_score_with_highlight(ls.pca, vois = c("PC1", "PC2"), highlight = 'town == "Andaman"', hli.label = TRUE)
#' @export
plot_pca_score_with_highlight <- function(results, color_by_var = "region",
  vois = c("PC1", "PC2"), highlight = 'town == "Andaman"', hli.label = FALSE,
  ...)
{
  # highlight = params$highlight # results = results
  # extract df of pca score from the list
  df = results$score
  
  # extract pves from the list
  pves = results$pves |>
    set_names(paste0("PC", 1:length(results$pves)))
  
  # choose appropriate corner for legend
  theme.corner = choose_corner(df, vois)
  
  # ---- qual/seq appearance ----
  # label function need to be inside the plot function
  label_fun <- function(values) {
    sapply(values, function(value) {
      stringr::str_c(value, " (", sum(df[[color_by_var]] == value), ")")
    })
  }
  
  # df of colors # rank by priority
  df.color = df |> count(.data[[color_by_var]]) |>
    arrange(desc(n))
  
  if (color_by_var == "region") {
    # use the same .colors vector
    # create scale color layer
    list_layer_ggplot = scale_color_manual(
      name = paste0(make_clean_names(color_by_var, "title"), " (all=", nrow(df), ")"),
      values = .colors, breaks = df.color[[color_by_var]],
      labels = ~ label_fun(.x)
    )
  } else {
    # qualitative/sequential color
    # qualitative is limited by number of color and do not note the difference in the number of sample
    # SPECTRAL palette function 
    brewer_pal_spectral <- brewer_pal(palette = "Spectral", direction = 1)
    
    # maximum 11 color
    if (nrow(df.color) <= 11) {
      color_brewed <- brewer_pal_spectral(nrow(df.color))
    } else {
      # otherwise, brew more colors
      color_brewed <- brewer_pal_spectral(11)
      color_brewed <- pal_gradient_n(color_brewed)(seq(0, 1, length.out = nrow(df.color)))
    }
    
    # apply the brewed color
    df.color = df.color |>
      mutate(color = color_brewed)
    
    # use df.color to transform df so that the pop grouped_by color appear the least will appear on top.
    df = df |>
      mutate(!!sym(color_by_var) := factor(.data[[color_by_var]], levels = df.color[[color_by_var]])) |>
      arrange(!!sym(color_by_var))
    
    # create scale color layer
    list_layer_ggplot = scale_color_manual(
        name = paste0(make_clean_names(color_by_var, "title"), " (all=", nrow(df), ")"),
        values = df.color$color, breaks = df.color[[color_by_var]], 
        labels = ~ label_fun(.x)
      )
  }
  
  # add theme.corner and lab to aesthetic layer
  list_layer_ggplot = c(
    list_layer_ggplot,
    list(
      theme.corner,
      labs(x = paste0(vois[1], " (", percent(pves[vois[1]], 0.01), ")"),
          y = paste0(vois[2], " (", percent(pves[vois[2]], 0.01), ")"))
    )
  )
  
  # check if highlight argument is supplied
  if (missing(highlight)) {
    # If the argument is missing, DONT highlight with SHAPE # missing still works even with declared default
    fig = ggplot(
      data = df,
      mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]],
                    color = .data[[color_by_var]])
    ) +
      geom_point(size = 3, alpha = 0.6) +
      list_layer_ggplot
  } else {
    # if highlight argument is supplied, create the data frame for the highlight layer
    df.hli = df |>
      filter(eval(parse(text = highlight)))
    
    # and df for the label of highlighted sample
    df.label.hli = df.hli %>%
      group_by(ethnicity, country, region) %>%
      # count number of sample for each group and calculate the median for each vois
      summarize(count = n(), x.median = median(!!sym(vois[1])), y.median = median(!!sym(vois[2])),
                .groups = "drop") %>%
      # select top N label for each region based on sample count
      group_by(region) %>%
      arrange(region, desc(count)) %>% filter(count > 1) %>%
      mutate(label = paste0(ethnicity, "(", count, ")"))
    
    # plot scatter plot
    fig = ggplot(
      data = df,
      mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]],
                    color = .data[[color_by_var]], shape = hli)
    ) +
      geom_point(size = 3, alpha = 0.6) +
      # custom parameter for the shape of the plot and LEGEND
      scale_shape_manual(name = make_clean_names(str_remove(highlight, " .*"), "big_camel"),
                         breaks = unique(df$hli),
                         values = c(16, 18)) +
      list_layer_ggplot
    
    # add layer of only the shape with NO legend
    fig = fig +
      ggnewscale::new_scale_color() +
      geom_point(data = df.hli,
        mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]], fill = region),
        shape = 23, size = 4, color = "black"
      ) +
      scale_fill_manual(guide = "none", values = .colors)
  }
  
  # annotate with N ethnicity label for each region
  N = 4
  df.label = df %>%
    group_by(ethnicity, town, country, region) %>%
    # count number of sample for each group and calculate the median for each vois
    summarize(count = n(), x.median = median(!!sym(vois[1])), y.median = median(!!sym(vois[2])),
              .groups = "drop") %>%
    # select top N label for each region based on sample count
    group_by(region) |>
    # replace top_n with slice_max
    slice_max(order_by = count, n = N) %>%
    arrange(region, desc(count))
  
  # annotate with ggrepel
  fig = fig +
    ggrepel::geom_text_repel(data = df.label,
      aes(x = x.median, y = y.median, label = paste0(ethnicity, "\n(", country, ")")), inherit.aes = FALSE,
      force = 1, box.padding = 0.5,
      segment.size = 0.2, segment.color = "grey50",
      max.overlaps = 40, size = 2)
  
  # check whether to label highlighted data with "hli.label"
  if (hli.label) {
    fig = fig +
      ggrepel::geom_label_repel(data = df.label.hli,
        aes(x = x.median, y = y.median, label = label), inherit.aes = FALSE,
        force = 1, box.padding = 0.5,
        # segment params
        segment.size = 0.5, segment.color = "black",
        max.overlaps = 40, size = 2)
  }
  return(fig)
}

#' validate dimension reduction result
#' 
#' @keyword Internal
validate_dm_result <- function(result) {
  if (comment(result) == "pca") {
    component_prefix = "PC"
    # assign("component_prefix", "PC", envir = parent.frame())
  } else if (comment(result) == "mds") {
    component_prefix = "C"
    # assign("component_prefix", "C", envir = parent.frame())
  } else {
    stop("Invalid result object")
  }
  return(component_prefix)
}

#' Plot dimension reduction score in 2d space
#' 
#' annotate number of sample in legend and highlight sample of interest
#' refactor to work with mds.
#' 
#' @param results list of object from the function `read_pca_result`, with score updated with metadata
#' @param vois vector of 2 variable name to plot on x and y axis
#' @param highlight a string equivalence of the R expression that can be used to filter the data frame for highlight
#' @param hli.label logical, add geom_label_repel to label the highlighted points
#' @param color_by_var variable to color by, default to region
#' @param ... other arguments to pass to custom plot
#' @return a ggplot object
#' @examples
#' # example code
#' plot_pca_score_with_highlight(ls.pca, vois = c("PC1", "PC2"), highlight = 'town == "Andaman"', hli.label = TRUE)
#' @export
plot_dimension_reduction <- function(results, color_by_var = "region",
  vois, highlight = 'town == "Andaman"', hli.label = FALSE,
  ...)
{
  # validate results to get component_prefix 
  component_prefix = validate_dm_result(results)
  # highlight = params$highlight # results = results
  # extract df of pca score from the list
  df = results$score
  
  # extract pves from the list
  pves = results$pves |>
    set_names(paste0(component_prefix, 1:length(results$pves)))
  
  # choose appropriate corner for legend
  theme.corner = choose_corner(df, vois)
  
  # ---- qual/seq appearance ----
  # label function need to be inside the plot function
  label_fun <- function(values) {
    sapply(values, function(value) {
      stringr::str_c(value, " (", sum(df[[color_by_var]] == value), ")")
    })
  }
  
  # df of colors # rank by priority
  df.color = df |> count(.data[[color_by_var]]) |>
    arrange(desc(n))
  
  if (color_by_var == "region") {
    # use the same .colors vector
    # create scale color layer
    list_layer_ggplot = scale_color_manual(
      name = paste0(make_clean_names(color_by_var, "title"), " (all=", nrow(df), ")"),
      values = .colors, breaks = df.color[[color_by_var]],
      labels = ~ label_fun(.x)
    )
  } else {
    # qualitative/sequential color
    # qualitative is limited by number of color and do not note the difference in the number of sample
    # SPECTRAL palette function 
    brewer_pal_spectral <- brewer_pal(palette = "Spectral", direction = 1)
    
    # maximum 11 color
    if (nrow(df.color) <= 11) {
      color_brewed <- brewer_pal_spectral(nrow(df.color))
    } else {
      # otherwise, brew more colors
      color_brewed <- brewer_pal_spectral(11)
      color_brewed <- pal_gradient_n(color_brewed)(seq(0, 1, length.out = nrow(df.color)))
    }
    
    # apply the brewed color
    df.color = df.color |>
      mutate(color = color_brewed)
    
    # use df.color to transform df so that the pop grouped_by color appear the least will appear on top.
    df = df |>
      mutate(!!sym(color_by_var) := factor(.data[[color_by_var]], levels = df.color[[color_by_var]])) |>
      arrange(!!sym(color_by_var))
    
    # create scale color layer
    list_layer_ggplot = scale_color_manual(
      name = paste0(make_clean_names(color_by_var, "title"), " (all=", nrow(df), ")"),
      values = df.color$color, breaks = df.color[[color_by_var]], 
      labels = ~ label_fun(.x)
    )
  }
  
  # add theme.corner and lab to aesthetic layer
  list_layer_ggplot = c(
    list_layer_ggplot,
    list(
      theme.corner,
      labs(x = paste0(vois[1], " (", percent(pves[vois[1]], 0.01), ")"),
           y = paste0(vois[2], " (", percent(pves[vois[2]], 0.01), ")"))
    )
  )
  
  # check if vois appear in df
  if (!all(vois %in% colnames(df))) {
    stop("vois must be column names in the data frame")
  }
  
  # check if highlight argument is supplied
  if (missing(highlight)) {
    # If the argument is missing, DONT highlight with SHAPE # missing still works even with declared default
    fig = ggplot(
      data = df,
      mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]],
                    color = .data[[color_by_var]])
    ) +
      geom_point(size = 3, alpha = 0.6) +
      list_layer_ggplot
  } else {
    # if highlight argument is supplied, create the data frame for the highlight layer
    df.hli = df |>
      filter(eval(parse(text = highlight)))
    
    # and df for the label of highlighted sample
    df.label.hli = df.hli %>%
      group_by(ethnicity, country, region) %>%
      # count number of sample for each group and calculate the median for each vois
      summarize(count = n(), x.median = median(!!sym(vois[1])), y.median = median(!!sym(vois[2])),
                .groups = "drop") %>%
      # select top N label for each region based on sample count
      group_by(region) %>%
      arrange(region, desc(count)) %>% filter(count > 1) %>%
      mutate(label = paste0(ethnicity, "(", count, ")"))
    
    # plot scatter plot
    fig = ggplot(
      data = df,
      mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]],
                    color = .data[[color_by_var]], shape = hli)
    ) +
      geom_point(size = 3, alpha = 0.6) +
      # custom parameter for the shape of the plot and LEGEND
      scale_shape_manual(name = make_clean_names(str_remove(highlight, " .*"), "big_camel"),
                         breaks = unique(df$hli),
                         values = c(16, 18)) +
      list_layer_ggplot
    
    # add layer of only the shape with NO legend
    fig = fig +
      ggnewscale::new_scale_color() +
      geom_point(data = df.hli,
                 mapping = aes(x = .data[[vois[1]]], y = .data[[vois[2]]], fill = region),
                 shape = 23, size = 4, color = "black"
      ) +
      scale_fill_manual(guide = "none", values = .colors)
  }
  
  # annotate with N ethnicity label for each region
  N = 4
  df.label = df %>%
    group_by(ethnicity, town, country, region) %>%
    # count number of sample for each group and calculate the median for each vois
    summarize(count = n(), x.median = median(!!sym(vois[1])), y.median = median(!!sym(vois[2])),
              .groups = "drop") %>%
    # select top N label for each region based on sample count
    group_by(region) |>
    # replace top_n with slice_max
    slice_max(order_by = count, n = N) %>%
    arrange(region, desc(count))
  
  # annotate with ggrepel
  fig = fig +
    ggrepel::geom_text_repel(data = df.label,
                             aes(x = x.median, y = y.median, label = paste0(ethnicity, "\n(", country, ")")), inherit.aes = FALSE,
                             force = 1, box.padding = 0.5,
                             segment.size = 0.2, segment.color = "grey50",
                             max.overlaps = 40, size = 2)
  
  # check whether to label highlighted data with "hli.label"
  if (hli.label) {
    fig = fig +
      ggrepel::geom_label_repel(data = df.label.hli,
                                aes(x = x.median, y = y.median, label = label), inherit.aes = FALSE,
                                force = 1, box.padding = 0.5,
                                # segment params
                                segment.size = 0.5, segment.color = "black",
                                max.overlaps = 40, size = 2)
  }
  return(fig)
}

#' Choose appropriate corner for legend
#' 
#' @param df dataframe for plotting
#' @param vois 2 string of variable of interest that represennt the x and y axis in that order
#' @return a theme function that can be added to ggplot object to customize the look
#' @export
choose_corner <- function(df, vois) {
  # prepare location for legend
  df.smm = df %>% select(any_of(vois)) %>%
    rstatix::get_summary_stats() %>%
    mutate(quart = (max - min) / 4)
  
  # set the 4 corner point on given coordinate
  ls.corner = list("lb" = c(df.smm$min[1] + df.smm$quart[1], df.smm$min[2] + df.smm$quart[2]),
                   "rb" = c(df.smm$max[1] - df.smm$quart[1], df.smm$min[2] + df.smm$quart[2]),
                   "lt" = c(df.smm$min[1] + df.smm$quart[1], df.smm$max[2] - df.smm$quart[2]),
                   "rt" = c(df.smm$max[1] - df.smm$quart[1], df.smm$max[2] - df.smm$quart[2]))
  
  # add names of element to the comment of element within the list
  for (i in 1:length(ls.corner)) {
    comment(ls.corner[[i]]) <- names(ls.corner)[i]
  }
  # compute number of point falling within corner of plot
  ls.count = lapply(ls.corner, function(X) {
    # X is the coordinate of the corner points c(x, y)
    # check the condition of corner within the comment of X # compute number of point falling within corner of plot depend on the condition of corner
    if (comment(X) == "lb") {
      n.point = df %>%
        filter(.data[[vois[1]]] < X[1] & .data[[vois[2]]] < X[2]) %>% nrow()
    } else if (comment(X) == "rb") {
      n.point = df %>%
        filter(.data[[vois[1]]] > X[1] & .data[[vois[2]]] < X[2]) %>% nrow()
    } else if (comment(X) == "lt") {
      n.point = df %>%
        filter(.data[[vois[1]]] < X[1] & .data[[vois[2]]] > X[2]) %>% nrow()
    } else if (comment(X) == "rt") {
      n.point = df %>%
        filter(.data[[vois[1]]] > X[1] & .data[[vois[2]]] > X[2]) %>% nrow()
    }
    return(n.point)
  })
  
  # choose the corner based on the number of point falling within the corner
  # if there is no corner of 0 points, then choose the corner with the least number of points
  if (all(unlist(ls.count) > 0)) {
    print("noempty corner")
    corner.select = which.min(unlist(ls.count)) |> names()
  } else {
    corner.select = which.min(unlist(ls.count)) |> names()
  }
  
  # turn the select corner into corresponding theme function to add to ggplot
  if (corner.select == "lb") {
    theme.corner = theme(legend.position.inside = c(0,0),
                         legend.justification = c(0,0))
  } else if (corner.select == "rb") {
    theme.corner = theme(legend.position.inside = c(1,0),
                         legend.justification = c(1,0))
  } else if (corner.select == "lt") {
    theme.corner = theme(legend.position.inside = c(0,1),
                         legend.justification = c(0,1))
  } else if (corner.select == "rt") {
    theme.corner = theme(legend.position.inside = c(1,1),
                         legend.justification = c(1,1))
  }
  return(theme.corner)
}

#' Read list of mds result into a list
#' 
#' Assuming files are avail with the same prefix
#' @param prefix prefix to the mds result
#' @param nc number of mds component to plot
#' @return a list of data frame
#' @export
read_mds_result <- function(prefix, nc = 3) {
  # prefix = "output/ga4k/ga4k_sel.autosomes.common_rare_phased.maf0.01.indep_500kb_0.2.pruned.ibd.mds"
  # read eigvals file
  file = paste0(prefix, ".mds.eigvals")
  eigvals = readLines(file) |> as.numeric() # length(eigvals) # 3948
  
  # read mds result # fread because the file have leading space 
  file = paste0(prefix, ".mds")
  df = data.table::fread(file) |> tibble::as_tibble() |>
    # select only col to plot # +3 for fid,iid,sol
    dplyr::select(1:{nc+3})
  
  # compute pves if all eigvals are available # or at least n-1 eigvals are available
  if (length(eigvals) >= nrow(df) - 1) {
    pves = eigvals / sum(eigvals)
  } else {
    pves = rep(NA, nrow(df))
  }
  
  # return
  ls.df = list(score = df, eigvals = eigvals, pves = pves)
  comment(ls.df) = "mds"
  return(ls.df)
}

#' Scale score of mds/pca based on pves
#' 
#' @param results list of mds/pca result that can be check in comment(results)
#' @return same format with results but with score scaled
#' @export
scale_score <- function(results) {
  # check comment of results to make sure it is mds/pca result
  if (! comment(results) %in% c("mds", "pca")) {
    stop("results must be mds/pca result")
  } else {
    col_prefix = if (comment(results) == "mds") "C" else "PC"
  }
  
  # scale value of column start with "C" to mean = 0 and sd = 1
  data = results$score |>
    mutate(across(starts_with(col_prefix), ~{scale(.x) |> as.vector()}))
  # rstatix::get_summary_stats(data, starts_with("C"))
  
  # multiply C1:C5 by pves[1:5]
  indices = grep(paste0("^", col_prefix), colnames(data))
  components = grep(paste0("^", col_prefix), colnames(data), value = TRUE) |> parse_number()
  
  # scale
  scales = results$pves[components] * 100
  data[,indices] = data[,indices] * matrix(rep(scales, each = nrow(data)), nrow = nrow(data), ncol = length(indices))
  results$score <- data
  #
  return(results)
}

#' plot mds and pca result with plotly
#' 
#' add feature to rotate the scatter3d plot; add arg to highlight some sample or NOT; 
#' add color_by_var to color the scatter3d plot by a variable; work with region but NOT ethnicity
#' BECAUSE: .colors only contain enough color for region. add a param to pass color vector
#' legend names associated with each trace rather than with the legend itself. So redefining names has to occur on each trace / the data itself
#' 
#' @param results list of mds/pca result that can be check in comment(results)
#' @param nc number of mds component to plot, assuming we want 1:nc component
#' @param highlight a string that can be evaluated inside df.plot to get a logical vector "sample_provider == 'PRISM'"
#' @param color_by_var a string that can be evaluated inside df.plot to get a factor vector "region"
#' @param .colors a vector of color to use for color_by_var, declared this way to differ from colors args of plotly
#' @param ... other argument to pass to plotly
#' @return plotly object that can be seen in interactive() or html
#' @examples
#' @export
plotly_3d <- function(results, nc = 3, highlight,
  color_by_var = "region", .colors)
{
  # df to plot
  df.plot = results$score
  
  # check if results is the output from read_mds_result or read_pca_result
  if (comment(results) == "mds") {
    # components to visualize # C for mds
    vois = paste0("C", 1:nc)
  } else if (comment(results) == "pca") {
    # components to visualize # PC for pca
    vois = paste0("PC", 1:nc)
  } else {
    stop("results should be the output from read_mds/pca_result")
  }
  
  # label to appear on hover info # IID for mds
  .labels = c("IID", "ethnicity", "country", "town", "sample_provider")
  
  # check if results$score contain all the .labels
  if (!all(.labels %in% colnames(df.plot))) {
    # if not, then .labels column is all vars that do not contain id or have a number at the end
    .labels = colnames(df.plot) %>%
      grep("id|\\d$", ., invert = TRUE, ignore.case = TRUE, value = TRUE)
    # count max character of each column
    nchar_label = apply(df.plot[,.labels], 2, function(x) max(nchar(x)))
    short_vars = names(nchar_label)[nchar_label < 10] |> rlang::syms()
    long_vars = names(nchar_label)[nchar_label >= 10] |> rlang::syms()
    
    # add label column to df.plot according to condition
    df.plot = df.plot |>
      mutate(label = paste0(
        str_c(!!!short_vars, sep = ";"), "\n",
        str_c(!!!long_vars, sep = ";")
      ))
    
    # compile list of arg for rlang::exec(add_trace, !!!ls.arg)
    ls.arg.add_trace = list(
      type = "scatter3d", mode = "markers",
      # data = df.plot, inherit = FALSE,
      # declare x, y, z axis
      x = ~.data[[vois[1]]], y = ~.data[[vois[2]]], z = ~.data[[vois[3]]],
      # switch to plot based on trace column
      color = ~trace,
      # add marker size independently using I()
      size = I(100),
      # add custom text to hover info
      text = ~label,
      hoverinfo = "text"
    )
  } else {
    # compile list of arg for rlang::exec(add_trace, !!!ls.arg)
    ls.arg.add_trace = list(
      type = "scatter3d", mode = "markers",
      # data = df.plot, inherit = FALSE,
      # declare x, y, z axis
      x = ~.data[[vois[1]]], y = ~.data[[vois[2]]], z = ~.data[[vois[3]]],
      # switch to plot based on trace column
      color = ~trace,
      # add marker size independently using I()
      size = I(100),
      # add custom text to hover info
      text = ~ paste0(
        .data[[.labels[1]]], ";",
        .data[[.labels[2]]], "\n(",
        .data[[.labels[3]]], ";",
        .data[[.labels[4]]], ";",
        .data[[.labels[5]]], ")"
      ),
      hoverinfo = "text"
    )
  }
  
  # generate dictionary of color_by_var and count of each factor first
  dict = df.plot |> count(!!sym(color_by_var)) |>
    arrange(desc(n)) |>
    # mutate(trace = str_glue("{ethnicity}({n})")) |>
    mutate(trace = paste0(!!sym(color_by_var), "(", n, ")")) |>
    select(-n)
  
  # create a trace column so that we can have customized legend key
  df.plot = df.plot |>
    left_join(dict, by = color_by_var) |>
    mutate(trace = factor(trace, levels = dict$trace)) |>
    arrange(trace)
  
  # qualitative/sequential color
  # qualitative is limited by number of color and do not note the difference in the number of sample
  # SPECTRAL palette function 
  brewer_pal_spectral <- brewer_pal(palette = "Spectral", direction = 1)
  
  # maximum 11 color
  if (nrow(dict) <= 11) {
    color_brewed <- brewer_pal_spectral(nrow(dict))
  } else {
    # otherwise, brew more colors
    color_brewed <- brewer_pal_spectral(11)
    color_brewed <- pal_gradient_n(color_brewed)(seq(0, 1, length.out = nrow(dict)))
  }
  
  # apply the brewed color to dictionary
  dict = dict |>
    mutate(color = color_brewed)
  
  # check if .colors is provided
  if (!missing(.colors)) {
    # if PROVIDED, .colors is read from file and correspond to region
    # this is just old code that assign
    ls.arg.add_trace$colors = .colors
  } else {
    # apply colors from dictionary
    ls.arg.add_trace$colors = dict$color
  }
  
  # compile list of arg for rlang::exec("add_trace", !!!ls.arg)
  ls.arg.layout = list(
    legend = list(
      title = list(text = str_glue("<b>{str_to_title(color_by_var)}</b>"))
    ),
    scene = list(
      # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
      camera = list(eye = list(x = 1.25, y = 1.25, z = 0.2),
                    center = list(x = 0, y = 0, z = 0)
      ),
      # add axis title
      xaxis = list(title = paste0(vois[1], " (", percent(results$pves[parse_number(vois[1])], 0.01), ")")),
      yaxis = list(title = paste0(vois[2], " (", percent(results$pves[parse_number(vois[2])], 0.01), ")")),
      zaxis = list(title = paste0(vois[3], " (", percent(results$pves[parse_number(vois[3])], 0.01), ")"))
    )
  )
  
  # check if highlight arg is provided
  if (!missing(highlight)) {
    # parse highlight # highlight = params$highlight
    # add condition to work with both double quoted and single quoted hli.value
    if (str_detect(highlight, "\'")) {
      hli.value = str_extract(highlight, "(?<=\\\').*(?=\\\')")
    } else {
      hli.value = str_extract(highlight, '(?<=\\").*(?=\\")')
    }
    hli.key = str_replace(highlight, ' .*', '')
    
    
    # check if the hli.value contain any special character
    if (str_detect(hli.value, "[[:alnum:]]", negate = TRUE)) {
      # if pattern contain special character # turn it into a working R expression for matching regular expression
      highlight = stringr::str_glue("grepl('{hli.value}', {hli.key})")
    }
    
    # update df.plot with a hli column
    df.plot = df.plot |>
      mutate(hli = case_when(
        eval(parse(text = highlight)) ~ hli.value,
        TRUE ~ "Others"
      ))

    # declare symbol for highlight
    .symbols = c("circle", "x")
    names(.symbols) = c("Others", hli.value)
    
    # add symbol to ls.arg.add_trace
    ls.arg.add_trace <- c(
      ls.arg.add_trace, 
      list(
        symbol = ~hli, symbols = .symbols
      )
    )
    
    # append hli.key to text of legend title in ls.arg.layout
    ls.arg.layout$legend$title$text = str_c(ls.arg.layout$legend$title$text, "/", str_to_title(hli.key))
  }
  
  # declare conflict prefer # this may be the reason why multiple child document failed while 1 success
  # DO NOT USE: conflicted::conflicts_prefer(plotly::layout, .quiet = TRUE)
  plotly_layout = get("layout", envir = asNamespace("plotly"))

  # pipe to exec with magrittr::%>% to add_trace
  fig = plotly::plot_ly(data = df.plot) %>%
    # the tilde in plot_ly args can: record R expression without evaluating it # so that ~col can be search in data env
    {rlang::exec("add_trace", p = ., !!!ls.arg.add_trace)} %>%
    {rlang::exec("plotly_layout", p = ., !!!ls.arg.layout)}
  # plot_ly() |>
  #   add_trace(
  #     type = "scatter3d", mode = "markers", data = df.plot, inherit = FALSE,
  #     # declare x, y, z axis
  #     x = ~.data[[vois[1]]], y = ~.data[[vois[2]]], z = ~.data[[vois[3]]],
  #     # declare color, symbol
  #     color = ~region, colors = .colors, symbol = ~hli, symbols = .symbols,
  #     # add marker size independently using I()
  #     size = I(100),
  #     # add custom text to hover info
  #     text = ~ paste0(
  #       .data[[.labels[1]]], ";",
  #       .data[[.labels[2]]], "\n(",
  #       .data[[.labels[3]]], ";",
  #       .data[[.labels[4]]], ";",
  #       .data[[.labels[5]]], ")"
  #     ),
  #     hoverinfo = "text"
  #   ) |>
  #   layout(
  #     legend = list(
  #       title = list(text = str_glue("Region/{hli.key}"))
  #     ),
  #     scene = list(
  #       # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
  #       camera = list(eye = list(x = 1.25, y = 1.25, z = 0.2),
  #                     center = list(x = 0, y = 0, z = 0)
  #       ),
  #       # add axis title
  #       xaxis = list(title = paste0(vois[1], " (", percent(results$pves[parse_number(vois[1])], 0.01), ")")),
  #       yaxis = list(title = paste0(vois[2], " (", percent(results$pves[parse_number(vois[2])], 0.01), ")")),
  #       zaxis = list(title = paste0(vois[3], " (", percent(results$pves[parse_number(vois[3])], 0.01), ")"))
  #     )
  #   )
  return(fig)
}

#' Rotate 3D plotly object
#' 
#' @fig the plotly object
#' @sigma the divisor to get angle to rotate in each frame, bigger sigma mean slower rotation speed, default to 900 (180*5)
#' @return html widget that can be displayed in html format
#' @export
rotate_plotly <- function(fig, sigma = 900) {
  # custom jscode
  .jsCode = glue::glue("
      function(el, x){
  var id = el.getAttribute('id');
  var gd = document.getElementById(id);
  Plotly.update(id).then(attach);
  function attach() {
    var cnt = 0;
    
    function run() {
      rotate('scene', Math.PI / {%sigma%});
      requestAnimationFrame(run);
    } 
    run();
    
    function rotate(id, angle) {
      var eye0 = gd.layout[id].camera.eye
      var rtz = xyz2rtz(eye0);
      rtz.t += angle;
      
      var eye1 = rtz2xyz(rtz);
      Plotly.relayout(gd, id + '.camera.eye', eye1)
    }
    
    function xyz2rtz(xyz) {
      return {
        r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
        t: Math.atan2(xyz.y, xyz.x),
        z: xyz.z
      };
    }
    
    function rtz2xyz(rtz) {
      return {
        x: rtz.r * Math.cos(rtz.t),
        y: rtz.r * Math.sin(rtz.t),
        z: rtz.z
      };
    }
  };
}", .open = "{%", .close = "%}")
  
  # render the animation
  fig %>%
    htmlwidgets::onRender(jsCode = .jsCode)
}

#' Read clumpak result and admixture result in one object
#'
#' refactor all the information we need in a result object with class "clumpak_admixture"
#' @param prefix prefix of the results of admixture run saved in admixture outdir
#' @return list object with class "clumpak_admixture"
#' @examples
#' wrapper_read_clumpak_admixture(prefix)
#' @export
wrapper_read_clumpak_admixture <- function(prefix) {
  # prefix = "output/adm_ga2k/ga2k_adm.autosomes.common_rare_phased.maf0.01.indep_500kb_0.2.pruned"
  # read cv error from list of log file of admixture run
  cv_error = wrapper_read_log(prefix)
  
  # read number of run and number of run in major cluster
  df_nrun = wrapper_read_nrun(prefix)
    
  # read clumpak major last
  clumpak_major = map(df_nrun$clumpak, \(file) {
    # file = qfiles[5]
    # read clump output with function because qlist convert require a matrix of proportion only
    cmd = paste("cat", file, "| tr -s ' ' | cut -d ' ' -f6-")
    df = fread(cmd = cmd, header = FALSE)
    # rename column to avoid adding leading zero later
    colnames(df) = paste0("Cluster", 1:ncol(df))
    return(df)
  }) |>
    as.qlist()
  
  # combine result
  result <- list(
    clumpak_major = clumpak_major,
    cv_error = cv_error,
    df_nrun = df_nrun
  )
  
  #
  class(result) <- "clumpak_admixture"
  return(result)
}

#' Wrap function to parse number of run and number of run in major cluster in clumpak result
#' 
#' @param prefix prefix of the admixture output, e.g. example.k5r1.Q
#' @return df of info from the run
#' @export
wrapper_read_nrun <- function(prefix) {
  # prefix = "output/adm_ga2k/ga2k_adm.autosomes.common_rare_phased.maf0.01.indep_500kb_0.2.pruned"
  # check if we have clumpak output saved anywhere in wd
  cmd = paste0(
    "find output/clumpak/clumpakOutput -type f -name 'clusterFiles' -exec grep -l '",
    basename(prefix),
    "' {} +"
  )
  files = system(cmd, intern = TRUE)
  
  # check if we got back any files from the search
  if (length(files) > 0) {
    # get dir from files
    dir = files[1] |> str_replace("K=\\d+.*", "")
    # parse avail k from files
    ks = str_extract(files, "K=\\d+") |> parse_number() |> unique() |> sort()
    
    # declare qfils from output of clumppak from MajorCluster folder
    qfiles = paste0(dir, "/K=", ks, "/MajorCluster/CLUMPP.files/ClumppIndFile.output")
    
    # read number of run for each k
    nrun = map_int(ks, \(k) {
      file = paste0(dir, "/K=", k, "/CLUMPP.files/FilesToIndex")
      if (file.exists(file)) {
        R.utils::countLines(file) |> as.integer()
      } else {
        return(0)
      }
    }) |> set_names(ks)
    
    # read major cluster support rate from number of line in file
    nrun_major = map_int(ks, \(k) {
      file =str_c(dir, "/K=", k, "/MajorCluster/clusterFiles")
      R.utils::countLines(file)
    })
    
    # combine result
    df_nrun = tibble(
      k = ks,
      nrun = nrun,
      nrun_major = nrun_major,
      clumpak = qfiles
    )
  } else {
    stop("stop for now, we need clumpak result to sync cluster across K")
  }
  #
  return(df_nrun)
}

#' wrap function to read log file of admixture run
#' 
#' change priority to log file
#' @param prefix prefix of the file output
#' @return data frame of cv error and supporting information
#' @examples
#' prefix = "output/adm_ga2k/ga2k_adm.autosomes.common_rare_phased.maf0.01.indep_500kb_0.2.pruned"
#' wrapper_read_log(prefix)
#' @export
wrapper_read_log <- function(prefix) {
  # suffix of interest parsed from avail log file # prefix = params$prefix_auto
  files = fs::dir_ls(path = dirname(prefix), regexp = str_glue("{basename(prefix)}.*log")) |>
    # sort by numeric version
    str_sort(numeric = TRUE)

  # map files to parse_adm_log
  df.cv = map_dfr(files, ~parse_adm_log(.x))
  
  # return the data frame
  return(df.cv)
}

#' Funtion to parse log file of admixture program
#' 
#' @param file path to log file
#' @return data frame of cv error and supporting information
#' @examples
#' file = "output/adm_ga2k/ga2k_adm.autosomes.common_rare_phased.maf0.01.indep_500kb_0.2.pruned.k9r1.log"
#' parse_adm_log(file)
#' @export
parse_adm_log <- function(file) {
  # read all lines
  lines = read_lines(file)
  
  # parse fold
  fold = str_subset(lines, "Folds") %>% str_extract("Folds=\\d+") |>
    parse_number()
  
  # parse seed
  seed = str_subset(lines, "seed") |>
    parse_number()
  
  # parse cv error with condition
  df = str_subset(lines, "^CV") %>% fread(text = .)
  # check if the line start with CV exist # some run may error out with no CV line
  if (nrow(df) == 0) {
    return(NULL)
  } else {
    # if there is at least one row, then convert to tibble
    df = df |> select(3, 4) |>
      rename(K = 1, cv = 2) |>
      mutate(K = parse_number(K))
  }
  
  # create data frame to store output
  df.out = tibble(
    seed = seed,
    fold = fold
  ) |>
    bind_cols(df)
  #
  return(df.out)
}

#' Arrange ancestral components from admixture output
#' 
#' Add feature to rank the country by the second highest component within the region
#' Add feature to rank the ethnicity by the second highest component within the country
#' Add custom step to aggregate ethnicity with similar ancestral component together
#' @param df data frame of admixture output
#' @param df.meta2 data frame of metadata that have the same order as df
#' @param custom_swap conditional, default to FALSE now, can fix later 
#' @return data frame of all identifiers and sorted index column
#' @examples
#' arrange_ancestral_component(df = results$clumpak_major[[5]], meta = df.meta2)
#' @export
arrange_ancestral_component <- function(df, meta, custom_swap = FALSE) {
  # df = qlist[[9]] # qlist already worked, no need for debug
  # start with the known region
  if (!exists("region.known")) {
    stop("Declare 'region.known' first")
  }
  
  # assuming meta contain corresponding metadata for each sample
  # annotate df and mark the original index of sample
  df2 = bind_cols(df, meta) |> as_tibble() |>
    mutate(index = row_number())

  # declare identifiers for group
  identifiers = c("region", "country", "ethnicity", "index")
  
  # stat by levels of identifiers
  ls.stat = map(seq_along(identifiers), \(i) {
    # use a subset of the identifiers
    identifiers = identifiers[1:i]
    
    # prominent cluster grouped by identifiers
    df = helper_mean_across_rank(df2, identifiers)
    comment(df) = str_c(identifiers, collapse = "_")
    
    #
    return(df)
  })
  names(ls.stat) <- map_vec(ls.stat, ~comment(.x))
  
  #
  thres_prop = 0.1
  #
  results = map(seq_along(identifiers), \(i) {
    # use a subset of the identifiers
    identifiers = identifiers[1:i]
    
    # use the corresponding precomputed stat of the upper level
    if (i>1) {
      # corresponding component with upper level stats
      upper.stat = ls.stat[[i - 1]] |>
        filter(rank == 1) |>
        select(1:component) |>
        rename(sel_cluster = component)
      
      # keys to joint stat by
      upper.keys = colnames(upper.stat) |> str_subset("sel", negate = TRUE)
      
      #
      current.stat = ls.stat[[i]]
      
      # compute mean value across cluster group by identifiers
      result = current.stat |> ungroup() |>
        # annotate the selected cluster based on the upper.stat
        left_join(upper.stat, by = upper.keys) |>
        # consider only the row where component == sel_cluster
        filter(component == sel_cluster) |>
        # arrange group by upper identifiers and recompute rank column
        group_by_at(upper.keys) |>
        mutate(rank = dense_rank(value)) |>
        arrange(desc(rank), .by_group = TRUE)
      
      # add feature to ignore rank when value < thres_prop # thres_prop = 0.1
      result = result |>
        mutate(rank = if_else(value < thres_prop, 0, rank))
      
      need_follow_up = {filter(result, rank == 0) |> nrow()} > 0
      if (need_follow_up) {
        # follow up with the second rank of upper stat
        df.dict = ls.stat[[i - 1]] |>
          filter(rank == 2) |>
          select(1:component) |>
          rename(sel_cluster = component) |>
          distinct() |>
          # match with only rank = 0
          mutate(rank = 0)
        
        #
        df.dict = filter(result, rank == 0) |>
          # update sel_cluster based on the rest of the column in df.dict
          rows_update(df.dict, by = c(upper.keys, "rank"), unmatched = "ignore") |>
          mutate(component = sel_cluster) |>
          rows_update(current.stat, by = c(identifiers, "component"), unmatched = "ignore") |>
          # recompute rank
          ungroup() |> group_by_at(upper.keys) |>
          mutate(rank = -dense_rank(value))
        
        # change sel_cluster of rank == 0
        result = result |>
          rows_update(df.dict, by = identifiers, unmatched = "ignore")
      }
    } else {
      result = tibble(region = region.known) |>
        mutate(rank = row_number() |> rev())
    }
    return(result)
  })
  
  # right now, the results is a list of data frame # merge the ranking
  results = map(seq_along(results), \(i) {
    lookup = c("rank") |>
      set_names(stringr::str_glue("rank_l{i}"))
    results[[i]] |>
      dplyr::rename(all_of(lookup)) |>
      select(!any_of(c("component", "value", "sel_cluster")))
  })
  
  # merge the list of df into 1
  df.out = rev(results) |>
    purrr::reduce(.f = ~{
      identifiers = colnames(.y) |> str_subset("rank", negate = TRUE)
      left_join(.x, .y, by = identifiers)
    })
  
  #
  arrange_by_vars = colnames(df.out) |> str_subset("rank") |>
    str_sort(numeric = TRUE)
  df.out = df.out |>
    arrange(across(all_of(arrange_by_vars), desc))
  
  # assuming we want to add a custom step to change order of ethnicity here
  if (custom_swap) {
    # ls.stat[[3]] |>
    #   filter(country =="Russia", region == "NorthAsia") |>
    #   pivot_wider(id_cols = "ethnicity",names_from = component, values_from = value) |>
    #   arrange(desc(Cluster9))
    # select cluster that represent the Yupik
    component_PaleoSiberian = ls.stat[[3]] |> filter(ethnicity == "Yupik", value > 0.9) |>
      pull(component) |> first()
    # select ethnicity that have prominent component_PaleoSiberian
    PaleoSiberians = ls.stat[[3]] |> filter(component == component_PaleoSiberian, value > 0.5) |>
      pull(ethnicity)
    
    # from left to right: c("Khakass", "Shors", "TomskTatar", "Nivkh")
    # Turkic ancestry except the Nivkh
    df.turkish = df.out |> ungroup() |>
      filter(ethnicity %in% c("Khakass", "Shors", "TomskTatar", "Nivkh"))
    
    # remove the turkish ancestry from df.out and rows_insert them before all idx_PaleoSiberians
    df.tmp = df.out |> anti_join(df.turkish, by = "index") |> ungroup()
    
    # find the row index of NativeBeringian in df.tmp
    idx_PaleoSiberians = which(df.tmp$ethnicity %in% PaleoSiberians)
    
    # insert rows at a specific index in a dataframe
    df.out = df.tmp |>
      dplyr::slice(1:{min(idx_PaleoSiberians) - 1}) |>
      bind_rows(df.turkish) |>
      bind_rows(df.tmp |> slice(min(idx_PaleoSiberians):n()))
  }

  # return of function here
  return(df.out)
}

#' Helper function to compute across component
#' 
#' Assuming it is a df, grouped by variable of interest
#' @param df df that we want to compute
#' @param voi variable of interest
#' @param sel.rank select rank that we want to return
#' @return df
#' @export
helper_mean_across_rank <- function(df, voi = "country", sel.rank = 1:3) {
  result = df |>
    group_by_at(voi) |>
    # summarize mean across multiple column
    summarise(across(starts_with("Cluster"), mean), .groups = "drop") |>
    # pivot longer
    pivot_longer(cols = starts_with("Cluster"), names_to = "component", values_to = "value") |>
    # group by voi
    group_by_at(voi) |>
    # rank the component
    mutate(rank = rank(desc(value), ties.method = "first")) |>
    # arrange by value
    arrange(desc(value))
  
  #
  if (!missing(sel.rank)) {
    result = result |>
      # select the highest component
      slice(sel.rank)
  }
  return(result)
}

#' Wrapper for adm plot for qlist object
#' 
#' this wrapper is to select one k for plotting stacked bar chart, make "clumpak_admixture" the input
#' @param results take object of class "clumpak_admixture" as input
#' @param k assumed number of components
#' @param meta corresponding metadata for label
#' @param subset a list of named character vector to subset the data
#' @param custom_swap conditional on swapping cluster, default to FALSE
#' @param ... parameters to be passed to plotQ in list format
#' @return ggplot object
#' @examples
#' wrapper_adm_plot_k(results, 14:15, df.meta2)
#' @export
wrapper_adm_plot_k <- function(results, k = 11, meta, subset, custom_swap = FALSE, ...) {
  # check input
  if (class(results) != "clumpak_admixture") {
    stop("results must be of class 'clumpak_admixture'")
  }
  
  # extract the selected k
  ls.plot = {map_int(results$clumpak_major, ~ncol(.x)) %in% k} |> which() %>%
    results$clumpak_major[.]
  
  # custom sorting the sample based on the maximum k ----
  # max(k) is required to supplied appropriate legend_labs
  # get the sorted table # with sorted "index" column
  df.group.sort = arrange_ancestral_component(ls.plot[[length(ls.plot)]], df.meta2, custom_swap = custom_swap)

  # sorted indices of sample to plot
  indices.sort = df.group.sort$index
  
  # prepare the data frame of group label for plot qlist
  df.grplab = df.meta2 |>
    # select label of interest to plot
    select(any_of(c("ethnicity", "country", "region", "lngg.family"))) |>
    # sorted the table based on the index of df.group.sort
    slice(indices.sort)
  
  # apply new indices on ls.plot
  ls.plot = map(ls.plot, \(df) {
    # reorder the data frame based on index of df.group.sort
    df = df[indices.sort, ]
  })
  
  # ---- color and label vector ----
  # known colors from region map # just a reminder
  .colors.region = c(Africa = "black", America = "red", Europe = "lightsteelblue", 
                     Oceania = "yellow4", CentralAsia = "magenta3", EastAsia = "green4", 
                     NorthAsia = "royalblue", SouthAsia = "firebrick", SoutheastAsia = "yellow2", 
                     WestAsia = "lightslategrey", `NA` = "white")
  
  # declare colors vectors
  .colors = pals::cols25()
  
  # label vector for strip panel on the left side of the vertical axis
  strip_panel_labels = results$df_nrun |>
    filter(k %in% .env$k) |>
    mutate(label = str_c("K=", k, "\n", nrun_major, "/", nrun, "run")) |>
    pull(label)
  
  # declare label vectors for legend key
  legend_labs = paste0("Cluster", 1:max(k)) |>
    set_names(paste0("Cluster", 1:max(k)))
  
  # Switch the cluster order systematically ----
  # the column on the left side will appear first on legend and on top on the stacked bar chart
  # conditional swap
  if (custom_swap) {
    #
    ls.plot = map(ls.plot, \(.x) {
      # right now: KalmykNivkhKhakass and SenoiOrangAsli are swapped between 11:13 and 14:15 # swap them only when k>=14
      if (ncol(.x) >= 14) {
        # .x = ls.plot[[4]]
        # use df.cmp from k=15 to swap for both k=14 and k=15
        df.cmp = helper_reorder_component(ls.plot, meta = df.grplab)
        idx_to_swap = df.cmp %>%
          filter(as.character(key) %in% c("KalmykNivkhKhakass", "SenoiOrangAsli")) |>
          pull(component) |> parse_number() |> unique()
        
        # sanity check if we get 2 indices
        if (length(idx_to_swap) == 2) {
          # only do this if we got 2 indices # turn to tibble because data.table object is quirky
          .x = as_tibble(.x)
          .x[, idx_to_swap] <- .x[, rev(idx_to_swap)]
        }
        # this should not change the colnames
      }
      # return the .x
      return(.x)
    })
  }
  
  # legend labs need to be annotated again after the swap happen
  # pophelper draw the component by the order of column in the dataframe
  # compute the main cluster for each named component
  df.cmp = helper_reorder_component(ls.plot, meta = df.grplab)
  
  # reorder the legend_labs based on the order of df.cmp # legend key will be ordered based on legend_labs
  indices = df.cmp$component |> map_vec(~{which(legend_labs == .x)})
  legend_labs = c(legend_labs[indices], legend_labs[-indices])
  legend_labs[1:nrow(df.cmp)] <- df.cmp$label
  
  # declare 12 chosen .colors index for unique ancestral components
  ordered_components = levels(df.cmp$key)
  # color vector named by ordered ancestral components # using indices from the pals::cols25 itslef
  indices = c(19, 20, 1, 8, 2, 3, 5, 7, 4, 17, 10, 22, 6)
  # sometimes adding a named cluster break thing, make a check here
  if (length(indices) != length(ordered_components)) {
    stop("The number of ordered_components and indices must be the same")
  }
  color_names = {length(indices) + 1}:length(.colors) %>% str_c("Cluster", .) %>%
    c(ordered_components, .)
  # set colnames with color_names here will not work, vector names need to be "Cluster%d"
  .colors = c(.colors[indices], .colors[-indices])
  
  # remove the one that do not appear in df.cmp$key
  indices = which(!ordered_components %in% as.character(df.cmp$key))
  if (length(indices) > 0) {
    .colors = .colors[-indices]
  }
  
  # rewrite plotQ with rlang::exec ----
  # switch of single run or multi run
  if (length(k) > 1) {type = "join"} else {type = "sep"}
  
  # use ... to pass the arguments to custom_plotQ
  optional_arg = list(...)
  
  # Define the base arguments to pass to custom_plotQ
  args_list <- list(
    # imgoutput=join for multi k plotting
    qlist = ls.plot, imgoutput = type, returnplot = TRUE, exportplot = FALSE,
    # supply 25 colors from pals package # because we have up to k=20
    clustercol = .colors,
    # strip panel params on the y axis parameters
    sppos = "left", splab = strip_panel_labels, showyaxis = FALSE, splabsize = 6,
    # group label params
    showgrplab = TRUE, grplab = df.grplab, selgrp = "region", ordergrp = FALSE,
    grplabsize = 2, grplabpos = 0.6, grplabjust = 0.8, grplabheight = 10, grplabangle = 35,
    # panel ratio should be based on the number of annotated fields in grplab
    panelratio = c(ncol(df.grplab), 1), linesize = 0.8, pointsize = 4,
    # legend params
    showlegend = TRUE, legendlab = legend_labs, legendpos = "left",
    legendrow = 2,
    legendtextsize = 6, legendkeysize = 4
  )
  
  # If optional_arg is supplied, modify the args_list
  if (!is.null(optional_arg)) {
    # update value from args_list with optional_arg
    args_list <- modifyList(args_list, optional_arg)
  }
  
  # xlim for coord cartesian ----
  # move here because subset before reorder is not helpful
  if (!missing(subset)) {
    # subset = list(region = c("Europe", "CentralAsia", "NorthAsia"))
    key = names(subset)
    values = unlist(subset)
    
    # df.grplab and ls.plot will have the same index now
    idx.keep = {df.grplab[[key]] %in% values} |> which()
    
    #
    xlim = c(min(idx.keep), max(idx.keep))
    
    # update value of args_list with xlim
    args_list$xlim = xlim
  }
  
  # Use rlang::exec to dynamically call custom_plotQ with the args_list
  result <- rlang::exec(custom_plotQ, !!!args_list)
  # test grid.arrange(result$plot[[1]])
  # and add another element to the result list
  result$df_sort = df.group.sort
  
  # plotQ want to return a list of result, we also need it for inspecting the location of the label
  return(result)
}

#' helper function to reorder the component in admixture results
#' 
#' @param ls.plot list of data frame of admixture results
#' @param meta data frame of corresponding to qlist
#' @return df of label ordered with the prominent component
#' @examples
#' (df.cmp = helper_reorder_component(ls.plot, meta = df.grplab))
#' @export
helper_reorder_component <- function(ls.plot, meta) {
  # rank based on the highest K
  df = ls.plot[[length(ls.plot)]]
  
  # find the component based on table of criteria
  df.crt = tibble(
    label = c("WestEurasia", "NativeNorthRussia", "NativeBeringian", "HanChinese",
              "Austronesian", "Dravidian", "OrangAsli", "Andamanese", "Oceania",
              "TibetoBurman", "AustroAsiatic", "SenoiOrangAsli", "KalmykNivkhKhakass"),
    region = c("WestAsia", NA, NA, NA, NA, NA, NA, NA, "Oceania", NA, NA, NA, NA),
    ethnicity = c(NA, "Evenk", "Yupik", "NorthernHan", NA, NA, "Kensiu", "Onge", NA, NA, "Bru", "MahMeri", "Kalmyk"),
    country = c(NA, NA, NA, NA, "Malaysia", "India", NA, NA, NA, "India", NA, NA, NA),
    lngg.family = c(NA, NA, NA, NA, "Austronesian", "Dravidian", NA, NA, NA, "Tibeto-Burman", NA, NA, NA)
  )
  
  # find the most prominent component according to criteria
  result = map_dfr(split(df.crt, 1:nrow(df.crt)), \(.row){
    # .row = split(df.crt, 1:nrow(df.crt))[[1]]
    # remove column with NA value
    crt = .row |> select(!where(is.na) & !label)
    
    # get the index of the corresponding metadata to match calculation with df / ls.plot
    df.idx = meta |> mutate(index = row_number()) |>
      semi_join(crt, by = setNames(colnames(crt), colnames(crt)))
    
    # check if the criteria return any row
    if (nrow(df.idx) == 0) {
      return(NULL)
    } else {
      # find the prominent component
      cmp = df |> slice(df.idx$index) |>
        summarise(across(starts_with("Cluster"), mean)) |> unlist() |>
        sort(decreasing = TRUE) |> names() |> first()
    }
    
    # return the label with the component
    df.out = select(.row, label) |> mutate(component = cmp)
    return(df.out)
  })
  
  # reorder the result based on the order of label
  ordered_components = c("WestEurasia", "NativeNorthRussia", "NativeBeringian", "KalmykNivkhKhakass",
    "Dravidian", "HanChinese", "TibetoBurman",
    "AustroAsiatic", "OrangAsli", "SenoiOrangAsli", "Austronesian", "Andamanese", "Oceania")
  result = result |>
    mutate(key = factor(str_replace(label, "_.*", ""), levels = ordered_components)) |>
    arrange(key)
  
  # double check if there are identical component
  if (any(duplicated(result$component))) {
    # collapse duplicated component
    result = result |>
      group_by(component) |>
      summarise(label = paste(label, collapse = "_")) |>
      mutate(key = factor(str_replace(label, "_.*", ""), levels = ordered_components)) |>
      arrange(key)
  }
  return(result)
}

#' a

