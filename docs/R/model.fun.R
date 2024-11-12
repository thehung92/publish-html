#' aggregate model metrics of interest
#'
#' assuming the predicted class is stored in the validation set data frame
#' @param df data frame of the validation set
#' @param foi field of interest
#' @param foi.pred predicted field of interest based on the current model
#' @return a list of metrics
#' @examples
#' # example code
#' aggregate_model_metrics(df = valid)
#' @export
aggregate_model_metrics <- function(df, foi = "super_population", foi.pred = "cls.predicted") {
  # confusion matrix
  ls.cm = caret::confusionMatrix(
    data = df[[foi.pred]],
    reference = df[[foi]]
  )

  # df.metrics
  fois = c("Prevalence","Precision", "Recall", "F1")
  df.metrics = as.data.frame(ls.cm$byClass) %>%
    select(any_of(fois))

  # return
  list(cm = ls.cm$table,
       accuracy = ls.cm$overall["Accuracy"],
       metrics = df.metrics)
}

#' extract dfplot for tidymodels
#'
#' @param results results in tibble format from workflow_map
#' @param mode a string, either "tune" or "fit"
#' @return data frame to plot
#' @export
extract_dfplot <- function(results, mode = "tune") {
  # extract workflow ids
  ids = select(results, any_of("wflow_id")) |> unlist()

  # check mode
  if (mode == "tune") {
    map2(results$result, ids,~ {
      .x$.metrics[[1]] %>%
        mutate(model = .y) %>%
        # remove redundant cols
        select(-.estimator) %>%
        # separate preproc and model
        separate(model,
                 into = c("preproc", "model"),
                 sep = "_", remove = FALSE) %>%
        # pivot wider # the estimate metrics
        pivot_wider(names_from = .metric,
                    values_from = .estimate) %>%
        # arrange the data
        arrange(desc(accuracy))
    }) %>%
      bind_rows()
  } else if (mode == "fit") {
    # start from the data of autoplot function for tuning result
    apply(results, 1, function(row) {
      autoplot(row$result)
    }) %>%
      # combine the results into a single dataframe
      purrr::map2(., ids,~ {
        .x$data %>%
          mutate(model = .y)
      }) %>% bind_rows() %>%
      separate(model,
               into = c("preproc", "model"),
               sep = "_", remove = FALSE)
  } else {
    stop("mode not recognized")
  }
}

#' wrapper for ggforce plot for dimension reduction result
#'
#' Add option to choose the outcome variable
#' @param df.plot data frame of the result
#' @param outcome string for the outcome variable
#' @return ggplot object
#' @export
wrapper_ggforce_dimension_reduction <- function(df.plot, outcome = "super_population") {
  # elapse time to print
  elapse = attributes(df.plot)$runtime["elapsed"] %>%
    round(0) %>% hms::as_hms()

  # PCA score with ggforce # and visualize the 1st 4 pc
  ggplot(df.plot, aes(x = .panel_x, y = .panel_y, color = !!sym(outcome), fill = !!sym(outcome))) +
    geom_point(alpha = 0.4, size = 0.5) +
    ggforce::geom_autodensity(alpha = .3) +
    ggforce::facet_matrix(
      vars( -!!sym(outcome)),
      layer.diag = 2,
    ) +
    scale_color_brewer(palette = "Dark2", guide = guide_legend(
      title = str_c("Runtime (hms):\n", elapse, "\n\n", outcome),
      override.aes = list(size = 2, alpha = 1, fill = NA, linewidth = 0))
    ) +
    scale_fill_brewer(palette = "Dark2", guide = "none")
}

#' Basic outlier detection from select columns
#'
#' sample is outlier if at least one of the columns is outlier
#' outlier is defined as 6 mean absolute deviation away from the median
#' OR 6 standard deviation away from the mean
#'
#' @param df data frame
#' @param vois variable of interest
#' @param method mean or median, default to median
#' @return a logical vector of outlier, same order as the input data frame
#' @export
basic_outlier_detect <- function(df, vois, method = "median") {
  # check if df colnames contain vois
  if (!all(vois %in% colnames(df))) {
    stop("Not all vois found in the data frame")
  }

  # check method
  if (method == "mean") {
    vt.is_outlier = select(df, all_of(vois)) |>
      # apply computation for each column to compute outlier for each sample / PC
      apply(2, function(x) {
        {abs(x - mean(x)) / sd(x)} > 6
      }) |>
      # apply OR condition to each row
      apply(1, function(x) reduce(x, `|`)) |>
      unlist()
  } else if (method == "median") {
    vt.is_outlier = select(df, all_of(vois)) |>
      # apply computation for each column to compute outlier for each sample / PC
      apply(2, function(x) {
        {abs(x - median(x)) / mad(x)} > 6
      }) |>
      # apply OR condition to each row
      apply(1, function(x) reduce(x, `|`)) |>
      unlist()
  } else {
    stop("method not recognized")
  }

  #
  return(vt.is_outlier)
}
