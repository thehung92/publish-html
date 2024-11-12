#' wrapper around knn
#'
#' @param data data frame to be used for training
#' @param col_hoverinfo colnames to show on the hoverinfo besides outcome and pred_class
#' @param outcome colname of the outcome variable
#' @return htmlwidget
#' @examples
#' # example code
#' wrapper_mds_knn(data = data)
#' @export
wrapper_mds_knn <- function(data, outcome = "region", col_hoverinfo) {
  # recipe
  rec = recipes::recipe(head(data)) |>
    update_role(starts_with("C"), new_role = "predictor") |>
    update_role(!!sym(outcome), new_role = "outcome")

  # declare knn model
  knn_model = parsnip::nearest_neighbor(neighbors = tune()) |>
    parsnip::set_engine("kknn") |>
    parsnip::set_mode("classification")

  # Declare simple map of workflow to whole data to inspect discordance class
  wf = workflow(rec, knn_model |> set_args(neighbors = 10))

  # fit with data
  wf_fit = fit(wf, data)
  # extract the parsnip model fitted
  model_fit = wf_fit |> extract_fit_parsnip()

  # and use it to predict the same data
  df.predict = bind_cols(
    data,
    predict(model_fit, new_data = data) |> rename(pred_class = 1),
    predict(model_fit, new_data = data, type = "prob")
  ) |>
    # and highlight discordant, concordant between outcome and pred_class
    mutate(hli = if_else(!!sym(outcome) != pred_class, "discordant", "concordant")) |>
    # add hli low_prob
    mutate(conf = case_when(
      if_all(starts_with(".pred_"), ~.x < 0.8) ~ "low_prob",
      .default = "high_prob"
    )) |>
    # add type_knn column to work with plotly
    unite("type_knn", c("hli", "conf"), sep = ";", remove = FALSE)

  # add hoverinfo based on avail column exclude outcome and predictor
  if (missing(col_hoverinfo)) {
    col_hoverinfo = colnames(data) |>
      str_subset(outcome, negate = TRUE) |>
      str_subset("^C", negate = TRUE)
  }

  # make R expression for hoverinfo inside data
  expr = col_hoverinfo |>
    str_flatten(collapse = ", ';', ")
  expr = str_glue("paste0(
    outcome, ': ', .data[[outcome]], '\\n',
    'pred_class: ', pred_class, '\\n',
    {expr}
  )")
  #
  df.predict = df.predict |>
    mutate(hoverinfo = eval(parse(text = expr)))

  ## workaround: label for knn classification
  # annotate freq of sample hard to classify
  dict = df.predict |> filter(type_knn != "concordant;high_prob") |>
    count(type_knn) |> mutate(label_knn = str_glue("{type_knn}({n})")) |>
    select(-n)
  df.predict = df.predict |>
    left_join(dict, by = "type_knn")

  # annotate freq of concordant sample
  dict = df.predict |> filter(type_knn == "concordant;high_prob") |>
    count(!!sym(outcome)) |> mutate(label_high_conf = str_glue("{.data[[outcome]]}({n})")) |>
    select(-n)
  df.predict = df.predict |>
    left_join(dict, by = outcome)

  ## workaround: create one single color vector for all the layer
  colors_combine = colors |>
    c("knn" = "grey")

  ## workaround: axis title for percentage variance explained based on standard deviation
  sds = data |>
    summarise(across(matches("C\\d+"), ~sd(.x) |> round(2))) |>
    unlist()

  # improvised condition on sds # sds is sorted in descending order # OR all sds < 1
  if (is.unsorted(rev(sds)) | all(sds < 1)) {
    # the score is not scaled to mds yet
    axis_titles = str_glue("{names(sds)} (mean=0;sd={sds})")
  } else {
    # the score is scaled on pves*100
    axis_titles = str_glue("{names(sds)} ({sds}%)")
  }

  # plotly 3d only accept one color vector for all the traces
  widget = plot_ly(colors = colors_combine) |>
    # traces of symbols
    add_trace(
      data = df.predict |> filter(type_knn != "concordant;high_prob"),
      #
      legendgroup = "low_conf", legendgrouptitle = list(text = "Difficult to classify"),
      type = "scatter3d", mode = "markers",
      x = ~C1, y = ~C2, z = ~C3,
      # color = ~pred_class, colors = colors,
      symbol = ~type_knn, symbols = c("cross", "x"), opacity = 0.8,
      # size = 10,
      color = "knn",
      name = ~label_knn,
      text = ~hoverinfo, hoverinfo = "text"
    ) |>
    # overlay colors
    add_trace(
      data = df.predict,
      # groups legend entries so that clicking on one legend entry will hide or show all of the traces in that group
      legendgroup = "high_conf", legendgrouptitle = list(text = str_glue("{str_to_title(outcome)}\n(#sample high confidence)")),
      type = "scatter3d", mode = "markers",
      x = ~C1, y = ~C2, z = ~C3,
      color = ~.data[[outcome]], colors = colors, name = ~label_high_conf,
      text = ~hoverinfo, hoverinfo = "text", inherit = FALSE
    ) |>
    # custom layout
    plotly::layout(
      font = list(family = 'arial'),
      scene = list(
        # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
        camera = list(eye = list(x = -1, y = 1, z = 0.8),
                      center = list(x = 0, y = 0, z = 0)
        ),
        # add axis title # based the percentage of variance explained on data$C standard deviation
        xaxis = list(title = axis_titles[1]),
        yaxis = list(title = axis_titles[2]),
        zaxis = list(title = axis_titles[3])
      )
    )

  # consider return list of both widget and df.predict
  ls.result = list(
    plot = widget,
    data = df.predict
  )
  return(ls.result)
}
