#!/usr/bin/env Rscript
# ---- RATIONALE ----
# scatter3d plot for MDS on merge clean cross-impute with K-nearest-neighbors
# start with chr22
# ---- Setup env ----
# set seed for reproducibility
set.seed(42)

# library
packages <- c(
  "argparse", # command line argument parser, more features than base-R
  "tidyverse", # multi-purpose data-wrangling # and ggplot
  "data.table", # read and write data multi-threading # with easier to remember syntax
  "janitor", # clean up name
  "openxlsx", # handling excel files # with more feature than readxl
  "paletteer", # common API for multiple color palette
  "scales", # label (date,time,currency,percentage) and axis break
  "patchwork", # combine and annotate figure
  "plotly",
  # "reticulate",
  "tidymodels"
)

# load library quietly and stop if library can not be loaded
for (package in packages) {
  if (suppressPackageStartupMessages(require(package, character.only = TRUE))) {
  } else {
    stop("install required packages before running script")
  }
}

# resolve conflict
tidymodels_prefer()

# source functions if defined in separate file
files = c("R/structure_analysis.fun.R", "R/mds_knn.fun.R", "R/model.fun.R")
for (file in files) {
  if (file.exists(file)) {
    source(file)
  } else {
    warning(paste(file, "not found."))
  }
}

# ---- Parse cmd line args ----
# parse command line arguments
parser <- ArgumentParser(description='Program role')
parser$add_argument('prefix', nargs='?', help = 'input file path')
parser$add_argument('--colors', type='character', nargs='?', default = "data/select_color_population.tsv",
                    help = 'file contain color palette')
parser$add_argument('--meta', type='character', nargs='?',
                    help = 'files contain info from label cleaning process')
parser$add_argument('--sqlite', type='character', nargs='?', default = "/Users/hung/Data/main-project/tasks/20241019_sqlite/data/ga5k-v1.2.sqlite",
                    help='sqlite database that contain info for ga100k samples')
parser$add_argument('-o', '--outfile', type='character',
                    help='output file path')

# parse arguments
if (interactive()) {
  # INTERACTIVE TEST WHEN NEEDED:
  args <- parser$parse_args(c("/Users/hung/Data/main-project/tasks/classificationModel/output/merge/merge_clean.chr22.bisnp.phased.impute_cross.maf0.001.exclhighld.indep-pairwise-500kb-0.2.pruned.ibd.mds",
                              "-o", "man/figs/mds3d-merge-clean-cross-impute",
                              "--meta", "/Users/hung/Data/main-project/tasks/classificationModel/manuscript/figs/*clean_knn.tsv"))
} else {
  args <- parser$parse_args()
}

# ---- metadata ----
# fam of pruned file set corresponding to mds file
fam = args$prefix |>
  str_replace("pruned.*", "pruned.fam")
fam = read_tsv(fam, show_col_types = FALSE,
         col_names = c("fid", "iid", "pid", "mid", "sex", "pheno"))

# read metadata from label clean process
files = str_glue("ls -1v {args$meta}") |> system(intern = TRUE)
df.meta = map_dfr(files, ~{
  # filter high confidence sample
  df = read_tsv(.x, show_col_types = FALSE) |>
    filter(grepl("high_prob", type_knn))
  
  # relabel KHV and CDX to SoutheastAsia
  if ("pop" %in% colnames(df)) {
    df = df |> 
      mutate(region = case_when(
        pop %in% c("KHV", "CDX") ~ "SoutheastAsia",
        .default = region
      )) |>
      mutate(source = "1kgp3")
  } else {
    df = df |>
      mutate(source = "ga100k") |>
      # remove tibeto burman ancestry in india
      filter(! ethnicity %in% c("Chakma", "Jamatia", "Mog", "Toto", "RanaTharu"))
  }
  
  # select only IID and region to plot
  df |>
    select(IID, region, source)
})
# let's check if fixing region of KHV and CDX is enough for separate cluster
# annotate fam$iid with metadata
df.meta = tibble(IID = fam$iid) |>
  left_join(
    df.meta,
    by = "IID"
  )

# read color for region
df.color = read_tsv(args$colors, show_col_types = FALSE) |>
  filter(type == "region")

# ---- MDS result ----
# default to read only the first 3 components but compute pves for all
results = read_mds_result(args$prefix) |>
  # then scale the score to pves*100
  scale_score()

# tibble(x=1:10, y = results$pves[1:10]*100) |>
#   ggplot(aes(x=as.factor(x), y=y)) + geom_point()
  

# annotate with ga100k metadata
results$score = results$score |>
  left_join(df.meta, by = c("IID" = "IID")) |>
  # and remove sample without region
  filter(!is.na(region))

# data frame to plot
df.plot = results$score |>
  # create hoverinfo column for plotly
  mutate(hoverinfo = paste0(
    IID, ";", region, ";", source
  ))

# annotate region column with nsample
dict = df.plot |>
  count(region) |>
  mutate(region_n = str_glue("{region}({n})"))
df.plot = df.plot |>
  left_join(dict, by = "region")

# use plotly to plot 3d scatter with C1,C2,C3 and colored by "region" with colors converted from df.color
colors = df.color |> filter(label %in% unique(df.plot$region)) |> pull(color)
names(colors) = df.color |> filter(label %in% unique(df.plot$region)) |> pull(label)

# plot and save as html
# choose component to visualize
if (comment(results) == "pca") {
  vois = c("PC1", "PC2", "PC3")
} else {
  vois = c("C1", "C2", "C3")
}
# compute pves to annotate axis title
axis_titles = map_vec(vois, ~{
  idx = .x |> parse_number()
  pct = results$pves[idx] |> percent(accuracy = 0.01)
  str_c(.x, " (", pct, ")")
})
#
widget = plot_ly(data = df.plot) |>
  add_trace(
    type = "scatter3d", mode = "markers",
    x = ~.data[[vois[1]]], y = ~.data[[vois[2]]], z = ~.data[[vois[3]]],
    color = ~region, colors = colors,
    text = ~hoverinfo, hoverinfo = "text", name = ~region_n
  ) |>
  plotly::layout(
    legend = list(
      title = list(text = str_glue("<b>Region</b>"))),
    scene = list(
      # # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
      # camera = list(eye = list(x = 1.25, y = -1.1, z = 0.8),
      #               center = list(x = 0, y = 0, z = 0)
      # ),
      # add axis title
      xaxis = list(title = axis_titles[1]),
      yaxis = list(title = axis_titles[2]),
      zaxis = list(title = axis_titles[3])
    )
  ) ; widget

# Simulate a conda environment to use Kaleido # path=/Users/hung/Library/r-miniconda-arm64
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_python('/Users/hung/Library/r-miniconda-arm64/envs/r-reticulate/bin/python', required = TRUE)
reticulate::use_miniconda('r-reticulate')

# Save the figure pdf
pdf = paste0(args$outfile, ".pdf")
save_image(widget, file = pdf, height=500, width=700) # browseURL(pdf)

# save html widget
html = paste0(args$outfile, ".html")
htmlwidgets::saveWidget(widget, html) # browseURL(html)
# save RData
rdata = paste0(args$outfile, ".RData")
obj = plotly_build(widget); save(obj, file = rdata)

# --- PCA results ----
# default to 10 PCs
prefix_pca = args$prefix |> str_replace("ibd.mds", "pca")
results = read_pca_result(prefix_pca) |>
  # then scale the score to pves*100
  scale_score()

# annotate with ga100k metadata
results$score = results$score |>
  left_join(df.meta, by = c("IID" = "IID"))

# data frame to plot
df.plot = results$score |>
  # create hoverinfo column for plotly
  mutate(hoverinfo = paste0(
    IID, ";", region, ";", source
  ))

# annotate region column with nsample
dict = df.plot |>
  count(region) |>
  mutate(region_n = str_glue("{region}({n})"))
df.plot = df.plot |>
  left_join(dict, by = "region")

# plot and save as html
# choose component to visualize
if (comment(results) == "pca") {
  vois = c("PC1", "PC2", "PC3")
} else {
  vois = c("C1", "C2", "C3")
}
# compute pves to annotate axis title
axis_titles = map_vec(vois, ~{
  idx = .x |> parse_number()
  pct = results$pves[idx] |> percent(accuracy = 0.01)
  str_c(.x, " (", pct, ")")
})
#
widget = plot_ly(data = df.plot) |>
  add_trace(
    type = "scatter3d", mode = "markers",
    x = ~.data[[vois[1]]], y = ~.data[[vois[2]]], z = ~.data[[vois[3]]],
    color = ~region, colors = colors,
    text = ~hoverinfo, hoverinfo = "text", name = ~region_n
  ) |>
  plotly::layout(
    legend = list(
      title = list(text = str_glue("<b>Region</b>"))),
    scene = list(
      # # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
      # camera = list(eye = list(x = 1.25, y = -1.1, z = 0.8),
      #               center = list(x = 0, y = 0, z = 0)
      # ),
      # add axis title
      xaxis = list(title = axis_titles[1]),
      yaxis = list(title = axis_titles[2]),
      zaxis = list(title = axis_titles[3])
    )
  ) ; widget



# ---- Supervised classification model ----
# build a K nearest neighbours on 5 dimension of mds # and scale score based on pves*100
results = read_mds_result(args$prefix, 5) |>
  scale_score()

# data for knn
data = results$score |>
  left_join(df.meta, by = c("IID" = "gaid")) |>
  select(c("region", starts_with("C"), "IID", "ethnicity", "country", "lngg_family")) |>
  # remove ethnicity == "unknown" from consideration
  filter(ethnicity != "unknown")

# # resample object with split of 4/5 train test
# data_rsplit = rsample::vfold_cv(data, v=5, strata = region)

# recipe
rec = recipes::recipe(head(data)) |>
  update_role(starts_with("C"), new_role = "predictor") |>
  update_role(region, new_role = "outcome")

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
  # and highlight discordant, concordant
  mutate(hli = if_else(region != pred_class, "discordant", "concordant")) |>
  # add hli low_prob
  mutate(conf = case_when(
    if_all(starts_with(".pred_"), ~.x < 0.8) ~ "low_prob",
    .default = "high_prob"
  )) |>
  # add hoverinfo
  mutate(hoverinfo = paste0(
    IID, ";", ethnicity, ";", country, ";", region, "\n",
    "pred_class: ", pred_class
  )) |>
  # add type_knn column to work with plotly
  unite("type_knn", c("hli", "conf"), sep = ";", remove = FALSE)

# ---- Active label cleaning ----
## workaround: label for knn classification
# annotate freq of sample hard to classify
dict = df.predict |> filter(type_knn != "concordant;high_prob") |>
  count(type_knn) |> mutate(label_knn = str_glue("{type_knn}({n})")) |>
  select(-n)
df.predict = df.predict |>
  left_join(dict, by = "type_knn")

# annotate freq of concordant sample
dict = df.predict |> filter(type_knn == "concordant;high_prob") |>
  count(region) |> mutate(label_high_conf = str_glue("{region}({n})")) |>
  select(-n)
df.predict = df.predict |>
  left_join(dict, by = "region")

## workaround: create one single color vector for all the layer
colors_combine = colors |>
  c("knn" = "grey")


# plotly 3d only accept one color vector for all the traces
widget = plot_ly(colors = colors_combine) |>
  # traces of symbols
  add_trace(
    data = df.predict |> filter(type_knn != "concordant;high_prob"),
    legendgroup = "knn", legendgrouptitle = list(text = "Difficult to classify"),
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
    legendgroup = "region", legendgrouptitle = list(text = "Region\n(#sample high confidence)"),
    type = "scatter3d", mode = "markers",
    x = ~C1, y = ~C2, z = ~C3,
    color = ~as.factor(region), colors = colors, name = ~label_high_conf,
    text = ~hoverinfo, hoverinfo = "text", inherit = FALSE
  ) |>
  # custom layout
  plotly::layout(
    scene = list(
      # z axis will control the camera eye on the vertical axis # this is important to rotate the camera
      camera = list(eye = list(x = -1, y = 1, z = 0.8),
                    center = list(x = 0, y = 0, z = 0)
      ),
      # add axist title
      xaxis = list(title = "C1"),
      yaxis = list(title = "C2"),
      zaxis = list(title = "C3")
    )
  ); widget

# Save the figure pdf
pdf = paste0(args$outfile, "_knn.pdf")
save_image(widget, file = pdf, height=500, width=700) # browseURL(pdf)

# save html widget
html = paste0(args$outfile, "_knn.html")
htmlwidgets::saveWidget(widget, html) # browseURL(html)

# save tsv
tsv = paste0(args$outfile, "_knn.tsv")
write_tsv(df.predict, tsv)


# ---- Resample to take care of imbalance class frequency ----
# manual inspection show that
# kazakh, uyghur, uzbek can be considered central asia
# 3 taiwan indigenous should be SoutheastAsia ancestry
# buriat are north asia ancestry
# Taz:Russia have genetic ancestry close to east asia # less than 1000 people left # descedant of intermarriage between Han and tungusic people

# data cleaning based on the previous run
# build a K nearest neighbours on 5 dimension of mds # and scale score based on pves*100
results = read_mds_result(args$prefix, 5) |>
  scale_score()

# data for knn along with annotation from metadata
col_hoverinfo = c("IID", "ethnicity", "country", "lngg_family")
data = results$score |>
  left_join(df.meta, by = c("IID" = "gaid")) |>
  select(c("region", starts_with("C"), all_of(col_hoverinfo))) |>
  # remove ethnicity == "unknown" from consideration
  filter(ethnicity != "unknown") |>
  # relabel Malay ethnic minority
  mutate(ethnicity = case_when(
    # Semang: dark skinned hunter gatherer #
    country == "Malaysia" & ethnicity == "Kintak" ~ "Kintaq-Semang",
    country == "Malaysia" & grepl("Kensiu|Batek|Mendriq|Jahai|Lanoh", ethnicity) ~ str_c(ethnicity, "-Semang"),
    # Senoi: largest group within OrangAsli, slash-burn agriculture
    country == "Malaysia" & ethnicity == "Jah Hut" ~ "JahHut-Senoi",
    country == "Malaysia" & ethnicity == "SenoiSemai" ~ "Semai-Senoi",
    country == "Malaysia" & ethnicity == "CheWong" ~ "CheqWong-Senoi",
    country == "Malaysia" & grepl("MahMeri|SemaqBeri|Temiar", ethnicity) ~ str_c(ethnicity, "-Senoi"),
    # ProtoMalay: lighter, farmer, malayic language
    country == "Malaysia" & grepl("Temuan|Jakun", ethnicity) ~ str_c(ethnicity, "-ProtoMalay"),
    # academic contexts: "KadazanDusun" is preferred as it recognizes the unity of these two closely related ethnic groups
    country == "Malaysia" & grepl("Kadazan|Dusun", ethnicity) ~ "KadazanDusun",
    # Dayak people # umbrella term for non-muslim native in Borneo
    country == "Malaysia" & grepl("Iban|Bidayuh", ethnicity) ~ str_c(ethnicity, "-Dayak"),
    .default = ethnicity
  ))

# data |> filter(country == "Malaysia") |> count(ethnicity) |> arrange(desc(n)) -> X

# count sample with prevalence < threshold
identifiers = c("ethnicity", "country", "region")
thres = 5
dict_prev_lt_thres = data |>
  group_by_at(identifiers) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(n) |>
  filter(n < thres)

# sample with low_prob from the first run of mds knn
iid_low_conf = df.predict |>
  filter(grepl("low_prob", type_knn)) |>
  pull(IID)

# run knn on scaled mds # 5 1st component and return 3d scatter plot with annotation
data_clean = data |>
  filter(country != "Taiwan") |>
  # no west asia
  filter(region != "WestAsia") |>
  # reclassify Uyghur into CentralAsia
  mutate(region = case_when(
    ethnicity == "Uyghur" ~ "CentralAsia",
    grepl("^Dai$|^CDX$", ethnicity) ~ "SoutheastAsia",
    ethnicity == "Buriat" ~ "NorthAsia",
    ethnicity == "Khalkh" ~ "NorthAsia",
    .default = region
  )) |>
  # remove sample with ethnicity, country prevalence less than threshold
  anti_join(dict_prev_lt_thres, by = identifiers) |>
  # also remove Onge Jarwa and Taz
  filter(ethnicity != "Jarwa") |> filter(ethnicity != "Taz") |>
  filter(ethnicity != "Onge") |> filter(ethnicity != "Nivkh") |>
  # Hmong have genetic structure very close to east asia
  filter(ethnicity != "Hmong") |>
  # mongolian have cluster between EA and NA. Closer to NA
  filter(ethnicity != "Mongolian") |>
  # Dungan Central Asia people have cluster close to eastasia
  filter(ethnicity != "Dungan") |>
  # exclude ambiguous Russian # which are close to european
  filter(ethnicity != "Russian") |>
  # exclude ethnic minority in north china
  filter(! grepl("^Oroqen$|^Xibo$|^Hezhen$", ethnicity)) |>
  # exclude Pakistan
  filter(country != "Pakistan") |>
  # exclude Lahu-China with 8 sample # more SEA
  filter(ethnicity != "Lahu") |>
  # Myanmar ethnicity are close with EA cluster # remove all except Burmese
  filter(!{country == "Myanmar" & ethnicity != "Burmese"}) |>
  # remove ethnic minority in Malaysia that can not represent the population
  filter(!{country == "Malaysia" & grepl("ProtoMalay|Rungus|Senoi|Semang", ethnicity)}) |>
  # remove ethnic minority that cannot represent the population of Thailand
  filter(!{country == "Thailand" & grepl("Bru|HtinPray|Moklen", ethnicity)}) |>
  # remove ethnic minority from China that contribute to SEA
  filter(!{country == "China" & grepl("CDX|Dai", ethnicity)}) |>
  # remove Bangladesh from consideration
  filter(country != "Bangladesh") |>
  # remove EA sample that can overlap with 1kgp3
  filter(!grepl("^Japanese$|CHB|CHS", ethnicity)) |>
  # remove sample of low_confidence from the initial df.predict
  filter(! IID %in% iid_low_conf)

# resample population with high prevalence
criteria = list(
  tibble(ethnicity = "Kinh", country = "Vietnam", n = 50),
  tibble(ethnicity = "Malay", country = c("Malaysia"), n = 50),
  tibble(ethnicity = "Malay", country = c("Singapore"), n = 50),
  tibble(ethnicity = "Chinese", country = "Singapore", n = 50),
  tibble(ethnicity = "Chinese", country = "Malaysia", n = 50),
  tibble(ethnicity = c("Kadazan", "KadazanDusun", "Dusun"), country = "Malaysia", n = 50),
  tibble(ethnicity = "Bajau", n = 50),
  tibble(ethnicity = "Iban-Dayak", country = "Malaysia", n = 50),
  tibble(ethnicity = "Bidayuh-Dayak", country = "Malaysia", n = 50),
  tibble(country = "Thailand", n = 100),
  tibble(region = "CentralAsia", n = 100),
  tibble(ethnicity = "Mansi", n = 10)
)

# map criteria to exclude outlier # then sample 100 sample
ids_exclude = map(criteria, \(df) {
  #
  df.sel = data_clean |> semi_join(df, by = colnames(df) |> str_subset("^n$", negate = TRUE))
  # quick check if any row is selected
  if (nrow(df.sel) == 0) {
    return(NULL)
  }

  # outlier if any of the 3 components is greater than 6 time median absolute deviation from median
  df.not_outlier = df.sel |>
    filter(! basic_outlier_detect(df.sel, str_c("C", 1:3)))

  # check if n > nrow(df.sel)
  if (df$n[1] > nrow(df.sel)) {
    return(NULL)
  } else {
    # sample_n
    df.resample = df.not_outlier |>
      sample_n(df$n[1])

    # return exclude list
    setdiff(df.sel$IID, df.resample$IID)
  }
}) |>
  unlist()

# resample population of high prevalence
data_resample = data_clean |>
  filter(! IID %in% ids_exclude)

#
ls.result = wrapper_mds_knn(data_resample, col_hoverinfo = col_hoverinfo)
ls.result$plot
ls.result$data |>
  filter(!grepl("low_prob", type_knn)) |>
  filter(!region %in% c("NorthAsia", "CentralAsia")) |>
  count(ethnicity, country, region) |>
  arrange(region ,desc(n)) -> X

# save dataset after clean sample of admixed population, resample high prevalence
prefix = args$outfile |> str_c("_clean_knn")
# Save the figure pdf
pdf = paste0(prefix, ".pdf")
save_image(ls.result$plot, file = pdf, height=500, width=700) # browseURL(pdf)
# save html widget
html = paste0(prefix, ".html")
htmlwidgets::saveWidget(ls.result$plot, html) # browseURL(html)
# save tsv
tsv = paste0(prefix, ".tsv")
write_tsv(ls.result$data, tsv)

# write sample list for concordant sample
txt = str_glue("tmp/sample_{basename(prefix)}.list")
ls.result$data |>
  filter(hli == "concordant", conf == "high_prob") |> pull(IID) |>
  write_lines(txt)
