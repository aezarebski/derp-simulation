library(hdf5r)
library(xtable)

sim_hdf5 <- "out/sim-charizard/dataset-charizard.hdf5"
db_conn <- H5File$new(sim_hdf5, mode = "r")
sim_names <- db_conn$ls()$name

## 1. read the config_json string attribute from the db_conn
## 2. parse it as JSON to get a list
## 3. write the list to a JSON file with the name "prior-simulation-config.json"
config_json_str <- h5attr(db_conn, "config_json")
config_list <- jsonlite::fromJSON(config_json_str, simplifyVector = TRUE)
jsonlite::write_json(config_list, path = "./src/prior-simulation-config.json", auto_unbox = TRUE, pretty = TRUE)

r0_vals_q <- "output/parameters/r0/values"
net_rem_rate_q <- "output/parameters/net_removal_rate/values"
sampling_prop_q <- "output/parameters/sampling_prop/values"
epi_dur_q <- "input/present"
final_prev_q <- "output/present_prevalence"
final_cuminf_q <- "output/present_cumulative"
wall_time_q <- \(g) h5attr(g,"simulation_wall_time")

r0_vals_vec <- c()
net_rem_rate_vec <- c()
sampling_prop_vec <- c()
epi_dur_col <- c()
final_prev_col <- c()
final_cuminf_col <- c()
wall_times_col <- c()

pb <- progress::progress_bar$new(format = " [:bar] :percent eta: :eta",
                                 total = length(sim_names),
                                 clear = FALSE, width= 60)
for (sx in sim_names) {
  pb$tick()
  rec <- db_conn[[sx]]
  r0_vals_vec <- c(r0_vals_vec, rec[[r0_vals_q]][])
  net_rem_rate_vec <- c(net_rem_rate_vec, rec[[net_rem_rate_q]][])
  sampling_prop_vec <- c(sampling_prop_vec, rec[[sampling_prop_q]][])
  epi_dur_col <- c(epi_dur_col, rec[["input/present"]][] )
  final_prev_col <- c(final_prev_col, rec[[final_prev_q]][] )
  final_cuminf_col <- c(final_cuminf_col, rec[[final_cuminf_q]][] )
  wall_times_col <- c(wall_times_col, wall_time_q(rec))
}

var_summary <- function(x) {
  list(
    lower = quantile(x, 0.025),
    median = median(x),
    upper = quantile(x, 0.975)
  )
}

sim_summary <- list(
  r0 = var_summary(r0_vals_vec),
  net_removal_rate = var_summary(net_rem_rate_vec),
  sampling_prop = var_summary(sampling_prop_vec),
  epi_duration = var_summary(epi_dur_col),
  final_prevalence = var_summary(final_prev_col),
  final_cumulative_infections = var_summary(final_cuminf_col),
  total_wall_time_secs = sum(wall_times_col)
)

summary_str <- function(var, n=2, as_e=FALSE) {
  if (as_e) {
    format_str <-gsub(pattern = "X", replacement = n, x = "$ %.Xe ~(%.Xe, %.Xe) $")
  } else {
    format_str <-gsub(pattern = "X", replacement = n, x = "$ %.Xf ~(%.Xf, %.Xf) $")
  }
  result <- sprintf(
    format_str,
    sim_summary[[var]][["median"]],
    sim_summary[[var]][["lower"]],
    sim_summary[[var]][["upper"]]
  )
  if (as_e) {
    result <- gsub(pattern = "e\\+00", replacement = "", x = result)
    result <- gsub(pattern = "e\\+01", replacement = " \\\\\\\\\\times 10^{1}", x = result)
    result <- gsub(pattern = "e\\+02", replacement = " \\\\\\\\\\times 10^{2}", x = result)
    result <- gsub(pattern = "e\\+03", replacement = " \\\\\\\\\\times 10^{3}", x = result)
    result <- gsub(pattern = "e\\+04", replacement = " \\\\\\\\\\times 10^{4}", x = result)
  }
  return(result)
}

## read the lines of the file "prior-latex-template.tex"
template_lines <- readLines("./src/prior-latex-template.tex")

## replace the patterns in the template with the summary statistics
template_lines <- gsub("summ1", summary_str("r0"), template_lines)
template_lines <- gsub("summ2", summary_str("net_removal_rate"), template_lines)
template_lines <- gsub("summ3", summary_str("sampling_prop", n = 3), template_lines)
template_lines <- gsub("summ4", summary_str("epi_duration"), template_lines)
template_lines <- gsub("summ5", summary_str("final_prevalence", n = 1, as_e = TRUE), template_lines)
template_lines <- gsub("summ6", summary_str("final_cumulative_infections", n = 1, as_e = TRUE), template_lines)

## write the modified lines to a new file "prior-latex-table.tex"
writeLines(template_lines, "./src/prior-latex-table.tex")
