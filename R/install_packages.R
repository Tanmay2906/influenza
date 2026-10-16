# Run once to install the R packages this project needs.
pkgs <- c("dplyr", "tidyr", "readr", "ggplot2", "stringr", "ape")
missing <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
