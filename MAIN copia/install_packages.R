list_of_packages <- c("readr", "ParallelPC","parallel","pcalg")

new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]

if(length(new_packages)){
  install.packages(new_packages)
} 