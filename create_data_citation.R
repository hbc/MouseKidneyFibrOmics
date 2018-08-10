library(MultiAssayExperiment)
library(tidyverse)
load("data/me.rda")

format = . %>% as.data.frame %>% rownames_to_column("id")

write_csv(experiments(obj$fa)[["gene"]] %>% format, "data_citation/fa_mrna.csv.gz")
write_csv(experiments(obj$fa)[["protein"]] %>% format, "data_citation/fa_protein.csv.gz")
write_csv(experiments(obj$fa)[["mirna"]] %>% format, "data_citation/fa_mirna.csv.gz")

write_csv(experiments(obj$uuo)[["gene"]] %>% format, "data_citation/uuo_mrna.csv.gz")
write_csv(experiments(obj$uuo)[["protein"]] %>% format, "data_citation/uuo_protein.csv.gz")
write_csv(experiments(obj$uuo)[["mirna"]] %>% format, "data_citation/uuo_mirna.csv.gz")
