# https://www.dropbox.com/sh/bmxjb06oi6eei3c/AABbhmXB7Qw-CTwAdqzRfTh9a?dl=1
# library(r2dropSmart)
# token = readRDS("~/.droptoken.rds")
# sync("browser_app", "fmm", dry = F, token = token)

library(ggplot2)
library(MultiAssayExperiment)
library(tibble)
library(dplyr)
library(tidyr)
options(shiny.maxRequestSize = -1)
load("se.rda")
assign('obj', obj, envir=.GlobalEnv)
df <- NULL
fa <- obj[["fa"]]
uuo <- obj[["uuo"]]
possible_genes <- c(
    row.names(experiments(fa)[["gene"]]),
    row.names(experiments(fa)[["protein"]]),
    row.names(experiments(uuo)[["gene"]]),
    row.names(experiments(uuo)[["protein"]])
) %>% as.character() %>% unique()

possible_mirnas <- c(
    row.names(experiments(fa)[["mirna"]]),
    row.names(experiments(uuo)[["mirna"]])
) %>% as.character() %>% unique()
