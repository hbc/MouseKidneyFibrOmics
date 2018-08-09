library(EnsDb.Mmusculus.v79)
library(MultiAssayExperiment)
library(tidyverse)
load("data/me.rda")

fa = obj[["fa"]]
uuo = obj[["uuo"]]

fa_genes = unique(row.names(experiments(fa)[["gene"]]), 
                  row.names(experiments(fa)[["protein"]]))
map_key = select(EnsDb.Mmusculus.v79, fa_genes,
       columns = "SYMBOL", keytype = "GENEID") %>% 
    mutate(id = paste(GENEID, SYMBOL, sep = "_"))

idx = match(row.names(experiments(fa)[["gene"]]), map_key[["GENEID"]])
row.names(experiments(fa)[["gene"]]) = map_key[idx, "id"]
idx = match(row.names(experiments(fa)[["protein"]]), map_key[["GENEID"]])
row.names(experiments(fa)[["protein"]]) = map_key[idx, "id"]
idx = match(metadata(fa)[["targets"]][["gene"]], map_key[["GENEID"]])
metadata(fa)[["targets"]][["gene"]] = map_key[idx, "id"]

uuo_genes = unique(row.names(experiments(uuo)[["gene"]]), 
                  row.names(experiments(uuo)[["protein"]]))
map_key = select(EnsDb.Mmusculus.v79, fa_genes,
                 columns = "SYMBOL", keytype = "GENEID") %>% 
    mutate(id = paste(GENEID, SYMBOL, sep = "_"))

idx = match(row.names(experiments(uuo)[["gene"]]), map_key[["GENEID"]])
row.names(experiments(uuo)[["gene"]]) = map_key[idx, "id"]
idx = match(row.names(experiments(uuo)[["protein"]]), map_key[["GENEID"]])
row.names(experiments(uuo)[["protein"]]) = map_key[idx, "id"]
idx = match(metadata(uuo)[["targets"]][["gene"]], map_key[["GENEID"]])
metadata(uuo)[["targets"]][["gene"]] = map_key[idx, "id"]

obj = list(fa=fa, uuo=uuo)
save(obj, file = "browser_app/se.rda")
