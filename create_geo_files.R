library(bcbioRNASeq)
library(tidyverse)

load("data/uuo.rda")
counts = counts(uuo, "tpm")
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, ".tpm")))
}

library(bcbioSmallRna)
load("data/suuo.rda")
counts = mirna(suuo, "rlog")
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, ".mirna")))
}

load("data/me.rda")
counts = experiments(obj[["uuo"]])[["protein"]]
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, "_uuo.tsv")))
}

counts = experiments(obj[["fa"]])[["protein"]]
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, "_fa.tsv")))
}

lapply(list.files("geo", full.names = T, pattern = "tsv"), function(fn){
    data.frame(sample = tools::file_path_sans_ext(basename(fn)),
               md5 = tools::md5sum(fn),
               fn = basename(fn),
               stringsAsFactors = T)
}) %>% bind_rows()


lapply(list.files("geo", full.names = T, pattern = "tpm"), function(fn){
    data.frame(sample = tools::file_path_sans_ext(basename(fn)),
               md5 = tools::md5sum(fn),
               fn = basename(fn),
               stringsAsFactors = T)
}) %>% bind_rows()

lapply(list.files("geo", full.names = T, pattern = "mirna"), function(fn){
    data.frame(sample = tools::file_path_sans_ext(basename(fn)),
               md5 = tools::md5sum(fn),
               fn = basename(fn),
               stringsAsFactors = T)
}) %>% bind_rows()

lapply(list.files("geo", full.names = T, pattern = "gtf"), function(fn){
    data.frame(sample = tools::file_path_sans_ext(basename(fn)),
               md5 = tools::md5sum(fn),
               fn = basename(fn),
               stringsAsFactors = T)
}) %>% bind_rows()
