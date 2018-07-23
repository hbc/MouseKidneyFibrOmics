library(bcbioRNASeq)
library(tidyverse)

load("../data/uuo.rda")
counts = counts(uuo, "tpm")
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, ".tpm")))
}

library(bcbioSmallRna)
load("../data/suuo.rda")
counts = mirna(suuo, "rlog")
for(c in colnames(counts)){
    counts[,c, drop = FALSE] %>% 
        as.data.frame %>% 
        rownames_to_column("gene") %>% 
        write_csv(file.path("geo", paste0(c, ".mirna")))
}
