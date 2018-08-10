library(bcbioRNASeq)
library(tidyverse)
load("data/fa.rda")
metrics(fa) %>% select(totalReads, rrnaRate, x5x3Bias, mappedReads, exonicRate)

library(bcbioSmallRna)
load("data/sfa.rda")
metrics(sfa) %>% select(sample,reads_before_trimming, read_with_adapter, sequence_length) %>% 
    mutate(pct_w_adapter = read_with_adapter/reads_before_trimming*100)


load("data/suuo.rda")
metrics(suuo) %>% select(sample,reads_before_trimming, read_with_adapter, sequence_length) %>% 
    mutate(pct_w_adapter = read_with_adapter/reads_before_trimming*100)

