library(bcbioRNASeq)
library(tidyverse)
load("../data/fa.rda")
metrics(fa) %>% select(totalReads, rrnaRate, x5x3Bias, mappedReads, exonicRate)

library(bcbioSmallRna)
load("../data/sfa.rda")
metrics(sfa) %>% select(sample,reads_before_trimming)

