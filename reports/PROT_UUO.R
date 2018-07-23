library(tidyverse)
library(rio)
library(janitor)

# root_path = "~/orch/scratch/vishal_mirna_kidney/publish/FA-model"
result_files = file.path("..", "tables/puuo/results/counts")
dir.create(result_files, showWarnings = F, recursive = T)

protein_file = file.path("~/orch/hbc/PIs/vishal_vaidya/vishal_mirna_kidney/scientific_data/uuo/protein/UUO_Proteomics_RawData.csv")

dat = read_csv(protein_file, skip = 1) %>% 
    clean_names() %>%
    as_tibble() %>% 
    separate(protein_id, c("sp", "uniprot_id", "type"), sep = "[::|::]", extra = "merge")


library(EnsDb.Mmusculus.v79)
mapping = select(EnsDb.Mmusculus.v79, unique(dat$uniprot_id), columns = "GENEID", keytype = "UNIPROTID") %>% 
    clean_names()

dat_ensg = inner_join(dat, mapping, by = c("uniprot_id" = "uniprotid"))

counts = dat_ensg %>% dplyr::select(normal_1:day14_2) %>% as.matrix()
row.names(counts) = dat_ensg$geneid
write_csv(counts %>% 
              as.data.frame() %>%
              rownames_to_column("id"),
          file.path(result_files, "uuo_model_counts.csv"))

dge = edgeR::DGEList(counts)
dge = edgeR::calcNormFactors(dge, method="TMM")
norm_counts = edgeR::cpm(dge, log=TRUE)

write_csv(norm_counts %>% 
              as.data.frame() %>%
              rownames_to_column("id"),
          file.path(result_files, "uuo_model_log2_counts.csv"))

prot_col =  data.frame(row.names=colnames(norm_counts), samples=colnames(norm_counts)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day))
rownames(prot_col) = colnames(norm_counts)

DEGreport::degPCA(norm_counts, prot_col, "day")
