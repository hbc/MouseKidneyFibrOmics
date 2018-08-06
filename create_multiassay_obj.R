library(tidyverse)
library(janitor)
library(MultiAssayExperiment)
library(isomiRs)
root = "."

# FA model
order_group=c("normal", "day1", "day2", "day3", "day7", "day14")
mrna_path = "tables/fa/results/counts/vst.csv.gz"
mrna_matrix = read_csv(mrna_path)[,c(1, 17:19, 2:4, 8:16, 5:7)] %>% 
    .[,c(1:8,10:19)] %>% 
    as.data.frame() %>% 
    column_to_rownames("rowname") %>% 
    as.matrix()
mrna_col = data.frame(row.names=colnames(mrna_matrix), samples=colnames(mrna_matrix)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day, levels=order_group))
rownames(mrna_col) = colnames(mrna_matrix)

prot_path = "tables/pfa/results/counts/fa_model_log2_counts.csv"
prot_matrix =  read_csv(prot_path) %>% 
    clean_names() %>% 
    set_names(gsub("control", "normal", colnames(.[]))) %>% 
    group_by(id) %>% 
    summarise_all(funs(sum)) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    column_to_rownames("id") %>% 
    as.matrix()
prot_col =  data.frame(row.names=colnames(prot_matrix), samples=colnames(prot_matrix)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day, levels=order_group))
rownames(prot_col) = colnames(prot_matrix)

mirna_path = "tables/sfa/results/counts/mirna_rlog.csv"
mirna_matrix =   read_csv(mirna_path)[,c(1, 6:8, 2, 9:10, 3:4, 11, 5, 12:13 )] %>% 
    clean_names() %>% 
    as.data.frame() %>% 
    column_to_rownames("mirna") %>% 
    as.matrix()
mirna_col =  data.frame(row.names=colnames(mirna_matrix), samples=colnames(mirna_matrix)) %>%
    mutate(day=gsub(".*_", "", samples)) %>%
    mutate(day=factor(day, levels=order_group)) %>% 
    .[,c("day"), drop = FALSE]
rownames(mirna_col) = colnames(mirna_matrix)

library(org.Mm.eg.db)
library(targetscan.Mm.eg.db)
clean_mir = DEGreport::degFilter(mirna_matrix, mirna_col, "day", min = 1, minreads = 7)
clean_mrna = DEGreport::degFilter(mrna_matrix, mrna_col, "day", min = 1, minreads = 7)
pairs = mirna2targetscan(mirna = rownames(clean_mir), species = "mmu", org = org.Mm.eg.db, keytype = "ENSEMBL")
targets = findTargets(SummarizedExperiment(assays = SimpleList(norm = clean_mir),
                                           colData = mirna_col,
                                           metadata = list(sign = rownames(clean_mir))),
                      SummarizedExperiment(assays = SimpleList(norm = clean_mrna),
                                           colData = mrna_col,
                                           metadata = list(sign = rownames(clean_mrna))),
                      pairs[,c("ENSEMBL", "mir")], "day")
targets_ready = reshape::melt(targets) %>%
    dplyr::filter(value != 0) %>%
    set_names(c("gene", "mir", "cor"))

fa_exp = ExperimentList = list(
    mirna = mirna_matrix, 
    gene = mrna_matrix,
    protein = prot_matrix
)

fa_map = listToMap(list(mirna = data.frame(primary = row.names(mirna_col), 
                                           colname = row.names(mirna_col),
                                           stringsAsFactors = FALSE),
                        gene = data.frame(primary = row.names(mrna_col), 
                                          colname = row.names(mrna_col),
                                          stringsAsFactors = FALSE),
                        protein = data.frame(primary = row.names(prot_col), 
                                          colname = row.names(prot_col),
                                          stringsAsFactors = FALSE)))
fa_model = MultiAssayExperiment(experiments = fa_exp,
                                sampleMap = fa_map,
                                colData = rbind(mirna_col, mrna_col, prot_col),
                                metadata = list(targets = targets_ready) 
)

# UUO model
order_group=c("normal", "day1", "day2", "day3", "day7", "day14")
mrna_path = "tables/uuo/results/counts/vst.csv.gz"
mrna_matrix = read_csv(mrna_path)[,c(1, 14:16,  6:9, 10:13, 2:5)] %>% 
    as.data.frame() %>% 
    column_to_rownames("rowname") %>% 
    as.matrix()
mrna_col = data.frame(row.names=colnames(mrna_matrix), samples=colnames(mrna_matrix)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day, levels=order_group))
rownames(mrna_col) = colnames(mrna_matrix)

prot_path = "tables/puuo/results/counts/uuo_model_log2_counts.csv"
prot_matrix =  read_csv(prot_path) %>% 
    clean_names() %>% 
    set_names(gsub("control", "normal", colnames(.[]))) %>% 
    group_by(id) %>% 
    summarise_all(funs(sum)) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    column_to_rownames("id") %>% 
    as.matrix()
prot_col =  data.frame(row.names=colnames(prot_matrix), samples=colnames(prot_matrix)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day, levels=order_group))
rownames(prot_col) = colnames(prot_matrix)

mirna_path = "tables/suuo/results/counts/mirna_rlog.csv"
mirna_matrix =   read_csv(mirna_path)[,c(1, 14:16, 6:9, 10:13, 2:5)] %>% 
    clean_names() %>% 
    as.data.frame() %>% 
    column_to_rownames("mirna") %>% 
    as.matrix()
mirna_col =  data.frame(row.names=colnames(mirna_matrix), samples=colnames(mirna_matrix)) %>%
    separate(samples, into = c("day"), extra = "drop", sep = "_") %>%
    mutate(day=factor(day, levels=order_group)) %>% 
    .[,c("day"), drop = FALSE]
rownames(mirna_col) = colnames(mirna_matrix)

library(org.Mm.eg.db)
library(targetscan.Mm.eg.db)
clean_mir = DEGreport::degFilter(mirna_matrix, mirna_col, "day", min = 1, minreads = 7)
clean_mrna = DEGreport::degFilter(mrna_matrix, mrna_col, "day", min = 1, minreads = 7)
pairs = mirna2targetscan(mirna = rownames(clean_mir), species = "mmu", org = org.Mm.eg.db, keytype = "ENSEMBL")
targets = findTargets(SummarizedExperiment(assays = SimpleList(norm = clean_mir),
                                           colData = mirna_col,
                                           metadata = list(sign = rownames(clean_mir))),
                      SummarizedExperiment(assays = SimpleList(norm = clean_mrna),
                                           colData = mrna_col,
                                           metadata = list(sign = rownames(clean_mrna))),
                      pairs[,c("ENSEMBL", "mir")], "day")
targets_ready = reshape::melt(targets) %>%
    dplyr::filter(value != 0) %>%
    set_names(c("gene", "mir", "cor"))


uuo_exp = ExperimentList = list(
    mirna = mirna_matrix, 
    gene = mrna_matrix,
    protein = prot_matrix
)

uuo_map = listToMap(list(mirna = data.frame(primary = row.names(mirna_col), 
                                           colname = row.names(mirna_col),
                                           stringsAsFactors = FALSE),
                        gene = data.frame(primary = row.names(mrna_col), 
                                          colname = row.names(mrna_col),
                                          stringsAsFactors = FALSE),
                        protein = data.frame(primary = row.names(prot_col), 
                                             colname = row.names(prot_col),
                                             stringsAsFactors = FALSE)))
uuo_model = MultiAssayExperiment(experiments = uuo_exp,
                                sampleMap = uuo_map,
                                colData = rbind(mirna_col, mrna_col, prot_col),
                                metadata = list(targets = targets_ready)
)


obj = list(fa=fa_model, uuo=uuo_model)

save(obj, file=file.path("../data", "me.rda"))
