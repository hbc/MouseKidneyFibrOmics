library(tidyverse)
library(janitor)
library(MultiAssayExperiment)
library(DEGreport)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(cowplot)
library(ggrepel)
theme_set(
    theme_light(base_size = 12L))
theme_update(
    legend.justification = "center",
    legend.position = "bottom")


uuo_cols = RColorBrewer::brewer.pal(8, "Dark2")[c(1,4,5,6)]
fa_cols = RColorBrewer::brewer.pal(8, "Dark2")[c(1,2,3,5,6)]
fa_mirna_cols = RColorBrewer::brewer.pal(8, "Dark2")[c(1,5,6)]

# > ·         For fibrosis aSMA=acta2; col1a1, fn1
# > ·         For kidney injury: ngal=lcn2, kim-1, clusterin
# > ·         Housekeeping genes: Ppia, Gapdh, Pgk1, Gusb, Polr2a, Actb, Tbp, Hprt, B2M, Rp2, Rpl32 (here, I chose some more, since some of them might not be good controls for UUO as a heavy fibrosis model) > ·         miRNAs: mir-21, miR-155, miR-192

load("data/me.rda")
cold = colData(obj$uuo) %>% as.data.frame() %>% 
    rownames_to_column("sample")
coldfa = colData(obj$fa) %>% as.data.frame() %>% 
    rownames_to_column("sample")

plot_grid(
    degPCA(experiments(obj$uuo)[["mirna"]],
           colData(obj$uuo)[colnames(experiments(obj$uuo)[["mirna"]]),,drop=F],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = uuo_cols) +
        ggtitle("A) UUO miRNA"),
    degPCA(experiments(obj$uuo)[["gene"]],
           colData(obj$uuo)[colnames(experiments(obj$uuo)[["gene"]]),,drop=F],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = uuo_cols) +
        ggtitle("B) UUO mRNA"),
    degPCA(experiments(obj$uuo)[["protein"]],
           colData(obj$uuo)[colnames(experiments(obj$uuo)[["protein"]]),,drop=F],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = uuo_cols) +
        ggtitle("C) UUO protein"),
    degPCA(experiments(obj$fa)[["protein"]],
           colData(obj$fa)[colnames(experiments(obj$fa)[["protein"]]),,drop=F],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = fa_cols) +
        ggtitle("D) FA protein")
) + ggsave("figures/uuo-pca.pdf", width = 9, height = 9)


fa_gene = experiments(obj$fa)[["gene"]]
fa_cd = colData(obj$fa)[colnames(experiments(obj$fa)[["gene"]]),,drop=F]
keep_gene = fa_cd$day != "day3"
fa_mirna = experiments(obj$fa)[["mirna"]]
fa_cd_mirna = colData(obj$fa)[colnames(experiments(obj$fa)[["mirna"]]),,drop=F]
keep_mirna = fa_cd_mirna$day != "day3"
plot_grid(
    degPCA(fa_gene[,keep_gene],
           fa_cd[keep_gene,,drop=FALSE],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = fa_cols) +
        ggtitle("B) FA mRNA"),
    degPCA(fa_mirna[, keep_mirna],
           fa_cd_mirna[keep_mirna,,drop=FALSE],
           "day") +
        geom_text_repel(aes(label=day)) + 
        scale_color_manual("", guide = FALSE, values = fa_mirna_cols) +
        ggtitle("A) FA miRNA")
) + ggsave("figures/fa_pca.pdf", width = 9, height = 5)


fibrosis = select(EnsDb.Mmusculus.v79, c("Acta2", "Col1a1", "Fn1"), columns = "GENEID", keytype = "SYMBOL")
injury = select(EnsDb.Mmusculus.v79, c("Lcn2", "Havcr1", "Clu"), columns = "GENEID", keytype = "SYMBOL")
hk = select(EnsDb.Mmusculus.v79, c("Ppia", "Gapdh", "Pgk1", "Gusb", "Polr2a", "Actb", "Tbp", "Hprt", "B2m"), columns = "GENEID", keytype = "SYMBOL")

# miR-155, miR-192

## Table with markers
df = bind_rows(
    experiments(obj$uuo)[["protein"]][fibrosis$GENEID,] %>%
        reshape::melt() %>%
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(fibrosis, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "fibrosis", data = "UUO protein"),
    experiments(obj$fa)[["protein"]][fibrosis$GENEID,] %>%
        reshape::melt() %>%
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(fibrosis, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "fibrosis", data = "FA protein"),
    experiments(obj$uuo)[["gene"]][fibrosis$GENEID,] %>%
        reshape::melt() %>% 
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(fibrosis, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "fibrosis", data="UUO mRNA"),
    experiments(obj$uuo)[["gene"]][injury$GENEID,] %>%
        reshape::melt() %>% 
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(injury, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "injury", data="UUO mRNA"),
    experiments(obj$uuo)[["protein"]][injury$GENEID[c(1,3)],] %>%
        reshape::melt() %>%  
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(injury, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "injury", data = "UUO protein"),
    experiments(obj$fa)[["protein"]][injury$GENEID[c(1,3)],] %>%
        reshape::melt() %>%
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(injury, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "injury", data = "FA protein"),
    experiments(obj$uuo)[["gene"]][hk$GENEID,] %>%
        reshape::melt() %>% 
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(hk, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "house-keeping", data="UUO mRNA"),
    experiments(obj$uuo)[["protein"]][hk$GENEID,] %>%
        reshape::melt() %>%  
        left_join(cold, by = c("X2" = "sample")) %>% 
        left_join(hk, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "house-keeping", data = "UUO protein"),
    experiments(obj$fa)[["protein"]][hk$GENEID,] %>%
        reshape::melt() %>%
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(hk, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "house-keeping", data = "FA protein")
) %>% 
    mutate(group = paste(SYMBOL, data))

## Plot
theme_set(
    theme_light(base_size = 9L))
plot_grid(
    ggplot(dplyr::filter(df, marker == "fibrosis"),
           aes(x = day, y = value,
               color = SYMBOL, group = group,
               shape = data, linetype = data)) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("A) fibrosis markers"),
    ggplot(dplyr::filter(df, marker == "injury"),
           aes(x = day, y = value,
               color = SYMBOL, group = group,
               shape = data, linetype = data)) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("B) injury markers"),
    ggplot(dplyr::filter(df, marker == "house-keeping"),
           aes(x = day, y = value,
               color = SYMBOL, group = group,
               shape = data, linetype = data)) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("C) house-keeping markers"),
    experiments(obj$uuo)[["mirna"]][c("mmu-miR-21a-5p",
                                      "mmu-miR-192-5p"),,drop=F] %>%
        reshape::melt() %>% 
        left_join(cold, by = c("X2" = "sample")) %>%
        ggplot(aes(x = day, y = value, color= X1, group = X1)) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        geom_smooth(se = FALSE) +
        ylab("log2 normalized counts") +
        ggtitle("D) miRNA fibrosis markers")
    
) + ggsave("figures/markers.pdf", width = 11, height = 9)

# Same figures for FA publised data
df = bind_rows(
    fa_gene[fibrosis$GENEID,keep_gene] %>%
        reshape::melt() %>% 
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(fibrosis, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "fibrosis", data="FA mRNA"),
    fa_gene[injury$GENEID,keep_gene] %>%
        reshape::melt() %>% 
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(injury, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "injury", data="FA mRNA"),
    fa_gene[hk$GENEID,keep_gene] %>%
        reshape::melt() %>% 
        left_join(coldfa, by = c("X2" = "sample")) %>% 
        left_join(hk, by = c("X1" = "GENEID")) %>% 
        mutate(marker = "house-keeping", data="FA mRNA")
)
plot_grid(
    ggplot(dplyr::filter(df, marker == "fibrosis"),
           aes(x = day, y = value,
               color = SYMBOL, group = SYMBOL,
               )) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("A) fibrosis markers"),
    ggplot(dplyr::filter(df, marker == "injury"),
           aes(x = day, y = value,
               color = SYMBOL, group = SYMBOL,
               )) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("B) injury markers"),
    ggplot(dplyr::filter(df, marker == "house-keeping"),
           aes(x = day, y = value,
               color = SYMBOL, group = SYMBOL,
               )) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        scale_linetype_manual(guide = FALSE, values = c(1,3,6)) +
        geom_smooth(se = FALSE)  +
        ylab("log2 normalized counts") +
        ggtitle("C) house-keeping markers"),
    fa_mirna[c("mmu-miR-21a-5p",
                                      "mmu-miR-192-5p"),keep_mirna,drop=F] %>%
        reshape::melt() %>% 
        left_join(coldfa, by = c("X2" = "sample")) %>%
        ggplot(aes(x = day, y = value, color= X1, group = X1)) +
        geom_point() +
        scale_color_brewer("", palette = "Set2") +
        geom_smooth(se = FALSE) +
        ylab("log2 normalized counts") +
        ggtitle("D) miRNA fibrosis markers")
    
) + ggsave("figures/fa_markers.pdf", width = 11, height = 9)


df %>% dplyr::filter(marker == "house-keeping") %>% group_by(data, SYMBOL) %>% summarise(cv = sd(value)/mean(value)*100) %>%  group_by(data) %>% summarise(cv_mean = mean(cv)) %>% rio::export("figures/cv-hk-genes.csv")
