
shinyServer(function(input, output, session) {
    updateSelectizeInput(session, "gene",
                         choices  = possible_genes,
                         selected = "ENSMUSG00000000420_Galnt1",
                         server = TRUE)
    updateSelectizeInput(session, "mirna",
                         choices  = possible_mirnas,
                         selected = "mmu-miR-27b-5p",
                         server = TRUE)
    
    datasetInput <- reactive({
        gene <- input$gene
        mirna <- input$mirna

        cfa <- as.data.frame(colData(fa)) %>%
            rownames_to_column("sample")
        cuuo <- as.data.frame(colData(uuo)) %>%
            rownames_to_column("sample")
        fa_g = fa_p = fa_m = data.frame()
        uuo_g = uuo_p = uuo_m = data.frame()

        if (gene %in% rownames(experiments(fa)[["gene"]])){
            fa_g = experiments(fa)[["gene"]][gene,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            fa_g = gather(fa_g, "sample", "expression", -feature) %>%
                left_join(cfa, by = "sample") %>%
                mutate(molecule = "mRNA", model = "FA")
        }
        if (gene %in% rownames(experiments(fa)[["protein"]])){
            fa_p = experiments(fa)[["protein"]][gene,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            fa_p = gather(fa_p, "sample", "expression", -feature) %>%
                left_join(cfa, by = "sample") %>%
                mutate(molecule = "protein", model = "FA")
        }
        if (mirna %in% rownames(experiments(fa)[["mirna"]])){
            fa_m = experiments(fa)[["mirna"]][mirna,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            fa_m = gather(fa_m, "sample", "expression", -feature) %>%
                left_join(cfa, by = "sample") %>%
                mutate(molecule = "miRNA", model = "FA")
        }


        if (gene %in% rownames(experiments(uuo)[["gene"]])){
            uuo_g = experiments(uuo)[["gene"]][gene,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            uuo_g = gather(uuo_g, "sample", "expression", -feature) %>%
                left_join(cuuo, by = "sample") %>%
                mutate(molecule = "mRNA", model = "UUO")
        }
        if (gene %in% rownames(experiments(uuo)[["protein"]])){
            uuo_p = experiments(uuo)[["protein"]][gene,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            uuo_p = gather(uuo_p, "sample", "expression", -feature) %>%
                left_join(cuuo, by = "sample") %>%
                mutate(molecule = "protein", model = "UUO")
        }
        if (mirna %in% rownames(experiments(uuo)[["mirna"]])){
            uuo_m = experiments(uuo)[["mirna"]][mirna,,drop = F] %>%
                as.data.frame() %>%
                rownames_to_column("feature")
            uuo_m = gather(uuo_m, "sample", "expression", -feature) %>%
                left_join(cuuo, by = "sample") %>%
                mutate(molecule = "miRNA", model = "UUO")
        }


        fa_df = bind_rows(fa_g, fa_p, fa_m) %>%
            group_by(feature, day, model, molecule) %>%
            summarise(expression = median(expression)) %>%
            ungroup() %>%
            group_by(feature, model, molecule) %>%
            mutate(zscore = scale(expression)) %>%
            ungroup()

        uuo_df = bind_rows(uuo_g, uuo_p, uuo_m)%>%
            group_by(feature, day, model, molecule) %>%
            summarise(expression = median(expression)) %>%
            ungroup()%>%
            group_by(feature, model, molecule) %>%
            mutate(zscore = scale(expression)) %>%
            ungroup()

        targets = bind_rows(metadata(fa)[["targets"]] %>%
            dplyr::rename(feature = gene) %>%
            dplyr::filter(feature == gene),
            metadata(uuo)[["targets"]] %>%
                dplyr::rename(feature = gene) %>%
                dplyr::filter(feature == gene)
        ) %>% distinct()

        list(fa=fa_df, uuo=uuo_df, targets = targets)
    })
    output$targets <- renderDataTable({
        df = datasetInput()
        df[["targets"]]
    })
    output$fa <- renderTable({
        df = datasetInput()
        dplyr::select(df[["fa"]], -zscore) %>%
            spread("day", "expression") %>% as.data.frame()
    })
    output$uuo <- renderTable({
        df = datasetInput()
        dplyr::select(df[["uuo"]], -zscore) %>%
            spread("day", "expression") %>% as.data.frame()
    })
    output$distPlot <- renderPlot({
        df = datasetInput()
        bind_rows(df[["fa"]], df[["uuo"]]) %>%
            ggplot(aes(day, zscore, color = molecule, group = molecule)) +
            geom_point(size = 3) +
            geom_line(size = 2) +
            # geom_smooth(method = "lm",formula = y~poly(x,3),
                                       # alpha=0.2, se = FALSE) +
            theme_bw(base_size = 16) +
            facet_wrap(~model, nrow = 2)

    })
})

# 
# Replace normal with day0
# Replace transcript with mRNA
# Maybe put ensemble number in brackets behind the gene name
# 
