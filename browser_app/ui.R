shinyUI(fluidPage(pageWithSidebar(
    headerPanel("Mouse Kidney FibrOmics browser"),
    sidebarPanel(
        selectizeInput("gene", "Gene ID:", choices = "ENSMUSG00000000420_Galnt1"),
        selectizeInput("mirna", "miRNA ID:", choices = "mmu-miR-27b-5p"),
        p("Click on the text box and delete the current text. Start typing to search a gene or miRNA."),
        
        
        submitButton(text = "Update View", icon("refresh"))
    ),
    mainPanel(
        tabsetPanel(
            tabPanel("Plot", 
                     br(),
                     p(strong("Normalized abundance of mRNA, protein and miRNA:")),
                     p("Values are z-scores of the log2 normalized abundance for each molecule."),
                     plotOutput("distPlot"),
                     br(),
                     p(strong("Table values:")),
                     p("Normalized abundance (log2 counts) of molecules in each mouse model."),
                     tableOutput("fa"),
                     tableOutput("uuo"),
                     br(),
                     p(strong("miRNAs targetting this gene:")),
                     p("Predicted miRNA regulating this gene according targetscan DB."),
                     dataTableOutput("targets")),
            tabPanel("About", 
                     withTags({
                         div(class="header", checked=NA,
                             p("Multi Omics Analysis of Fibrotic kidneys in two Mouse Models"),
                             br(),
                             p("Mira Pavkovic, Lorena Pantano, Cory V. Gerlach, Sergine Brutus, Sarah A. Boswell, Robert A. Everley, Jagesh Shah, Shannan H Sui,  Vishal S. Vaidya")
                             )
                     }),

                     HTML("Kidney fibrosis represents an urgent unmet clinical need due to the lack of effective therapies and inadequate understanding of the molecular pathogenesis. 
We have generated a comprehensive and integrated multi-omics data set (proteomics, mRNA and small RNA transcriptomics) of fibrotic kidneys that is searchable through a user-friendly web application. Two commonly used mouse models were utilized: a reversible chemical-induced injury model (folic acid (FA) induced nephropathy) and an irreversible surgically-induced fibrosis model (unilateral ureteral obstruction (UUO)). mRNA and small RNA sequencing as well as 10-plex tandem mass tag (TMT) proteomics were performed with kidney samples from different time points over the course of fibrosis development. 
                          The bioinformatics workflow used to process, technically validate, and integrate the single data sets is described in the manuscript and raw data files can be found on GEO GSE118341 
                          In summary, we present temporal and integrated multi-omics data from fibrotic mouse kidneys that are accessible through an interrogation tool to provide a searchable transcriptome and proteome for kidney fibrosis researchers. 
                          ")
                     )
        )
        
    )
)))