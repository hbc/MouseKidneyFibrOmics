shinyUI(fluidPage(pageWithSidebar(
    headerPanel("MKF browser"),
    sidebarPanel(
        selectizeInput("gene", "Gene ID:", choices = "ENSMUSG00000026193_Fn1"),
        textInput("mirna", "miRNA ID:", "mmu-miR-192-5p"),
        
        submitButton("do","Update View")
    ),
    mainPanel(
        plotOutput("distPlot"),
        p(strong("Table values:")),
        tableOutput("fa"),
        tableOutput("uuo"),
        p(strong("miRNAs targetting this gene:")),
        dataTableOutput("targets")
    )
)))