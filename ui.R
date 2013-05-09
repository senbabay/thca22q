library(shiny)
# Define UI for slider demo application
shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("THCA chr 22q loss"),
  
  # Sidebar with sliders that demonstrate various available options
  sidebarPanel(
    
    # Choose mutation context
    selectInput("context", "Choose mutation context:", 
                choices = c("BRAF", "RAS")),
    
    # Specification of the range for log2 copy number values for HETLOSS 
    sliderInput("int1", "Choose interval for heterozygous loss (log2r values):",
                min = -1, max = 0, value = c(-0.95,-0.2),step=0.01),
    
    br(),
    
    textInput("genename", "Show plots for gene:", "CHEK2"),
    
    br(),
    
    # Choose scoring scheme
    selectInput("scheme", "Choose ranking scheme for genes:", 
                choices = c("Mann-Whitney test between HETLOSS and DIPLOID expression values",
                            "Expression vs copy number distance correlation",
                            "Expression unimodality (Hartigan's dip test)",
                            "Number of mutations",
                            "Expression vs methylation correlation")),
    
    numericInput("num2show", "Number of genes to view:", 10),
    
    h5("Gene ranklist"),
    tableOutput("view")
    
    
  ),
  
 
  # Show a table summarizing the values entered
  mainPanel(

    h4(textOutput("currentGene")),
    plotOutput("row1plots"),
    plotOutput("row2plots"),
    br(),
    tableOutput("geneStats")
    #verbatimTextOutput("summary")
  )
))
