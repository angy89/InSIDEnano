phenotypic_network_UI = function(input,output){
  output$items_net_page1 = renderUI({
    sidebarPanel(
      fluidRow(
        column(12, checkboxGroupInput("checkGroup", 
                                      label = "Object Network", 
                                      choices = list("Nanomaterials" = "nano", 
                                                     "Drugs" = "drugs", "Chemical" = "chemical",
                                                     "Disease" = "disease"),selected = "nano"))
      ),
      fluidRow(column(12,selectInput("Edges",
                                     label="Select the type of edges",
                                     choices = list("Positive" = "P",
                                                    "Negative" = "N",
                                                    "Both" = "B"),selected="B"))
      ),
      fluidRow(column(12,uiOutput("input_dis"))),
      
      fluidRow(
        column(12, selectInput("ATCGroup", 
                               label ="Drugs ATC code", multiple = TRUE,
                               choices = ATC_choice_list,selected = "A"))
      ),
      
      fluidRow(
        column(12, selectInput("ChemicalGroup", 
                               label ="Chemical Classes", multiple = TRUE,
                               choices = chemical_choice_list,
                               selected = "Amino acids"))
      ),
      fluidRow(
        column(9, sliderInput("slider1", label = " Connection Strength",
                              min = 1, max = 99, value = 99,step=1))
      ),
      fluidRow(
        column(9, sliderInput("repulserad", label = "Repulseration Strenght",
                              min = 100, max = 1000, value = 100,step=1))
      ),
      fluidRow(
        column(9, sliderInput("length", label = "Edge length",
                              min = 2, max = 10, value = 2,step=1))
      )
    ) #end Sidebar
  })
  
  output$items_net_page2 = renderUI({
    
    mainPanel(
      fluidRow(
        column(9,wellPanel(forceNetworkOutput("itemNetwork")))
      ),
      fluidRow(
        column(9,wellPanel(plotOutput("pie_chart")))
      )
    )
  })
}