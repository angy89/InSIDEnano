gene_network_UI = function(input,output){
  output$gene_net_page1 = renderUI({
    sidebarPanel(
      fluidRow(column(12,uiOutput("Patway"))),
      fluidRow(
        column(9, sliderInput("gene_repulserad", label = "Repulseration Strenght",
                              min = 100, max = 1000, value = 100,step=1))
      ),
      fluidRow(
        column(9, sliderInput("gene_length", label = "Edge length",
                              min = 2, max = 10, value = 2,step=1))
      )
    )
  })
  
  output$gene_net_page2 = renderUI({
    mainPanel(
      fluidRow(column(12, wellPanel(forceNetworkOutput("geneNetwork"))))   
    )
  })

  
#   output$gene_net_page3 = renderUI({
#     sidebarPanel(
#       fluidRow(
#         column(9, 
#                selectizeInput(
#                  'gene_nano_query', label = "Select the nanomaterial to start the analysis", choices = c(nano), multiple = FALSE,selected="MWCNT",
#                  options = list(create = TRUE)
#                )       
#         )
#       ),
#       fluidRow(
#         column(9, 
#                selectizeInput(
#                  'gene_disease_query', label = "Select the disease to start the analysis", choices = c(disease), multiple = FALSE,selected="Asthma",
#                  options = list(create = TRUE)
#                )       
#         )
#       ),
#       fluidRow(
#         column(9, sliderInput("th_gene_query", label = "Percentage of drugs",
#                               min = 1, max = 100, value = 99,step=1))
#       )
#    
#    )
#    })
#    
# output$gene_net_page4 = renderUI({
#      mainPanel(
#        fluidRow(column(12, wellPanel(plotOutput("geneNetwork_net3"))))   
#      )
#      
#    })
#   
 
}
  