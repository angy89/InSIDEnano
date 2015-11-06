library(shiny)
shinyUI(navbarPage("",
                   tabPanel("INSIdE nano", icon=icon("home",lib = "glyphicon"),
                            fluidRow(
                              column(2,img(src='logo.png', align = "left",width = 150)),
                              column(7,titlePanel("\n Integrated Network of Systems Biology Effects of Nanos"))
                            ),
                            
                            fluidRow( 
                              column(4, 
                                wellPanel(
                                      fluidRow(
                                          column(7, textInput("username", "Username",value="", width = NULL))
                                      ),
                                      fluidRow(
                                        column(7, passwordInput("password", "Password",value="", width = NULL))
                                      ),
                                      fluidRow(
                                        column(7,actionButton(inputId = "login",label  = "Login",icon("circle-arrow-right", lib = "glyphicon"))),
                                          singleton(
                                            tags$head(tags$script(src = "message-handler.js"))
                                          )
                                      )
                                ) #end login panel    
                              ),
                              column(8, 
                                       tags$head(tags$script('
                                                             Shiny.addCustomMessageHandler("myCallbackHandler",
                                                             function(typeMessage) {console.log(typeMessage)
                                                             if(typeMessage == 1){
                                                             console.log("got here");
                                                             $("a:contains(Free Query)").click();
                                                             }
                                                             if(typeMessage == 5){
                                                             console.log("got here");
                                                             $("a:contains(Conditional Query)").click();
                                                             }
                                                             if(typeMessage == 2){
                                                             $("a:contains(Phenotypic Network)").click();
                                                             }
                                                             if(typeMessage == 3){
                                                             $("a:contains(Gene Network)").click();
                                                             }
                                                             if(typeMessage == 4){
                                                             $("a:contains(Demo)").click();
                                                             }
                                                             });
                                                             ')),
                                     wellPanel(
                                       fluidRow(column(6, div(style="display:inline-block",actionButton("action1", label = "Free Query",width = 400, icon=icon("search",lib = "glyphicon")), style="align:center")),
                                                column(6, div(style="display:inline-block",actionButton("action4", label = "Conditional Query",width = 400, icon=icon("search",lib = "glyphicon")), style="align:center")))
                                     ),
                                     wellPanel(
                                        fluidRow(column(6, div(style="display:inline-block",actionButton("action2", label = "Browse Phenotypic Network",width = 400,icon=icon("search",lib = "glyphicon")), style="align:center")),
                                                column(6, div(style="display:inline-block",actionButton("action3", label = "Browse Gene Network",width = 400,icon=icon("search",lib = "glyphicon")), style="align:center")))
                                     )  
                                     
                              ) #end buttonAction column
                             
                            ),
                            fluidRow(
                              column(12,
                                     wellPanel(
                                       fluidRow(
                                         p("INSIdE nano is a graphical tool  
                                           that  highlights connections between phenotypic entities including :"),
                                         #uiOutput("myList"),
                                         HTML("<ul>
                                              <li>Nanomaterials exposures</li>
                                              <li>Drugs treatments</li>
                                              <li>Chemicals exposures</li>
                                              <li>Diseases</li>
                                              </ul>"),
                                         p("Based on their effects on the genes.")
                                         ),
                                       fluidRow(
                                         actionButton(inputId = "Demo",label  = "Read More",icon("circle-arrow-right", lib = "glyphicon"))
                                       )
                                         ) #end description wellPanel
                                       ) 
                              
                            )
                                      
                            
                   ),

                              tabPanel("Free Query",
                                      fluidRow(
                                        column(8,
                                               wellPanel(
                                                 fluidRow(
                                                   column(6, uiOutput("disease")),
                                                   column(6, uiOutput("th_slider"))
                                                 ),
                                                 fluidRow(
                                                   column(4, uiOutput("Go")),
                                                   column(4, uiOutput("refresh_free")),
                                                   column(4, uiOutput("ok"))
                                                   
                                                 )
                                               )
                                        ),
                                        column(4,
                                          wellPanel(
                                            HTML("</br>The Free Analysis allows to extrapolate subnetworks of nodes connected to one or more query items.</br>"),
                                            HTML("The results of the query depends on the strenght of similarity. </br>
                                                 In case of multiple query items, the final subnetwork will be the union of all the subtenetwork related to each item.
                                                 </br></br></br>"),
                                            htmlOutput("checkLOGIN_free")
                                          ))
                                      ),
                                      
                                       wellPanel(
                                         fluidRow(column(8,uiOutput("plotTripel_total"))),
                                         fluidRow(uiOutput("NodesOfInterest_totale")),
                                         fluidRow( plotOutput('ggplot_totale',height=900))
                                       )     
                                       
                                       
                              ),
                              tabPanel("Conditional Query",
                                       fluidRow(
                                         column(8,
                                                wellPanel(
                                                  fluidRow(
                                                    column(3, uiOutput("nano_input")),
                                                    column(3, uiOutput("drug_input")),
                                                    column(3, uiOutput("disease_input")),
                                                    column(3, uiOutput("chemical_input"))
                                                  )
                                                ),
                                                wellPanel(
                                                  fluidRow(
                                                    column(3, uiOutput("th_slider2")),
                                                    column(3, uiOutput("percentuale_somma")),
                                                    column(3, uiOutput("nroCliques")),
                                                    column(3, uiOutput("clique_type"))
                                                  
                                                  )
                                                ),
                                                wellPanel(
                                                  fluidRow(
                                                    column(2,uiOutput("Go2")),
                                                    column(2,uiOutput("Refresh"))
                                                    
                                                  )
                                                )
                                         ),
                                         column(4,
                                                wellPanel(
                                                          fluidRow(htmlOutput("info2_1"))
                                                ),
                                                wellPanel(
                                                  htmlOutput("extimatedTime"),
                                                  HTML("<br/>")
                                                )
                                         )
                                       ),
                                       fluidRow(
                                         wellPanel(
                                           fluidRow(
                                             
                                             HTML("</br>The Conditional Analysis allows to extrapolate subnetworks of nodes connected to multiple query items and to analyze the pattern of items (cliques) in the subnetwork.</br>"),
                                             HTML("The results of the query depends on three paramenters: </br>
                                                  <ul>
                                                  <li>The Strength of similarity between query items and their neighbors</li>
                                                  <li>How many query items are connected at the same time with other selected nodes</li>
                                                  <li>How many query items are expected to be in the same cliques</li>
                                                  </ul>
                                                   
                                                  </br>"),
                                             htmlOutput("checkLOGIN")
                                             
                                             )
                                         )
                                       )
                              ),
                             
         navbarMenu("Query Results",
                    tabPanel("Items Subnetwork", 
                             sidebarPanel(
                               fluidRow(column(12,uiOutput("NodesOfInterest"))),
                               fluidRow(
                                 column(12, checkboxGroupInput("sub_checkGroup", 
                                                               label = "Object Network", 
                                                               choices = list("Nanomaterials" = "nano", 
                                                                              "Drugs" = "drugs", "Chemical" = "chemical",
                                                                              "Disease" = "disease"),selected = c("nano","drugs","chemical","disease")))
                               ),
                               fluidRow(column(12,selectInput("sub_Edges",
                                                              label="Select the type of edges",
                                                              choices = list("Positive" = "P",
                                                                             "Negative" = "N",
                                                                             "Both" = "B"),selected="B"))
                               ),
                               fluidRow(
                                 column(12, selectInput("sub_ATCGroup", 
                                                        label ="Drugs ATC code", multiple = TRUE,
                                                        choices = ATC_choice_list,selected = "A"))
                               ),
                               fluidRow(
                                 column(12, selectInput("sub_ChemicalGroup", 
                                                        label ="Chemical Classes", multiple = TRUE,
                                                        choices = chemical_choice_list,
                                                        selected = "Amino acids"))
                               ),
                               fluidRow(
                                 column(9, sliderInput("sub_repulserad", label = "Repulseration Strenght",
                                                       min = 100, max = 10000, value = 100,step=1))
                               ),
                               fluidRow(
                                 column(9, sliderInput("sub_length", label = "Edge length",
                                                       min = 2, max = 10, value = 2,step=1))
                              )
                             ), #end Sidebar
                             mainPanel("",
                                tabsetPanel( 
                                       tabPanel("Items Subnetwork", 
                                               wellPanel(forceNetworkOutput("Subnetwork_plot")),
                                               wellPanel(plotOutput("dig_dist")),       
                                               wellPanel(plotOutput("Subnetwork_plot_statistic"))
                                                
                                       ), 
                                       tabPanel("Genes", 
                                                wellPanel(forceNetworkOutput("gene_Subnetwork_plot")),
                                                wellPanel(plotOutput("gene_Subnetwork_plot_statistics"))
                                                
                                       )
                                )
                             )#End mainPanel
                    ),#end tabPanel subnetwork
                    tabPanel("Patterns",
                             fluidRow(column(4,uiOutput("NetworkPattern"))),
                             fluidRow(column(4,wellPanel(plotOutput('xx', height = 500))),
                                      column(8,wellPanel(DT::dataTableOutput('clique_data_table')))
                             ),
                             fluidRow(column(4,selectInput("plotTripel",label = "Plot of Association Frequencies",
                                                           choices = list("Disease-Nano-Drug" = 1,
                                                                          "Disease-Nano-Chemical"= 2,
                                                                          "Disease-Nano"=3,
                                                                          "Disease-Drug"=4,
                                                                          "Disease-Chemical"=5,
                                                                          "Chemical-Nano" = 6,
                                                                          "Chemical-Drug" = 7,
                                                                          "Nano-Drug" = 8),selected = 1)),
                                       column(4,numericInput(inputId = "percentuale",label = "% of elements to show",
                                                            value = 10,min = 1,max = 100,step=5)),
                                      column(4,uiOutput("NodesOfInterest_items"))
                             ),
                             fluidRow(
                               column(4,"  "),
                               column(8,wellPanel(plotOutput('ggplot',height=900)))
                             )
                             
                            )#End tabpanel Pattern
         ),
          tabPanel("Phenotypic Network",
                   wellPanel(
                     p("The network of interaction between phenotypic entities is shown in this page"),
                     htmlOutput("checkLOGIN_entities")
                   ),
                   uiOutput("items_net_page1"),
                   uiOutput("items_net_page2")
          ),
          tabPanel("Gene Network",
                   wellPanel(
                       p("The network of interaction of the genes is shown."),
                       p("Two genes are connected in the network if they share a significantly amount of phenotipic entities."),
                       htmlOutput("checkLOGIN_gene")
                   ),
                   fluidRow(uiOutput("gene_net_page1"),
                            uiOutput("gene_net_page2"))
                   
#                    fluidRow(uiOutput("gene_net_page3"),
#                             uiOutput("gene_net_page4"))
                   
                   
          ),
          tabPanel("Demo",
                   wellPanel(
                     fluidRow(
                       h4("Description")
                     ),
                     fluidRow(
                       p("INSIdE nano is a graphical tool  
                         that  highlight connections between objects of biological relevance like :"),
                       #uiOutput("myList"),
                       HTML("<ul>
                            <li>Nanomaterials</li>
                            <li>Drugs</li>
                            <li>Chemicals</li>
                            <li>Diseases</li>
                            <li>Genes</li>
                            </ul>"),
                       p("Data are stored in a network based structure, where the nodes are nanomaterials, drugs, 
                         chemicals and diseases. The similarity between the nodes is 
                         evaluated based on gene expression value."),
                        p("When loaded "),
                       img(src='home.png',width = "80%"),
                        p("From the home page the following operation can be performed:"),
                      img(src="start_button.png",width = "80%"),
                       HTML("<ul>
                            <li>Button Query or tab Start Analisys allows the user to query the network in order to find interactions between items</li>
                            <li>Button Browse Item Network or tab Networks allows the user to visualize the full database of apply some filter</li>
                            <li>Button Browse Gene Network or tab Gene Network allows the user to visualiza gene network interaction. Genes are grouped by KEGG Pathways.</li>
                          </ul>"),
                       p("In the query section the user can choose between two different kind of query. The former, called free query, allows the user to choose which intems he wants to analyze and a threshold that indicates how strongly the items are connected"),
                       img(src = "free_analysis_multiple.png",width = "80%"),
                       p("The system will find out all the connection between the items and their neigboorg in the network. Then it will evaluate the different kinds of connection between nanomaterials, drugs, diseases and chemicals."),
                       p("The latter, called Filter Analysis, allows the user to query the network by applying different filters. The user can specify each items for each data type, can specify how strongly they are connected and how many of them must be interacting each other in the final results."),
                       img(src = "start_filter_analisys.png",width = "80%"),
                       p("Query results are available in tab menu Network Analysis. This menu is divided into two sub-menu. The former (Items Subnetwor) is a visualization tool where the subnetwork related to the query is displayed."),
                       img(src = "network_analisys.png",width = "80%"),
                       #img(src = "network_analisys_stat.png",width = "80%"),
                       p("The latter (Patterns) shows how object are connected each other. Patterns are categorized based on connection properties (positive/negative connection) between the items. The tool will give a graphics representation of the pattern and a list of all the different structure in the network."),
                       img(src = "patterns.png",width = "80%"),
                       p("The tool displayed provides statistics about connection frequencies between couple or triples of object."),
                       img(src = "patterns_statistics.png",width = "80%"),
                       p("In the Network section the user can visualize the interaction network of all nanomaterials, drugs, chemicals and disease. He can apply filter to the network in order to visualize a subnetwork of his interest."),
                       img(src = "complete_network.png",width = "80%"),
                       p("In the Gene Network section the user can visualize the network of all genes. He can apply filter to the network in order to visualize a subnetwork of his interest. Connection by the genes are evaluated based on how many items do they share in their neighborhood."),
                      img(src = "gene_network.png",width = "80%")
                     )
                     
                     )#end well panel
                   
          )
))
 
    
