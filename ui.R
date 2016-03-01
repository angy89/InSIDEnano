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
                                         column(7, textInput("username", "Username",value="user2015", width = NULL))
                                       ),
                                       fluidRow(
                                         column(7, passwordInput("password", "Password",value="helsinki2015", width = NULL))
                                       ),
                                       fluidRow(
                                         column(7,actionButton(inputId = "login",label  = "Login",icon("circle-arrow-right", lib = "glyphicon"))),
                                         singleton(
                                           tags$head(tags$script(src = "message-handler.js"))
                                         ),
                                         column(5,textOutput("loglog"))
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
                                                             if(typeMessage == 6){
                                                             $("a:contains(Gene Query)").click();
                                                             }
                                                             if(typeMessage == 4){
                                                             $("a:contains(Demo)").click();
                                                             }
                                                             });
                                                             ')),
                                     wellPanel(
                                       fluidRow(column(6, div(style="display:inline-block",actionButton("action1", label = "Free Query",width = 350, icon=icon("search",lib = "glyphicon")), style="align:center")),
                                                
                                                # column(4, div(style="display:inline-block",actionButton("action5", label = "Gene Query",width = 400, icon=icon("search",lib = "glyphicon")), style="align:center")),
                                                
                                                column(6, div(style="display:inline-block",actionButton("action4", label = "Conditional Query",width = 350, icon=icon("search",lib = "glyphicon")), style="align:center")))
                                     ),
                                     wellPanel(
                                       fluidRow(column(6, div(style="display:inline-block",actionButton("action2", label = "Browse Phenotypic Network",width = 350,icon=icon("search",lib = "glyphicon")), style="align:center")),
                                                column(6, div(style="display:inline-block",actionButton("action3", label = "Browse Gene Network",width = 350,icon=icon("search",lib = "glyphicon")), style="align:center")))
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
                   #                    tabPanel("Load Query",
                   #                             fluidRow(
                   #                               wellPanel(
                   #                                 HTML("</br>This area allow the user to load results of the past query</br>"),
                   #                                 htmlOutput("checkLOGIN_load")
                   #                               )  
                   #                             ),
                   #                             fluidRow(
                   #                               wellPanel(
                   #                                 fluidRow(actionButton(inputId = "Refresh_query_table",label  = "Refresh Table",icon("circle-arrow-right", lib = "glyphicon"))),
                   #                                 fluidRow(
                   #                                   column(10,DT::dataTableOutput('previous_query')),
                   #                                   column(2,uiOutput("LoadQuery"))
                   #                                   
                   #                                 )
                   #                               )          
                   #                             )
                   #                             
                   #                    ),
                   tabPanel("Single Query",
#                             fluidRow(
#                               column(6,
#                                      fluidRow(
#                                        wellPanel(
#                                          HTML("</br>The Free Analysis allows to extrapolate subnetworks of nodes connected to one or more query items.</br>"),
#                                          HTML("The results of the query depends on the strenght of similarity. </br>
#                                         In case of multiple query items, the final subnetwork will be the union of all the subtenetwork related to each item.
#                                         </br></br></br>")
#                                          # htmlOutput("checkLOGIN_free")
#                                           )
#                                      ),
#                                      fluidRow(
#                                         wellPanel(
#                                             fluidRow(
#                                               column(6, uiOutput("disease_couple",width=500,height=900)),
#                                               column(6, uiOutput("th_slider_couple"))
#                                             ),
#                                             fluidRow(
#                                                column(4, uiOutput("Go_couple")),
#                                                column(4, uiOutput("refresh_couple"))                                         
#                                             )
#                                         )
#                                      )
#                               )
#                             ),
                            fluidRow(
                              wellPanel(
                                HTML("</br>The Free Analysis allows to extrapolate subnetworks of nodes connected to one or more query items.</br>"),
                                HTML("The results of the query depends on the strenght of similarity. </br>
                                                                    In case of multiple query items, the final subnetwork will be the union of all the subtenetwork related to each item.
                                                                    </br></br></br>")
                                # htmlOutput("checkLOGIN_free")
                              )
                            ),
                            fluidRow(
                              wellPanel(
                                fluidRow(
                                  column(6, uiOutput("disease_couple",width=500,height=900)),
                                  column(6, uiOutput("th_slider_couple"))
                                ),
                                fluidRow(
                                  column(4, uiOutput("Go_couple")),
                                  column(4, uiOutput("refresh_couple"))                                         
                                )
                              )
                            ),
                             
#                             fluidRow(
#                               column(3,wellPanel(plotOutput(("nano_couple_boxplot")))),
#                               column(3,wellPanel(plotOutput(("drug_couple_boxplot")))),
#                               column(3,wellPanel(plotOutput(("chemical_couple_boxplot")))),
#                               column(3,wellPanel(plotOutput(("disease_couple_boxplot"))))
#                             ),
                            tabsetPanel(
                              tabPanel("Statistic",
                                       fluidRow(
                                        column(6,wellPanel(plotOutput("boxplot_statistics"))),
                                        column(6,wellPanel(plotOutput("barplot_statistics")))
                                        )
                              ),
                              tabPanel("Tables",
                                       fluidRow(
                                         column(3,"Nanomaterials"),
                                         column(3,"Drugs"),
                                         column(3,"Chemicals"),
                                         column(3,"Diseases")
                                       ),
                                       fluidRow(
                                         column(3,wellPanel(DT::dataTableOutput(("nano_couple_table")))),
                                         column(3,wellPanel(DT::dataTableOutput(("drug_couple_table")))),
                                         column(3,wellPanel(DT::dataTableOutput(("chemical_couple_table")))),
                                         column(3,wellPanel(DT::dataTableOutput(("disease_couple_table"))))
                                       )
                              )
                            )
                            
                   ),
                   tabPanel("Free Query",
                            fluidRow(
                              wellPanel(
                                HTML("</br>The Free Analysis allows to extrapolate subnetworks of nodes connected to one or more query items.</br>"),
                                HTML("The results of the query depends on the strenght of similarity. </br>
                                                 In case of multiple query items, the final subnetwork will be the union of all the subtenetwork related to each item.
                                                 </br></br></br>"),
                                htmlOutput("checkLOGIN_free")
                              )  
                            ),
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
                              column(4,wellPanel(htmlOutput("infoFree")))
                            ),
                            wellPanel(
                              fluidRow(column(8,uiOutput("plotTripel_total"))),
                              fluidRow(uiOutput("NodesOfInterest_totale")),
                              fluidRow( plotOutput('ggplot_totale',height=900))
                            )        
                   ),
                   #                               tabPanel("Gene Query",
                   #                                        fluidRow(
                   #                                          column(8,
                   #                                                 fluidRow(
                   #                                                 wellPanel(
                   #                                                   fluidRow(
                   #                                                     column(6, uiOutput("gene_input")),
                   #                                                     column(6, uiOutput("Go3"))
                   #                                                   )
                   #                                                 )),
                   #                                                 fluidRow(
                   #                                                   column(12,
                   #                                                          wellPanel(
                   #                                                            fluidRow(
                   #                                                              column(3, uiOutput("gene_query_nano_input")),
                   #                                                              column(3, uiOutput("gene_query_drug_input")),
                   #                                                              column(3, uiOutput("gene_query_disease_input")),
                   #                                                              column(3, uiOutput("gene_query_chemical_input"))
                   #                                                            )
                   #                                                          ),
                   #                                                          wellPanel(
                   #                                                            fluidRow(
                   #                                                              column(3, uiOutput("th_slider3")),
                   #                                                              column(3, uiOutput("gene_query_percentuale_somma")),
                   #                                                              column(3, uiOutput("gene_query_nroCliques")),
                   #                                                              column(3, uiOutput("gene_query_clique_type"))
                   #                                                              
                   #                                                            )
                   #                                                          ),
                   #                                                          wellPanel(
                   #                                                            fluidRow(
                   #                                                              column(2,uiOutput("Go4")),
                   #                                                              column(2,uiOutput("gene_query_Refresh"))
                   #                                                              
                   #                                                            )
                   #                                                          )
                   #                                                   )
                   #                                                   
                   #                                                 )
                   #                                          ),
                   #                                          column(4,
                   #                                                 wellPanel(
                   #                                                   HTML("</br>The Free Analysis allows to extrapolate subnetworks of nodes that have influence on a group of genes</br>"),
                   #                                                   HTML("The results of the query depends on the strenght of similarity. </br>
                   #                                                        In case of multiple query items, the final subnetwork will be the union of all the subtenetwork related to each gene
                   #                                                        </br></br></br>"),
                   #                                                   htmlOutput("checkLOGIN_geneQuery")
                   #                                                   ))
                   #                                          )
                   #                                        
                   #                                        wellPanel(
                   #                                          fluidRow(column(8,uiOutput("plotTripel_total"))),
                   #                                          fluidRow(uiOutput("NodesOfInterest_totale")),
                   #                                          fluidRow( plotOutput('ggplot_totale',height=900))
                   #                                        )    
                   #                              ),
                   tabPanel("Conditional Query",
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
                            ),
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
                            #                             fluidRow(
                            #                               wellPanel(
                            #                                 HTML("</br>This area allow the user to load results of the past query</br>"),
                            #                                 htmlOutput("checkLOGIN_load")
                            #                               )  
                            #                             ),
                            fluidRow(
                              wellPanel(
                                fluidRow(
                                  column(10,DT::dataTableOutput('previous_query')),
                                  column(2,
                                         fluidRow(uiOutput("LoadQuery")),
                                         fluidRow(actionButton(inputId = "Refresh_query_table",label  = "Refresh Table",icon("circle-arrow-right", lib = "glyphicon")))
                                  )
                                  
                                )
                              )          
                            )
                   ),
                   
                   navbarMenu("Query Results",
                              tabPanel("Items Subnetwork", 
                                       sidebarPanel(
                                         fluidRow(column(12,uiOutput("NodesOfInterest"))),
                                         #                                          fluidRow(
                                         #                                            column(12, checkboxGroupInput("sub_checkGroup", 
                                         #                                                                          label = "Object Network", 
                                         #                                                                          choices = list("Nanomaterials" = "nano", 
                                         #                                                                                         "Drugs" = "drugs", "Chemical" = "chemical",
                                         #                                                                                         "Disease" = "disease"),selected = c("nano","drugs","disease")))
                                         #                                          ),
                                         #                                          fluidRow(column(12,selectInput("sub_Edges",
                                         #                                                                         label="Select the type of edges",
                                         #                                                                         choices = list("Positive" = "P",
                                         #                                                                                        "Negative" = "N",
                                         #                                                                                        "Both" = "B"),selected="B"))
                                         #                                          ),
                                         fluidRow(
                                           column(12, selectInput("sub_ATCGroup", 
                                                                  label ="Drugs ATC code", multiple = TRUE,
                                                                  choices = ATC_choice_list,selected = "ALL"))
                                           
                                         ),
                                         fluidRow(
                                           column(12, selectInput("sub_ChemicalGroup", 
                                                                  label ="Chemical Classes", multiple = TRUE,
                                                                  choices = chemical_choice_list,
                                                                  selected = "ALL"))
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
                                                            #wellPanel(plotOutput("dig_dist")),       
                                                            wellPanel(plotOutput("Subnetwork_plot_statistic",width = "100%",height = 500))
                                                            
                                                   ), 
                                                   tabPanel("Genes", 
                                                            wellPanel(forceNetworkOutput("gene_Subnetwork_plot")),
                                                            wellPanel(plotOutput("gene_Subnetwork_plot_statistics",width = "100%",height = 800))
                                                            
                                                   )
                                                 )
                                       )#End mainPanel
                              ),#end tabPanel subnetwork
                              tabPanel("Patterns",
                                       fluidRow(uiOutput("NetworkPattern")),
                                       fluidRow(
                                         column(4,wellPanel(plotOutput('xx'))),
                                         column(8,wellPanel(DT::dataTableOutput('clique_data_table')))
                                       ),
                                       
                                       tabsetPanel(
                                         tabPanel("Statistic",
                                                  wellPanel(
                                                    fluidRow(
                                                      column(4,uiOutput("bubbleCoupleChoice")),
                                                      column(4,numericInput(inputId = "percentuale",label = "% of elements to show",value = 10,min = 1,max = 100,step=5))
                                                    ),
                                                    fluidRow(column(8,plotOutput('trendPlot',width ="1200px",height="700px")))#plotlyOutput   
                                                  )
                                         ),#End tabpanel statistic
                                         tabPanel("Enrichment",
                                                  wellPanel(
                                                    fluidRow(
                                                      #fluidRow(wellPanel(DT::dataTableOutput('drugs_chemical_DT'))),
                                                      column(6,
                                                             fluidRow(
                                                               wellPanel(
                                                                 sliderInput("EnrichTh", label = "Connection Strenght",
                                                                             min = 1, max = 99, value = 99,step=1)
                                                               )),
                                                             fluidRow(wellPanel(nanoClusterOutput('enriched_clique', height = 500)))
                                                      ),
                                                      column(6,wellPanel(DT::dataTableOutput('genes_data_table')))
                                                    )
                                                  )
                                         ),#end tabpanel Enrichment
                                         tabPanel("Pubmed Search",
                                                  fluidRow(
                                                    column(4,
                                                           wellPanel(
                                                             helpText("Type a word below and search PubMed to find documents that contain that word in the text.
                                                             You can even type multiple words. You can search authors, topics, any acronym, etc."),
                                                             textInput("text", label = h3("Keyord(s)"), value = "MWCNT Asthma"),
                                                             helpText("You can specify the start and end dates of your search, use the format YYYY/MM/DD"),
                                                             textInput("date1", label = h3("From"),value="1999/01/01"),
                                                             textInput("date2", label = h3("To"),  value = format(Sys.time(), "%Y/%d/%m")),
                                                             actionButton("wordButton","Search"),
                                                             plotOutput('wordPlot')
                                                           )  
                                                    ),
                                                    column(8,wellPanel(DT::dataTableOutput('articleTable')))
                                                  )
                                         )
                                       )#end tabset
                              )#End tabpanel Pattern
                   ),#end query menu
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
                            fluidRow(uiOutput("KEGG_net_page1"),
                                     uiOutput("KEGG_net_page2")),
                            fluidRow(uiOutput("gene_net_page1"),
                                     uiOutput("gene_net_page2"))
                            #                    fluidRow(uiOutput("gene_net_page3"),
                            #                             uiOutput("gene_net_page4"))                   
                   ),
                   tabPanel("Cluster Analysis",
                            tabsetPanel(
                              tabPanel("Nano based clustering",
                                       fluidRow(
                                         wellPanel(HTML("<p>The user can investigate what are the items strongly correlated to each nanomaterial.</p>"),
                                                   htmlOutput("checkLOGIN_NBC"))
                                       ),
                                       sidebarPanel(
                                         fluidRow(
                                           column(6, uiOutput("nano_cluster_input"))
                                         ),
                                         #                               fluidRow(
                                         #                                 column(6, uiOutput("nano_cluster_input")),
                                         #                                 column(6, uiOutput("drug_cluster_input"))
                                         #                               ),
                                         #                               fluidRow(
                                         #                                 column(6, uiOutput("disease_cluster_input")),
                                         #                                 column(6, uiOutput("chemical_cluster_input"))  
                                         #                               ),
                                         fluidRow(
                                           column(9, sliderInput("repulseration_cluster", label = "Repulseration Strenght",
                                                                 min = 100, max = 10000, value = 4000,step=1))
                                         ),
                                         fluidRow(
                                           column(9, sliderInput("edge_length_cluster", label = "Edge length",
                                                                 min = 2, max = 100, value =2 ,step=1))
                                         ),
                                         fluidRow(
                                           column(4, actionButton("action1_cluster", label = "Start", icon=icon("search",lib = "glyphicon")))
                                           
                                         )
                                       ), #end Sidebar
                                       mainPanel(
                                         wellPanel(nanoClusterOutput("cluster_output"))
                                       )        
                                       
                              ),
                              tabPanel("Nano clustering",
                                       fluidRow(
                                         wellPanel(HTML("<p>The user can investigate how nanomaterials are grouped each other and which genes are most important to each of them.</p>"),
                                                   htmlOutput("checkLOGIN_NC"))
                                       ),
                                       sidebarPanel(
                                         fluidRow(
                                           numericInput(inputId = "GP",label = "% of genes to show",
                                                        value = 0.05,min = 0,max = 100,step=0.01)),
                                         fluidRow(wellPanel(DT::dataTableOutput('node_DT')))
                                       ),
                                       mainPanel(
                                         #                                       tags$head(tags$style('.node {
                                         #                                                             cursor: pointer;
                                         #                                                            }
                                         #                                                            
                                         #                                                            .node circle {
                                         #                                                            fill: #fff;
                                         #                                                            stroke: steelblue;
                                         #                                                            stroke-width: 1.5px;
                                         #                                                            }
                                         #                                                            
                                         #                                                            .node text {
                                         #                                                            font: 12px sans-serif;
                                         #                                                            }
                                         #                                                            
                                         #                                                            .link {
                                         #                                                            fill: none;
                                         #                                                            stroke: #ccc;
                                         #                                                            stroke-width: 0.8px;
                                         #                                                            }')),
                                         collapsibleTreeOutput("treeNANO")
                                       )        
                              ),
                              #                      tabPanel("Drug clustering",radialNetworkOutput("drugs_radial",width = 4000,height = 4000)),
                              #                      tabPanel("Disease clustering",radialNetworkOutput("disease_radial",width = 4000,height = 4000)),
                              #                      tabPanel("Chemical clustering",radialNetworkOutput("chemical_radial",width = 6000,height = 6000))
                              
                              tabPanel("Drug clustering",collapsibleTreeOutput("treeDRUGS",width = 4000,height = 4000)),
                              tabPanel("Disease clustering",collapsibleTreeOutput("treeDISEASE",width = 4000,height = 4000)),
                              tabPanel("Chemical clustering",collapsibleTreeOutput("treeCHEMICAL",width = 6000,height = 6000))
                            )
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
                                p("When loaded the home page will appear. On the left user can login to the service. Whitout login analisy can't be performed."),
                                img(src='home.png',width = "80%"),
                                p("From the home page the following operation can be performed:"),
                                img(src="start_button.png",width = "80%"),
                                HTML("<ul>
                            <li>Free Query allows the user to query the network in order to find interactions between items</li>
                            <li>Conditional Query allows the user to query the network in order to find interactions between items depending on specific thresholds</li>
                            <li>Browse Phenotypic Network or tab Networks allows the user to visualize the full database of apply some filter</li>
                            <li>Browse Gene Network or tab Gene Network allows the user to visualiza gene network interaction. Genes are grouped by KEGG Pathways.</li>
                          </ul>"),
                                p("The Free query tool, allows the user to choose which intems he wants to analyze and a threshold that indicates how strongly the items are connected"),
                                img(src = "free_analysis_multiple.png",width = "80%"),
                                p("The system will find out all the connection between the items and their neigboorg in the network. Then it will evaluate the different kinds of connection between nanomaterials, drugs, diseases and chemicals."),
                                p("The Filter Analysis tool, allows the user to query the network by applying different filters. The user can specify each items for each data type, can specify how strongly they are connected and how many of them must be interacting each other in the final results."),
                                img(src = "start_filter_analisys.png",width = "70%"),
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


