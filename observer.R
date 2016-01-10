set_observer = function(input,output,session){
  observe({
    if(input$action1 > 0){
      print('1')
      session$sendCustomMessage("myCallbackHandler", "1")
    }
  })
  observe({
    if(input$action2 > 0){
      print('2')
      session$sendCustomMessage("myCallbackHandler", "2")
    }
  })
  observe({
    if(input$action3 > 0){
      print('3')
      session$sendCustomMessage("myCallbackHandler", "3")
    }
  })
  observe({
    if(input$action4 > 0){
      print('5')
      session$sendCustomMessage("myCallbackHandler", "5")
    }
  })
  
#   observe({ #gene query
#     if(input$action5 > 0){
#       print('6')
#       session$sendCustomMessage("myCallbackHandler", "6")
#     }
#   })
  
  observe({
    if(input$Demo > 0){
      print('4')
      session$sendCustomMessage("myCallbackHandler", "4")
    }
  })

}


