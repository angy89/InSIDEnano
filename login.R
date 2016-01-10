validate_login = function(input,output){
  validate(need(input$username != "", "Please insert login name"))
  validate(need(input$password != "", "Please insert password"))
}

login_successful = function(input,output,session){
  LOGGED_IN = TRUE
  if(DEBUGGING)
    cat("Login successfull... \n")
  output$loglog = renderText("Logged!")
  
  session$sendCustomMessage(type = 'testmessage',
                            message = "Login Successfull! ")
}

login_failed = function(input,output,session){
  if(DEBUGGING)
    cat("Login failed... \n")
  output$loglog = renderText("Login Failed!")
  session$sendCustomMessage(type = 'testmessage',
                            message = "Login Failed! ")
}