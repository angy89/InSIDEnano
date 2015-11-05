setwd("~/Scrivania/")

pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest("pbapply")
pkgTest("igraph")
pkgTest("plotrix")
pkgTest("plyr")
pkgTest("ggplot2")

if(!require(shiny)){
  if (!require("devtools"))
    install.packages("devtools")
  devtools::install_github("rstudio/shiny")
}

if(!require(DT)){
  if (!require("devtools"))
    install.packages("devtools")
  devtools::install_github('rstudio/DT')
}


if(!require(networkD3)){
  if (!require("devtools"))
    install.packages("devtools")
  devtools::install_github('christophergandrud/networkD3')}

library(shiny)
library(DT)
library(igraph)
library(networkD3)
library(plotrix) 
library(plyr)
library(ggplot2)

runApp("inside_nano")