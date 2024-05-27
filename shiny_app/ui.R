library(shiny)
library(shiny)
library(readr) #read file
library(dplyr) #data manipulation
library(tidyr) #overall
library(lubridate) #read date
library(tibble)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(gridExtra)
library(plotly)
library(tidyverse)


ui <- fluidPage(
    # App title ----
    titlePanel("UCLA practicum -  IBM Open Water Platform"),
    
    ### select bars ####
    sidebarLayout(
        sidebarPanel(
            selectInput("action", label = h4("Select Action"), 
                        choices = c("About","Githambara and Karurumo (Seasonally)", "Githambara and Karurumo (Yearly)", 
                                    "Mbogiti and Thika-Valley (Yearly)","Conclusion"), 
                        selected = 1),
            uiOutput('ui1'),
            uiOutput('ui2'),
            uiOutput('ui3'),
            uiOutput('ui_category'),
            hr()
        ),
        
        # end ####
        mainPanel(
            
            uiOutput('result1'),
            uiOutput('result2'),
            uiOutput('result3'),
            uiOutput('result4'),
            uiOutput('condition1'),
            uiOutput('condition2'),
            uiOutput('condition3'),
            uiOutput('result5'),
            uiOutput('result6'),
            uiOutput('result7'),
            uiOutput('condition4'),
            uiOutput('condition5'),
            uiOutput('condition6'),
            uiOutput('condition7'),
            uiOutput('condition8'),
            uiOutput('condition9'),
            uiOutput('condition10'),
            uiOutput('condition11'),
            uiOutput('condition12'),
            uiOutput('condition13'),
            uiOutput('condition14'),
            uiOutput('condition15'),
            uiOutput('condition16')
            
            
            
            
        )
    )
)
