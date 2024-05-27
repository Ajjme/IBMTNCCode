

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


ck1 <- read_csv("k1.csv")
cg1 <- read_csv("g1.csv")
ck2 <- read_csv("k2.csv")
cg2 <- read_csv("g2.csv")
ck3 <- read_csv("k3.csv")
cg3 <- read_csv("g3.csv")
cgka <- read_csv("gka.csv")
cm1 <- read_csv("m1.csv")
ct1 <- read_csv("t1.csv")
cm2 <- read_csv("m2.csv")
ct2 <- read_csv("t2.csv")
cm3 <- read_csv("m3.csv")
ct3 <- read_csv("t3.csv")
cmt1 <- read_csv("mt1.csv")
cmt23 <- read_csv("mt23.csv")
k1 <- ggplot(data=ck1, aes(week, wday, fill = waterLevel)) +
    geom_tile(color = "white") +
    facet_wrap(~ year, ncol = 1) +
    scale_fill_gradientn(colours = c("#F0EFF0", "#AAB6FB","#6096FD","#031B88")) +
    scale_x_continuous(name = "Week of the Year", breaks = 1:52, labels = 1:52) +
    labs(subtitle = "Change in Water Level (m) in Karurumo(C)", y = "Day of week") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8, face = "bold"))
g1 <- ggplot(data=cg1, aes(week, wday, fill = waterLevel)) +
    geom_tile(color = "white") +
    facet_wrap(~ year, ncol = 1) +
    scale_fill_gradientn(colours = c("#F0EFF0", "#AAB6FB","#6096FD","#031B88")) +
    scale_x_continuous(name = "Week of the Year", breaks = 1:52, labels = 1:52) +
    labs(subtitle = "Change in Water Level (m) in Githambara(T)", y = "Day of week") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8, face = "bold"))
k2 <- ggplot(data=ck2,aes(x=Date,y=TSS))+
    geom_point()+
    #labs(title="TSS change in Karurumo(C)") +
    ylab("Average TSS (mg/l)")+
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    xlab("Year")
g2 <- ggplot(data=cg2,aes(x=Date,y=TSS))+
    geom_point()+
    #labs(title="TSS change in Githambara(T)") +
    ylab("Average TSS (mg/l)")+
    xlab("Year")
k3 <- ggplot(data=ck3,aes(x=Date,y=Turbidity))+
    geom_point()+
    labs(title="Turbidity change in Karurumo(C)") +
    ylab("Average Turbidity (NTU)")+
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    xlab("Year")
g3 <- ggplot(data=cg3,aes(x=Date,y=Turbidity))+
    geom_point()+
    labs(title="Turbidity change in Githambara(T)") +
    ylab("Average Turbidity (NTU)")+
    xlab("Year")
gka <- ggplot(data=cgka,aes(x=Year,y=Average,color=watershed))+
    geom_point()+
    geom_line()+
    labs(title="Annual average water level change in Githambara and Karurumo") +
    ylab("Water level(m)")+
    xlab("year")
m1 <- ggplot(data =cm1, aes(week, wday, fill = waterLevel)) +
    geom_tile(color = "white") +
    facet_wrap(~ year, ncol = 1) +
    scale_fill_gradientn(colours = c("#F0EFF0", "#AAB6FB","#6096FD","#031B88")) +
    scale_x_continuous(name = "Week of the Year", breaks = 1:52, labels = 1:52) +
    labs(subtitle = "Change in Water Level (m) in Mbogiti(T)", y = "Day of week") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8, face = "bold"))
t1 <- ggplot(data=ct1, aes(week, wday, fill = waterLevel)) +
    geom_tile(color = "white") +
    facet_wrap(~ year, ncol = 1) +
    scale_fill_gradientn(colours = c("#F0EFF0", "#AAB6FB","#6096FD","#031B88")) +
    scale_x_continuous(name = "Week of the Year", breaks = 1:52, labels = 1:52) +
    labs(subtitle = "Change in Water Level (m) in Thika-Valley(C)", y = "Day of week") +
    theme(axis.text = element_text(size = 5), axis.title = element_text(size = 8, face = "bold"))
m2 <- ggplot(data=cm2,aes(x=Date,y=TSS))+
    geom_point()+
    #labs(title="TSS change in Mbogiti") +
    ylab("Average TSS (mg/l)")+
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    xlab("Year")
t2 <- ggplot(data=ct2,aes(x=Date,y=TSS))+
    geom_point()+
    #labs(title="TSS change in Thika-Valley") +
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    ylab("Average TSS (mg/l)")+
    xlab("Year")
m3 <- ggplot(data=cm3,aes(x=Date,y=Average_by_day_and_year))+
    geom_point()+
    #labs(title="Turbidity change in Mbogiti") +
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    ylab("Average Turbidity(NTU)")+
    xlab("Year")
t3 <- ggplot(data=ct3,aes(x=Date,y=Average_by_day_and_year))+
    geom_point()+
    #labs(title="Turbidity change in Thika-valley") +
    scale_x_date(date_breaks= "1 year",date_labels="%Y")+
    ylab("Average Turbidity(NTU)")+
    xlab("Year")
mt1 <- ggplot(data=cmt1,aes(x=Year,y=Average,color=watershed))+
    geom_line()+
    geom_point()+
    labs(title="Annual average water level change in Thika-Valley(C) and Mbogiti(T)") +
    ylab("Water level(m)")+
    xlab("year")
mt2 <- ggplot(data=cmt23,aes(x=year,y=TSS,color=watershed))+
    geom_line()+
    geom_point()+
    labs(title="Annual average TSS change (mg/l) in Thika-Valley(C) and Mbogiti(T)") +
    ylab("Average TSS (mg/l)")+
    xlab("year")
mt3 <- ggplot(data=cmt23,aes(x=year,y=Turbidity,color=watershed))+
    geom_line()+
    geom_point()+
    labs(title="Annual average Turbidity change (NTU) in Thika-Valley(C) and Mbogiti(T)") +
    ylab("Average Turbidity(NTU)")+
    xlab("year")
pk1 <- ggplotly(
    p = k1,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pg1 <- ggplotly(
    p = g1,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pk2 <- ggplotly(
    p = k2,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pg2 <- ggplotly(
    p = g2,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pk3 <- ggplotly(
    p = k3,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pg3 <- ggplotly(
    p = g3,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
gka <- ggplotly(
    p = gka,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pm1 <- ggplotly(
    p = m1,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pt1 <- ggplotly(
    p = t1,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pm2 <- ggplotly(
    p = m2,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pt2 <- ggplotly(
    p = t2,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pm3 <- ggplotly(
    p = m3,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pt3 <- ggplotly(
    p = t3,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pmt1 <- ggplotly(
    p = mt1,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pmt2 <- ggplotly(
    p = mt2,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
pmt3 <- ggplotly(
    p = mt3,
    width = NULL,
    height = NULL,
    tooltip = "all",
    dynamicTicks = FALSE,
    layerData = 1,
    originalData = TRUE,
    source = "A"
)
server <- function(input, output){
    
    output$ui1 <- renderUI({
        if (is.null(input$action))
            return()
        switch (input$action,
                'Githambara and Karurumo (Seasonally)' = selectInput("category", label = h4("Select Category"), 
                                                                     choices = c("Water Level vs. Precipitation Significance Testing", 
                                                                                 "Water Level vs. Precipitation Linear Regression", 
                                                                                 "Water Level vs. Others" ), 
                                                                     selected = 1) 
        ) 
    })
    output$ui2 <- renderUI({
        if (is.null(input$action))
            return()
        switch (input$action,
                'Githambara and Karurumo (Yearly)' = selectInput("category", label = h4("Select Category"), 
                                                                 choices = c("Water Level", "Total Suspended Solids", "Turbidity" ), 
                                                                 selected = 1) 
        ) 
    })
    output$ui3 <- renderUI({
        if (is.null(input$action))
            return()
        switch (input$action,
                'Mbogiti and Thika-Valley (Yearly)' = selectInput("category", label = h4("Select Category"), 
                                                                  choices = c("Water Level", "Total Suspended Solids", "Turbidity" ), 
                                                                  selected = 1) 
        ) 
    })
    output$result1 <- renderUI({
        switch (input$action,
                "About" = uiOutput('about'),
                "Githambara and Karurumo (Seasonally)" = uiOutput("gks_title"),
                "Githambara and Karurumo (Yearly)" = uiOutput("gk_title"),
                "Mbogiti and Thika-Valley (Yearly)" = uiOutput("mt_title"),
                "Conclusion" = uiOutput("conclusion_title")
        ) 
    })
    output$result2 <- renderUI({
        if(is.null(input$category)) return()
        if(is.null(input$action)) return()
        if (input$action != "Githambara and Karurumo (Seasonally)")
            return()
        switch (input$category,
                "Water Level vs. Precipitation Linear Regression"=tabsetPanel(tabPanel("Githambara (T)", imageOutput("gks21")),
                                                                              tabPanel("Karurumo (C)", imageOutput("gks22"))),
                "Water Level vs. Others"=tabsetPanel(tabPanel("Githambara (T)", imageOutput("gks31")),
                                                     tabPanel("Karurumo (C)", imageOutput("gks32")))
        )
        
    })
    output$condition1 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Precipitation Significance Testing") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= imageOutput("gks11")
        )
    })
    output$condition2 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Precipitation Significance Testing") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= uiOutput("gks12_title")
        )
    })
    output$condition3 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Precipitation Significance Testing") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= imageOutput("gks12")
        )
    })
    output$condition4 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Precipitation Significance Testing") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= imageOutput("gks13")
        )
    })
    output$condition5 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Precipitation Linear Regression") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= uiOutput("gks2_title")
        )
    })
    output$condition6 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level vs. Others") return()
        switch(
            input$action,
            'Githambara and Karurumo (Seasonally)'= uiOutput("gks3_title")
        )
    })
    output$result3 <- renderUI({
        if(is.null(input$category)) return()
        if(is.null(input$action)) return()
        if (input$action != "Githambara and Karurumo (Seasonally)")
            return()
        switch (input$category,
                "Water Level vs. Precipitation Linear Regression"=tabsetPanel(tabPanel("Githambara (T)", imageOutput("gks23")),
                                                                              tabPanel("Karurumo (C)", imageOutput("gks24")),
                                                                              tabPanel("Combined", imageOutput("gks25"))),
                "Water Level vs. Others"=tabsetPanel(tabPanel("Githambara (T)", imageOutput("gks33")),
                                                     tabPanel("Karurumo (C)", imageOutput("gks34")))
        )
        
    })
    output$result4 <- renderUI({
        if(is.null(input$category)) return()
        if(is.null(input$action)) return()
        if (input$action != "Githambara and Karurumo (Yearly)")
            return()
        switch (input$category,
                "Water Level"=tabsetPanel(tabPanel("Githambara (T)", plotlyOutput("g1_heatmap")),
                                          tabPanel("Karurumo (C)", plotlyOutput("k1_heatmap"))),
                "Total Suspended Solids"=tabsetPanel(tabPanel("Githambara (T)", plotlyOutput("g2_scatter")),
                                                     tabPanel("Karurumo (C)", plotlyOutput("k2_scatter"))),
                "Turbidity"=tabsetPanel(tabPanel("Githambara (T)", plotlyOutput("g3_scatter")),
                                        tabPanel("Karurumo (C)", plotlyOutput("k3_scatter")))
        )
        
    })
    output$result5 <- renderUI({
        if(is.null(input$category)) return()
        if(is.null(input$action)) return()
        if (input$action != "Mbogiti and Thika-Valley (Yearly)")
            return()
        switch (input$category,
                "Water Level"=tabsetPanel(tabPanel("Mbogiti (T)", plotlyOutput("m1_heatmap")),
                                          tabPanel("Thika-Valley (C)", plotlyOutput("t1_heatmap"))),
                "Total Suspended Solids"=tabsetPanel(tabPanel("Mbogiti (T)", plotlyOutput("m2_scatter")),
                                                     tabPanel("Thika-Valley (C)", plotlyOutput("t2_scatter"))),
                "Turbidity"=tabsetPanel(tabPanel("Mbogiti (T)", plotlyOutput("m3_scatter")),
                                        tabPanel("Thika-Valley (C)", plotlyOutput("t3_scatter")))
                
        )
    })
    output$condition7 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level") return()
        switch(
            input$action,
            'Githambara and Karurumo (Yearly)'= uiOutput("gkana_title")
        )
    })
    output$condition8 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level") return()
        switch(
            input$action,
            'Githambara and Karurumo (Yearly)'= plotlyOutput("gkana_line")
        )
    })
    output$condition9 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level") return()
        switch(
            input$action,
            'Githambara and Karurumo (Yearly)'= imageOutput("gkana_test")
        )
    })
    output$result6 <- renderUI({
        switch (input$action,
                "Mbogiti and Thika-Valley (Yearly)" = uiOutput("mtana_title"),
                "About" = uiOutput('about2')
        ) 
    })
    output$result7 <- renderUI({
        switch (input$action,
                "About" = imageOutput("about_image"),
                "Conclusion" = imageOutput("conclusion_image")
        )
    })
    output$condition10 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = plotlyOutput("mt1ana_line")
        )
    })
    output$condition11 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Water Level") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = imageOutput("mt1ana_test")
        )
    })
    output$condition12 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Total Suspended Solids") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = plotlyOutput("mt2ana_line")
        )
    })
    output$condition13 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Total Suspended Solids") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = imageOutput("mt2ana_test")
        )
    })
    output$condition14 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Turbidity") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = plotlyOutput("mt3ana_line")
        )
    })
    output$condition15 <- renderUI({
        if(is.null(input$category)) return()
        if(input$category != "Turbidity") return()
        switch(
            input$action,
            "Mbogiti and Thika-Valley (Yearly)" = imageOutput("mt3ana_test")
        )
    })
    
    
    output$about <- renderUI ({
        h2("Research Question 3 - Evaluating TNC interventions")
    })
    output$about2 <- renderUI ({
        h3("Goal: exploring water variability between sites and how interventions can contribute towards improving water availability.")
    })
    output$gks12_title <- renderUI({ 
        h2(" Analysis - ANOVA Test and T-test")
    })
    output$gks_title <- renderUI({ 
        h2(input$category," in Githambara and Karurumo")
    })
    output$gk_title <- renderUI({ 
        h3("Change in ",input$category," in Githambara and Karurumo")
    })
    output$mt_title <- renderUI({ 
        h3("Change in ",input$category," in Mbogiti and Thika-Valley")
    })
    output$gkana_title <- renderUI({ 
        h3("Analysis - Annual Trend and T-test")
    })
    output$mtana_title <- renderUI({ 
        h3("Analysis - Annual Trend and T-test")
    })
    output$conclusion_title <- renderUI ({
        h1("Conclusion")
    })
    output$gkana_line <- renderPlotly({
        gka
    })
    
    output$conclusion_image <- renderImage({
        list(src = "conclusion_image.png",
             contentType = 'png',
             width = 962,
             height = 528)
    },deleteFile = FALSE)
    
    output$about_image <- renderImage({
        list(src = "about_image.png",
             contentType = 'png',
             width = 660,
             height = 510
        )
    },deleteFile = FALSE)
    
    output$gks11 <- renderImage({
        list(src = "gks11.png",
             contentType = 'png',
             width = 580,
             height = 400
        ) 
    },  deleteFile = FALSE)
    output$gks12 <- renderImage({
        list(src = "gks12.png",
             contentType = 'png',
             width = 340,
             height = 240
        ) 
    },  deleteFile = FALSE)
    output$gks13 <- renderImage({
        list(src = "gks13.JPG",
             contentType = 'jpeg',
             width = 600,
             height = 133
        ) 
    },  deleteFile = FALSE)
    
    output$gks21 <- renderImage({
        list(src = "gks21.jpeg",
             contentType = 'jpeg',
             width = 680,
             height = 400
        ) 
    },  deleteFile = FALSE)
    output$gks22 <- renderImage({
        list(src = "gks22.jpeg",
             contentType = 'jpeg',
             width = 680,
             height = 400
        ) 
    },  deleteFile = FALSE)
    output$gks23 <- renderImage({
        list(src = "gks23.jpeg",
             contentType = 'jpeg',
             width = 800,
             height = 450
        ) 
    },  deleteFile = FALSE)
    output$gks24 <- renderImage({
        list(src = "gks24.jpeg",
             contentType = 'jpeg',
             width = 800,
             height = 450
        ) 
    },  deleteFile = FALSE)
    output$gks25 <- renderImage({
        list(src = "gks25.jpeg",
             contentType = 'jpeg',
             width = 750,
             height = 453
        ) 
    },  deleteFile = FALSE)
    output$gks31 <- renderImage({
        list(src = "gks31.png",
             contentType = 'png',
             width = 816.8,
             height = 256
        ) 
    },  deleteFile = FALSE)
    output$gks32 <- renderImage({
        list(src = "gks32.png",
             contentType = 'png',
             width = 820,
             height = 236
        ) 
    },  deleteFile = FALSE)
    output$gks33 <- renderImage({
        list(src = "gks33.png",
             contentType = 'png',
             width = 927.2,
             height = 293.6
        ) 
    },  deleteFile = FALSE)
    output$gks34 <- renderImage({
        list(src = "gks34.png",
             contentType = 'png',
             width = 930.4,
             height = 267.2
        ) 
    },  deleteFile = FALSE)
    
    
    
    output$gkana_test <- renderImage({
        list(src = "gk_ana.png",
             contentType = 'png',
             width = 750,
             height = 260
        ) 
    },  deleteFile = FALSE)
    output$mtana_title <- renderUI({ 
        h3("Analysis - Annual Trend and T-test")
    })
    output$mt1ana_line <- renderPlotly({
        pmt1
    })
    output$mt1ana_test <- renderImage({
        list(src = "mt1_ana.png",
             contentType = 'png',
             width = 750,
             height = 260
        ) 
    },  deleteFile = FALSE)
    output$mt2ana_line <- renderPlotly({
        pmt2
    })
    output$mt2ana_test <- renderImage({
        list(src = "mt2_ana.png",
             contentType = 'png',
             width = 750,
             height = 260
        ) 
    },  deleteFile = FALSE)
    output$mt3ana_line <- renderPlotly({
        pmt3
    })
    output$mt3ana_test <- renderImage({
        list(src = "mt3_ana.png",
             contentType = 'png',
             width = 750,
             height = 260
        ) 
    },  deleteFile = FALSE)
    
    
    
    output$g1_heatmap <- renderPlotly({
        pg1
    })
    output$k1_heatmap <- renderPlotly({
        pk1
    })
    output$g2_scatter <- renderPlotly({
        pg2
    })
    output$k2_scatter <- renderPlotly({
        pk2
    })
    output$g3_scatter <- renderPlotly({
        pg3
    })
    output$k3_scatter <- renderPlotly({
        pk3
    })
    output$m1_heatmap <- renderPlotly({
        pm1
    })
    output$t1_heatmap <- renderPlotly({
        pt1
    })
    output$m2_scatter <- renderPlotly({
        pm2
    })
    output$t2_scatter <- renderPlotly({
        pt2
    })
    output$m3_scatter <- renderPlotly({
        pm3
    })
    output$t3_scatter <- renderPlotly({
        pt3
    })
}



