heatmapSubplotUI <- function(id, dpi, width_val, height_val) {
    
    ns <- NS(id)
    
    mult <- as.numeric(dpi) / 100
    
    min.val <- 100 * mult
    max.val <- 3200 * mult
    width.val <- min(width_val * mult, max.val)
    height.val <- min(height_val * mult, max.val)
    
    ## insert a single heatmap
    fluidRow(
        column(
            4,
            
            sliderInput(ns('height'), 'Figure height (px)',
                        min.val, max.val, height.val, 10),
            
            sliderInput(ns('width'), 'Figure width (px)',
                        min.val, max.val, width.val, 10)
       ),
       column(
           8,
           uiOutput(ns('plot1'))
       )
       
    )
}

heatmapSubplot <- function(input, output, session, plot1, i, dpi) {
    force(i)
    
    output$plot1 <- renderUI({
        tagList(
            renderPlot({
                draw(plot1()[[i]], 
                     merge_legends = T,
                     align_heatmap_legend = 'heatmap_top',
                     legend_grouping = 'original'
                     )
            },
            res = dpi,
            width = input$width,
            height = input$height
            )
        )
    })
    
    ## https://stackoverflow.com/questions/48882427/how-to-store-the-returned-value-from-a-shiny-module-in-reactivevalues
    res <- reactiveValues()
    observe({ 
        res$width <- input$width
        res$height <- input$height
        })
    return(res)
}