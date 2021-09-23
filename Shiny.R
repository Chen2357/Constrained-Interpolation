ui <- fluidPage(
    sidebarLayout(
        sidebarPanel(
            textInput("xvalue", "Enter x values (comma delimited)", placeholder = "1,2,3,4,5,6"),
            textInput("yvalue", "Enter y values (comma delimited)", placeholder = "1,2,3,4,5,6"),
            textInput("upper", "Upper Bound", placeholder = "10"),
            textInput("lower", "Low Bound", placeholder = "-10"),
            actionButton("go" ,"Interpolate", class = "btn btn-primary"),
            actionButton("view" ,"View", class = "btn btn-primary"),
            actionButton("stop", "Stop")
        ),
        mainPanel(
            tags$h2("Range Restricted C^2 Interpolant"),
            plotOutput("plot")
        )
    )
)

server <- function(input, output) {

    interpolate <- eventReactive(input$go, {
        x <- as.numeric(unlist(strsplit(input$xvalue,",")))
        y <- as.numeric(unlist(strsplit(input$yvalue,",")))
        high <- as.numeric(input$upper)
        low <- as.numeric(input$lower)

        interpolation <- rrinterpolate(x, y, low, high)
        rrplot(interpolation, x = x, y = y, limits = c(low, high), autodiff = TRUE)
    })

    observeEvent(input$view, {
        x <- as.numeric(unlist(strsplit(input$xvalue,",")))
        y <- as.numeric(unlist(strsplit(input$yvalue,",")))
        high <- as.numeric(input$upper)
        low <- as.numeric(input$lower)

        interpolation <- rrinterpolate(x, y, low, high)
        View(interpolation)
    })

    output$plot <- renderPlot({
        interpolate()
    })

    observeEvent(input$stop, {
        stopApp()
    })
}

shinyApp(
    ui = ui,
    server = server
)
