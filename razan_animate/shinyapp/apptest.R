library(shiny)
source('animate_functions.R')

##winry's ip addr
options(shiny.host = "10.161.48.62")


ui <- fluidPage(

  ##App title ----
  titlePanel("Scatterplot animation"),

  ##Sidebar layout with input and output definitions ----
  sidebarLayout(

    ##Sidebar panel for inputs ----
    sidebarPanel(

      ##live cell data
      fileInput("live", "Live cell data",
                multiple = TRUE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),

      # Horizontal line ----
      tags$hr(),

      ##dead cell file
      fileInput("live", "Dead cell data",
                multiple = TRUE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),

      numericInput("totframes", "Total Frames", value = 0, min = 1),
      numericInput("lastframe", "Last Frame", value = 0, min = 1),
      numericInput("startslow", "Slow Motion Start (leave at 0 for no slow mo)", value = 0),
      numericInput("endslow", "Slow Motion End (leave at 0 of no slow mo)", value = 0)

      tags$hr(),
      actionButton(
          inputId = "submit_loc",
          label = "Submit"
      ),
      ),
    ##Main panel for displaying outputs ----
    mainPanel(
        
        'Hello world',
        tags$hr(),
        downloadButton("download", "Download .mp4")

    )

  )
)


server <- function(input, output, session) {
    output$plot1 <- renderImage({
    ## A temp file to save the output.
    ## This file will be removed later by renderImage
    outfile <- tempfile(fileext='.mp4')

    ##now make the animation
    p = ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, 
      color = continent)) + geom_point() + scale_x_log10() +
      transition_time(year) # New
}
shinyApp(ui, server)


