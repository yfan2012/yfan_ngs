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
      fileInput("dead", "Dead cell data",
                multiple = TRUE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),

      numericInput("totframes", "Total Frames (required)", value = 0, min = 1),
      numericInput("lastframe", "Last Frame (leave at 0 for all frames)", value = 0, min = 1),
      numericInput("startslow", "Slow Motion Start (leave at 0 for no slow mo)", value = 0),
      numericInput("endslow", "Slow Motion End (leave at 0 of no slow mo)", value = 0),

      tags$hr(),
      actionButton(
          inputId = "submit_loc",
          label = "Submit"),
      ),
    
    ##Main panel for displaying outputs ----
    mainPanel(
        'Hello world',
        tags$hr(),
        ##tags$video(id="videoID", type = "video/mp4",src = "output.mp4", controls = "controls")
        imageOutput("plot1"),
        tags$hr(),
        downloadButton("download", "Download .mp4")
        
    )

  )
)


server <- function(input, output, session) {
    deaddata=read_csv(input$dead$datapath, col_names=F)
    livedata=read_csv(input$live$datapath, col_names=F)
    
    deadpercell=sep.colors(deaddata, input$totframes, 'dead')
    livepercell=sep.colors(livedata, input$totframes, 'live')
    
    a=make_movie(livepercell, deadpercell)

    output$vid <- renderUI({
    
        anim_save("outfile.mp4", a)
        list(src = "outfile.mp4",
             contentType = 'image/gif'
             )}, deleteFile = TRUE)
}

shinyApp(ui, server)





    
    output$download <- downloadHandler(
        filename = function() {
            paste0(, ".mp4")
        },
        content = function(file) {
            write.csv(data(), file)
        }
    )
