library(shiny)
library(readr)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
library(ggplot2)
library(gtools)
require(ggtext)

### Functions ####
# Function for loading buscos
read_busco <- function(buscoFile){
    read_tsv(buscoFile,
             col_names = c("Busco_id", "Status", "Sequence",
                           "start", "end", "strand", "Score", "Length",
                           "OrthoDB_url", "Description"),
             col_types = c("ccciicdicc"),
             comment = "#") %>%
        filter(Status == "Complete") %>%
        select(Busco_id, Sequence, start, end)
}

cols <- c("A" = "#af0e2b", "B" = "#e4501e",
          "C" = "#4caae5", "D" = "#f3ac2c",
          "E" = "#57b741", "N" = "#8880be",
          "X" = "#81008b")

nigonDict <- read_tsv("gene2Nigon_busco20200927.tsv.gz",
                      col_types = c(col_character(), col_character()))

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("Paint your genome by Nigon units"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput("file1", "Choose BUSCO full table TSV file",
                      accept = c("tsv")),
            numericInput("windowSize", 
                         "bin size for BUSCO loci", 
                         value = 500000),
            numericInput("minLoci", 
                         "remove scaffolds with less than this many BUSCO loci", 
                         value = 10),
            numericInput("width", 
                         "width for plot in PDF", 
                         value = 5),
            numericInput("height", 
                         "height of plot in PDF", 
                         value = 6),
            textInput("species", "Plot title", value = "Genus_species"),
            selectInput("facetOrientation", 
                        "Orientation for facets", 
                        choices = c("Vertical" = "vertical", "Horizontal" = "horizontal"), 
                        selected = "vertical"),
            selectInput("plotTheme", 
                        "Plot theme", 
                        choices = c("Minimal" = "minimal", "Standard" = "standard"), 
                        selected = "minimal"),
            downloadButton('downloadPlot', 'Download Plot')
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot")
        )
    )
)



# Define server logic required to draw a histogram
server <- function(input, output) {
    dat <- reactive({
        infile <- input$file1
        if (is.null(infile)) {
            # User has not uploaded a file yet
            return(NULL)
        } else {
            opt = list(windowSize = input$windowSize, 
                       minimumGenesPerSequence = input$minimumGenesPerSequence, 
                       species = input$species)
            
            spName <- paste0("*", sub("_", " ", opt$species), "*")
            windwSize <- input$windowSize
            
            busco_tsv <- read_busco(infile$datapath)
            fbusco <- left_join(busco_tsv, nigonDict, by = c("Busco_id" = "Orthogroup")) %>%
                mutate(nigon = ifelse(is.na(nigon), "-", nigon),
                       stPos = start) %>%
                filter(nigon != "-")
            
            consUsco <- group_by(fbusco, Sequence) %>%
                mutate(nGenes = n(),
                       mxGpos = max(stPos)) %>%
                ungroup() %>%
                filter(nGenes > input$minLoci, mxGpos > windwSize * 2)
            
            
            grp_busco <- group_by(consUsco, Sequence) %>%
                mutate(ints = as.numeric(as.character(cut(stPos,
                                                          breaks = seq(0, max(stPos), windwSize),
                                                          labels = seq(windwSize, max(stPos), windwSize)))),
                       ints = ifelse(is.na(ints), max(ints, na.rm = T) + windwSize, ints)) %>%
                count(ints, nigon) %>%
                ungroup() %>%
                mutate(scaffold_f = factor(Sequence,
                                           levels = mixedsort(unique(Sequence))))
            results = list(grp_busco = grp_busco, spName = spName, windwSize = windwSize)
            return(results)
        }
    })
    
    plotInput <- reactive({
        if (is.null(input$file1)) {
            # User has not uploaded a file yet
            return(NULL)
        } else {
            mydata <- dat()
            windwSize <- mydata$windwSize
            if (input$facetOrientation == "vertical") {
                p <- ggplot(mydata$grp_busco, aes(fill=nigon, y=n, x=ints-windwSize)) + 
                    facet_grid(scaffold_f ~ ., switch = "y")
                legend_position <- "right"
                strip_angle <- 0
                gncol <- 1
            } else {
                p <- ggplot(mydata$grp_busco, aes(fill=nigon, y=n, x=ints-windwSize)) + 
                    facet_grid(. ~ scaffold_f)
                legend_position <- "bottom"
                strip_angle <- 90
                gncol <- 7
            }
            p <- p + geom_bar(position="stack", stat="identity") +
                ggtitle(mydata$spName) +
                # theme_minimal() +
                scale_y_continuous(breaks = scales::pretty_breaks(4),
                                   position = "right",
                                   expand = expansion(mult = c(0, 0.1))) +
                scale_x_continuous(labels = label_number_si()) +
                scale_fill_manual(values = cols) +
                guides(fill = guide_legend(ncol = gncol,
                                           title = "Nigon")) 
            # Apply selected theme
            if (input$plotTheme == "minimal") {
                p <- p + theme_minimal()
            } else {
                p <- p + theme(panel.background = element_rect(fill = "white", color = "black"),
                               panel.grid.major = element_line(color = "grey90"),
                               panel.grid.minor = element_line(color = "grey90"))
            }
            p <- p + theme(axis.title.y=element_blank(),
                           axis.title.x=element_blank(),
                           panel.border = element_blank(),
                           legend.position = legend_position,
                           strip.text.y.left = element_text(angle = strip_angle),
                           text = element_text(size=9),
                           plot.title = ggtext::element_markdown()
            )
            
            return(p)
        }
    })
    
    
    output$distPlot <- renderPlot({
        print(plotInput())
    })
    
    output$downloadPlot <- downloadHandler(
        filename = function() { paste("nema_nigon_painted.pdf", sep='') },
        content = function(file) {
            ggsave(file, plotInput(), width = input$width,
                   height = input$height)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
