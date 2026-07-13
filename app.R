# Breeding Strategy Sim App
## App code

# Load Packages and Source Scripts==================================================================
## Load packages
pacman:: p_load(shiny, tidyverse, DT)


## Load functions
source("utils.R")



# Create App========================================================================================
## UI
ui <- sidebarLayout(
  sidebarpanel(
    sliderInput("sld_pop_init", min=50, max=1000, value=100, step=50),
    sliderInput("sld_chromes", min=1, max=21, value=10, step=1),
    sliderInput("sld_genes_per_chrome", min=5, max=100, value=20, step=5),
    sliderInput("sld_herit_init", min=.05, max=.95, value=0.3, step=0.05),
    sliderInput("sld_n_cross", min=10, max=200, value=50, step=10),
    sliderInput("sld_n_select", min=5, max=50, value=20, step=5),
    actionButton("bt_run", "Run Simulation")
  ),
  mainPanel(
    tabsetPanel(
      tabPanel(
        "Genetic Gain Plots",
        plotOutput("plot_gen_gain")
        ),
      tabPanel(
        "Variance/Diversity",
        plotOutput("plot_var_diversity")
      ),
      tabPanel(
        "Comparison Table",
        plotOutput("")
      )
    )

  )
  
)
  
server <- function(input, output, session) {
  
  react_pop_sp <- eventReactive("btn_run", {
    create_founders(
      n_ind=input$sld_pop_init, 
      n_chr=input$sld_chromes, 
      n_gpchr=input$sld_genes_per_chrome, 
      h2=input$sld_herit_init
    )
  })
  
  output$plot_gen_gain <- renderPlot()
  
}
    
shinyApp(ui, server)
    



