
# Libraries and dependencies ----------------------------------------------


library(shiny)
library(GOexpress)


# Load required objects and set default values ----------------------------


# Instead of the entire GOexpress result and annotations,
# load pre-extracted list of genes and GO terms:

# Gene names represented in the data set
genes_choices = readRDS(file = 'data/external_gene_names.rds')

# GO terms represented in the data set
go_choices = readRDS(file = 'data/go_choices.rds')

# Ensembl gene identifiers
gene_ids = readRDS(file = 'data/ensemblIDs.rds')

# Time points
hours_choices = readRDS(file = 'data/timepoints.rds')

# Animal identifiers
animals_choices = readRDS(file = 'data/animalIDs.rds')


# Populate user controls and set default values ---------------------------


# Default gene displayed (SAA3)
default_geneId = 'ENSBTAG00000022396'

# Default gene displayed
default_gene = 'SAA3'

# Default GO term displayed in heatmap
# (positive regulation of interferon-alpha production)
default_go = 'GO:0032727'

# Default time point(s) displayed in the heatmap
hours_heatmap <- c('48 hrs')

# Default row label size in heatmap
cexRow_choice = 1.5

# Default right margin size in heatmap
right_margin = 10

# Prepare the individual selectable infection choices
infection_choices = list(
  'Control'='Control',
  'M. tuberculosis'='M. tuberculosis',
  'M. bovis'='M. bovis')

# Default infection(s) displayed in all plots
infection_selected = as.character(unlist(infection_choices))

# maximum filter value allowed for minimal count of total annotated
# genes
max.GO.total = 1E3

# Default filter for the minimal count of annotated genes
GO.filter_choice = 6


# Define the page components and layout -----------------------------------


shinyUI(fluidPage(
  
  titlePanel('AlvMac full app (updated 08/04/2015)'),
  
  navlistPanel(
    widths = c(2, 10),
    
    'Genes',
    
    tabPanel(
      'Expression profiles (Gene name)',
      h3('Expression profiles by gene name'),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = 'external_gene_name',
            label = 'Gene name:',
            choices = genes_choices,
            selected = default_gene),
          
          checkboxGroupInput(
            inputId = 'animals_symbol',
            label = 'Animal IDs:',
            choices = animals_choices,
            selected = animals_choices,
            inline = TRUE),
          
          checkboxGroupInput(
            inputId = 'infection_symbol',
            label = 'Infection:',
            choices = infection_choices,
            selected = infection_selected),
          
          checkboxGroupInput(
            inputId = 'hours_symbol',
            label = 'Hours post-infection:',
            choices = hours_choices,
            selected = hours_choices[-1],
            inline = TRUE),
          
          sliderInput(
            inputId = 'linesize',
            label = 'Line size',
            min = 0,
            max = 2,
            value = 1.5,
            step = 0.25
          ),
          
          numericInput(
            inputId = 'index',
            label = 'Plot index: (0 for all plots)',
            value = 0,
            min = 0,
            step = 1
          )
          
        ), # end of sidebarPanel
        
        mainPanel(
          tabsetPanel(
            type = 'pills',
            tabPanel(
              'Sample series',
              plotOutput(
                'exprProfilesSymbol',
                width = '100%', height = '600px'
              )
            ), 
            tabPanel(
              'Sample groups',
              plotOutput(
                'exprPlotSymbol',
                width = '100%', height = '600px'
              )
            )
            
          )
          
        ) #  end of mainPanel
        
      ) # end of sidebarLayout
      
      
    ), # end of tabPanel
    
    tabPanel(
      'Expression profiles (Ensembl ID)',
      h3('Expression profiles by Ensembl gene identifier'),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = 'ensembl_gene_id',
            label = 'Gene name:',
            choices = gene_ids,
            selected = default_geneId),
          
          checkboxGroupInput(
            inputId = 'animals',
            label = 'Animal IDs:',
            choices = animals_choices,
            selected = animals_choices,
            inline = TRUE),
          
          checkboxGroupInput(
            inputId = 'infection',
            label = 'Infection:',
            choices = infection_choices,
            selected = infection_selected),
          
          checkboxGroupInput(
            inputId = 'hours',
            label = 'Hours post-infection:',
            choices = hours_choices,
            selected = hours_choices[-1],
            inline = TRUE)
          
        ), # end of sidebarPanel
        
        mainPanel(
          tabsetPanel(
            type = 'pills',
            tabPanel(
              'Sample series',
              plotOutput(
                'exprProfiles',
                width = '100%', height = '600px'
              )
            ), 
            tabPanel(
              'Sample groups',
              plotOutput(
                'exprPlot',
                width = '100%', height = '600px'
              )
            )
            
          )
          
        ) #  end of mainPanel
        
      ) # end of sidebarLayout
      
      
    ), # end of tabPanel
    
    tabPanel(
      'Scoring table',
      h3('Scoring table'),
      dataTableOutput('genesScore')
    ),
    
    'Gene ontologies',
    
    tabPanel(
      'Heatmap',
      h3('Heatmap'),
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = 'go_id',
            label = 'GO ID:',
            choices = go_choices,
            selected = default_go
          ),
          
          checkboxGroupInput(
            inputId = 'hours.GO',
            label = 'Hours post-infection:',
            choices = hours_choices,
            selected = hours_heatmap,
            inline = TRUE),
          
          checkboxGroupInput(
            inputId = 'infection.GO',
            label = 'Infection:',
            choices = infection_choices,
            selected = infection_selected),
          
          checkboxGroupInput(
            inputId = 'animal.GO',
            label = 'Animals:',
            choices = animals_choices,
            selected = animals_choices,
            inline = TRUE),
          
          sliderInput(
            inputId = 'cexRow.GO',
            label = 'Row label size:',
            min = 0.5,
            max = 2,
            value = cexRow_choice,
            step = 0.1
          ),
          
          sliderInput(
            inputId = 'margins.heatmap.bottom',
            label = 'Bottom margin:',
            min = 3,
            max = 30,
            value = 7,
            step = 1
          ),
          
          sliderInput(
            inputId = 'margins.heatmap.right',
            label = 'Right margin:',
            min = 3,
            max = 30,
            value = right_margin,
            step = 1
          )
          
        ),
        mainPanel(
          plotOutput(
            'heatmap',
            width = '100%', height = '700px'
          )
        )
      )
    ),
    
    tabPanel(
      'Scoring table',
      h3('Scoring table'),
      fluidRow(
        column(width = 1,
               numericInput(
                 inputId = 'min.total',
                 label = 'Min. total:',
                 min = 0,
                 max = max.GO.total,
                 value = GO.filter_choice
               )
        )
      ),
      dataTableOutput('GOscores')
    ),
    
    '-----',
    
    tabPanel(
      'Samples info',
      h3('Sample phenotypic information'),
      dataTableOutput('Adataframe')
    )
    
  )
))
