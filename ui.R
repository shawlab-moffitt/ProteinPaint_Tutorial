####----Install and load packages----####

packages <- c("shiny","shinythemes","shinyjqui","dplyr","DT",
              "reshape2","tibble","readr","tidyr","shinycssloaders")

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))


####----User Data Input----####

#Parameter File
param_file <- 'Example1_Input_Files/AppParams_Rh18_SF3B1_Sudemycin.txt'


####----Read in files----####

params <- as.data.frame(read_delim(param_file,delim = '\t',col_names = c("Parameter","Value")))

ProjectName <- params[which(params[,1] == "ProjectName"),2]
meta_file <- params[which(params[,1] == "MetaFile"),2]
gene_coord_file <- params[which(params[,1] == "CoordFile"),2]
GenBuild <- params[which(params[,1] == "GenBuild"),2]
Stranded <- params[which(params[,1] == "TwoStranded"),2]


meta_withFiles <- as.data.frame(read_delim(meta_file, delim = '\t'))

if (Stranded == "TRUE") {
    meta_file_rev <- meta_withFiles[which(grepl(".rev",meta_withFiles$BigWig) == TRUE),]
    meta_file_fwd <- meta_withFiles[which(grepl(".fwd",meta_withFiles$BigWig) == TRUE),]
    meta <- meta_withFiles[,-c(2:4), drop = F]
    meta <- unique(meta)
}
if (Stranded == "FALSE") {
    meta_withFiles <- as.data.frame(read_delim(meta_file, delim = '\t'))
    meta <- meta_withFiles[,-c(2:4), drop = F]
}

coord <- as.data.frame(read_delim(gene_coord_file, delim = '\t', col_names = c("Gene","Coordinates")))


####----UI----####

# Set Shiny theme
shinytheme("sandstone")


ui <- 
    navbarPage(paste("{",ProjectName,"Protein Paint Analysis Tool }"),
               tabPanel("Overview",
                        fluidPage(
                            mainPanel(
                                p("")
                            )
                        )
               ),
               tabPanel("Three Sample Comparison",
                        fluidPage(
                            title = "Three Sample Comparison",
                            sidebarLayout(
                                sidebarPanel(
                                    tabsetPanel(
                                        id = "CompTwoTabs",
                                        tabPanel("Sample Selection",
                                                 p(),
                                                 uiOutput("SelectCondition1"),
                                                 h3("Select Sample 1"),
                                                 div(DT::dataTableOutput("SampleOptTable1"), style = "font-size:10px; height:300px"),
                                                 h3("Select Sample 2"),
                                                 div(DT::dataTableOutput("SampleOptTable2"), style = "font-size:10px; height:300px"),
                                                 h3("Select Sample 3"),
                                                 div(DT::dataTableOutput("SampleOptTable3"), style = "font-size:10px; height:300px")
                                        ),
                                        tabPanel("Position Selection",
                                                 p(),
                                                 uiOutput("CoordTableOut"),
                                                 uiOutput("coordTable1output"),
                                                 checkboxInput("UserCoordButton","User Upload of Gene-Coordinate File"),
                                                 value = 22
                                        ),
                                        tabPanel("File Type",
                                                 p(),
                                                 #checkboxGroupInput("filetype","Select File Type to View",choices = c("BAM","BigWig","Junction")),
                                                 #radioButtons("filetype","Select File Type to View", choices = c("BAM","BigWig","Junction")),
                                                 radioButtons("filetype","Select File Types to View", choices = c("Junction & BigWig","Junction & BAM")),
                                                 value = 33
                                        )
                                    )
                                ),
                                mainPanel(
                                    id = "TwoSampComp",
                                    uiOutput("GenomeProtPaint")
                                )
                            )
                        )
               )
    )






















