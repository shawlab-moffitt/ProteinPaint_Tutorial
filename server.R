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



####----Server----####

server <- function(input, output, session) {
    
    # Render Select condition for sample 1
    output$SelectCondition1 <- renderUI({
        
        if ((ncol(meta) > 2) == TRUE) {
            
            choiceColumns <- colnames(meta)
            choiceColumns <- choiceColumns[-1]
            selectInput("CondOne","Select Descriptor(s)",
                        choices = choiceColumns, multiple = T)
            
        }
        
    })
    
    # Sample Selection table 1
    output$SampleOptTable1 <- renderDataTable({
        
        if ((ncol(meta) > 2) == TRUE) {
            
            sample_list <- meta[,1]
            subMeta <- meta[,input$CondOne, drop = F]
            subMeta$Samples <- sample_list
            subMeta <- subMeta %>%
                relocate(Samples)
            colnames(subMeta)[1] <- colnames(meta)[1]
            DT::datatable(subMeta,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 1),
                          rownames = F)
            
        }
        else if ((ncol(meta) <= 2) == TRUE) {
            
            meta_tb <- as_tibble(meta)
            DT::datatable(meta_tb,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 1),
                          rownames = F)
            
        }
        
    })
    
    # Sample Selection table 2
    output$SampleOptTable2 <- renderDataTable({
        
        if ((ncol(meta) > 2) == TRUE) {
            
            sample_list <- meta[,1]
            subMeta <- meta[,input$CondOne, drop = F]
            subMeta$Samples <- sample_list
            subMeta <- subMeta %>%
                relocate(Samples)
            colnames(subMeta)[1] <- colnames(meta)[1]
            DT::datatable(subMeta,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 2),
                          rownames = F)
            
        }
        else if ((ncol(meta) <= 2) == TRUE) {
            
            meta_tb <- as_tibble(meta)
            DT::datatable(meta_tb,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 2),
                          rownames = F)
            
        }
        
    })
    
    # Sample Selection table 3
    output$SampleOptTable3 <- renderDataTable({
        
        if ((ncol(meta) > 2) == TRUE) {
            
            sample_list <- meta[,1]
            subMeta <- meta[,input$CondOne, drop = F]
            subMeta$Samples <- sample_list
            subMeta <- subMeta %>%
                relocate(Samples)
            colnames(subMeta)[1] <- colnames(meta)[1]
            DT::datatable(subMeta,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 3),
                          rownames = F)
            
        }
        else if ((ncol(meta) <= 2) == TRUE) {
            
            meta_tb <- as_tibble(meta)
            DT::datatable(meta_tb,
                          options = list(lengthMenu = c(5,10, 20, 100, 1000),
                                         pageLength = 5,
                                         scrollX = T),
                          selection=list(mode = 'single', selected = 3),
                          rownames = F)
            
        }
        
    })
    
    output$CoordTableOut <- renderUI({
        
        if (input$UserCoordButton == TRUE) {
            
            fileInput("UserCoordFile", "Gene-Coordinate File (no header)", accept = c(".txt",".tsv",".csv"))
            
        }
        
    })
    
    coord.u <- reactive({
        
        cd.u <- input$UserCoordFile
        ext <- tools::file_ext(cd.u$datapath)
        req(cd.u)
        validate(need(ext == c("txt","tsv","csv"), "Please upload .txt, .tsv, or .csv file (no header)"))
        if (ext == "csv") {
            
            as.data.frame(read_delim(cd.u$datapath, delim = ',', col_names = c("Gene","Coordinates")))
            
        }
        else {
            
            as.data.frame(read_delim(cd.u$datapath, delim = '\t', col_names = c("Gene","Coordinates")))
            
        }
        
    })
    
    output$coordTable1output <- renderUI({
        
        if (input$UserCoordButton == FALSE) {
            
            div(DT::dataTableOutput("coordTable1"), style = "font-size:10px; height:500px; overflow-X: scroll")
            
        }
        else if (input$UserCoordButton == TRUE) {
            
            req(input$UserCoordFile)
            div(DT::dataTableOutput("coordTable1"), style = "font-size:10px; height:500px; overflow-X: scroll")
            
        }
        
    })
    
    # Coordinate position table
    output$coordTable1 <- renderDataTable({
        
        if (input$UserCoordButton == FALSE) {
            
            coord_tb <- as_tibble(coord)
            DT::datatable(coord_tb, options = list(lengthMenu = c(10, 20, 100, 1000), pageLength = 10),
                          selection=list(mode = 'single', selected = 1),
                          rownames = F)
            
        }
        else if (input$UserCoordButton == TRUE) {
            
            req(input$UserCoordFile)
            user_coord <- coord.u()
            DT::datatable(user_coord, options = list(lengthMenu = c(10, 20, 100, 1000), pageLength = 10),
                          selection=list(mode = 'single', selected = 1),
                          rownames = F)
            
        }
        
    })
    
    
    
    
    output$GenomeProtPaint <- renderUI({
        
        #req(input$coordTable1_rows_selected)
        # Sample Names
        samp_name1 <- as.character(meta[input$SampleOptTable1_rows_selected,1])
        samp_name2 <- as.character(meta[input$SampleOptTable2_rows_selected,1])
        samp_name3 <- as.character(meta[input$SampleOptTable3_rows_selected,1])
        # Position
        if (length(as.character(input$coordTable1_rows_selected)) > 0) {
            position = as.character(coord[input$coordTable1_rows_selected,2])
            highlighted_position = as.character(coord[input$coordTable1_rows_selected,3])
        }
        else {
            position = as.character(coord[1,2])
            highlighted_position = as.character(coord[1,3])
        }
        # File Type
        file_type <- input$filetype
        
        if (Stranded == "TRUE") {
            if (file_type == "Junction & BAM") {
        
                Strand <- "_ONE"
                
                # Junction Files
                samp1_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name1),"Junction"]
                samp2_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name2),"Junction"]
                samp3_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name3),"Junction"]
                
                # BAM Files
                samp1_2 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name1),"BAM"]
                samp2_2 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name2),"BAM"]
                samp3_2 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name3),"BAM"]
                
                file_type1 <- "junction"
                file_type2 <- "bam"
                
                urllink = paste("http://biostools/4472414/Shiny/Dev/ProteinPaintTest/HTML_Files/", # This should be a URL path to your HTML file hosted on a server
                                "ProtPaint_3samp_Shiny_",GenBuild,Strand,".html?", #HTML Script
                                "sample_name1_1=",samp1_1,        #Sample 1 path - FT 1
                                "&sample_name2_1=",samp2_1,       #Sample 2 path - FT 1
                                "&sample_name3_1=",samp3_1,       #Sample 3 path - FT 1
                                "&sample_name1_2=",samp1_2,        #Sample 1 path - FT 2
                                "&sample_name2_2=",samp2_2,       #Sample 2 path - FT 2
                                "&sample_name3_2=",samp3_2,       #Sample 3 path - FT 2
                                "&position=", position,       #Position Chosen
                                "&highlight=", highlighted_position,       #Position Chosen
                                "&file_type1=",file_type1,    #File type
                                "&file_type2=",file_type2,    #File type
                                "&samp1=",samp_name1,            #Sample 1 Name
                                "&samp2=",samp_name2,            #Sample 2 Name
                                "&samp3=",samp_name3,            #Sample 3 Name
                                sep="")
                my_test <- tags$iframe(src=urllink, frameBorder="0", height=2000, width=2035)
                print(my_test)
                
            }
            else if (file_type == "Junction & BigWig") {
                Strand <- "_TWO"
                
                # Junction Files
                samp1_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name1),"Junction"]
                samp2_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name2),"Junction"]
                samp3_1 <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name3),"Junction"]
                
                # Junction Files
                #samp1_1f <- meta_file_fwd[which(meta_file_fwd$SampleName == samp_name1),"Junction"]
                #samp2_1f <- meta_file_fwd[which(meta_file_fwd$SampleName == samp_name2),"Junction"]
                #samp3_1f <- meta_file_fwd[which(meta_file_fwd$SampleName == samp_name3),"Junction"]
                #samp1_1r <- meta_file_rev[which(meta_file_rev$SampleName == samp_name1),"Junction"]
                #samp2_1r <- meta_file_rev[which(meta_file_rev$SampleName == samp_name2),"Junction"]
                #samp3_1r <- meta_file_rev[which(meta_file_rev$SampleName == samp_name3),"Junction"]
                
                # BW Files
                samp1_2f <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name1),"BigWig"]
                samp2_2f <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name2),"BigWig"]
                samp3_2f <- meta_file_fwd[which(meta_file_fwd[,1] == samp_name3),"BigWig"]
                samp1_2r <- meta_file_rev[which(meta_file_rev[,1] == samp_name1),"BigWig"]
                samp2_2r <- meta_file_rev[which(meta_file_rev[,1] == samp_name2),"BigWig"]
                samp3_2r <- meta_file_rev[which(meta_file_rev[,1] == samp_name3),"BigWig"]
                
                file_type1 <- "junction"
                file_type2 <- "bigwigstranded"
                
                urllink = paste("http://biostools/4472414/Shiny/Dev/ProteinPaintTest/HTML_Files/", # This should be a URL path to your HTML file hosted on a server
                                "ProtPaint_3samp_Shiny_",GenBuild,Strand,".html?", #HTML Script
                                "sample_name1_1f=",samp1_1,        #Sample 1 path - FT 1
                                "&sample_name2_1f=",samp2_1,       #Sample 2 path - FT 1
                                "&sample_name3_1f=",samp3_1,       #Sample 3 path - FT 1
                                "&sample_name1_2f=",samp1_2f,       #Sample 1 path - FT 2
                                "&sample_name2_2f=",samp2_2f,       #Sample 2 path - FT 2
                                "&sample_name3_2f=",samp3_2f,       #Sample 3 path - FT 2
                                "&sample_name1_2r=",samp1_2r,       #Sample 1 path - FT 2
                                "&sample_name2_2r=",samp2_2r,       #Sample 2 path - FT 2
                                "&sample_name3_2r=",samp3_2r,       #Sample 3 path - FT 2
                                "&position=", position,             #Position Chosen
                                "&highlight=", highlighted_position,#Position Chosen
                                "&file_type1=",file_type1,          #File type
                                "&file_type2=",file_type2,          #File type
                                "&samp1=",samp_name1,               #Sample 1 Name
                                "&samp2=",samp_name2,               #Sample 2 Name
                                "&samp3=",samp_name3,               #Sample 3 Name
                                sep="")
                my_test <- tags$iframe(src=urllink, frameBorder="0", height=2000, width=2035)
                print(my_test)
                
                
            }
            
        }
        else if (Stranded == "FALSE") {
            Strand <- "_ONE"
            if (file_type == "Junction & BAM") {
                
                # Junction Files
                samp1_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name1),"Junction"]
                samp2_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name2),"Junction"]
                samp3_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name3),"Junction"]
                
                # BAM Files
                samp1_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name1),"BAM"]
                samp2_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name2),"BAM"]
                samp3_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name3),"BAM"]
                
                file_type1 <- "junction"
                file_type2 <- "bam"
                
            }
            else if (file_type == "Junction & BigWig") {
                
                # Junction Files
                samp1_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name1),"Junction"]
                samp2_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name2),"Junction"]
                samp3_1 <- meta_withFiles[which(meta_withFiles[,1] == samp_name3),"Junction"]
                # BW Files
                samp1_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name1),"BigWig"]
                samp2_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name2),"BigWig"]
                samp3_2 <- meta_withFiles[which(meta_withFiles[,1] == samp_name3),"BigWig"]
                
                file_type1 <- "junction"
                file_type2 <- "bigwig"
                
            }
            
            urllink = paste("http://biostools/4472414/Shiny/Dev/ProteinPaintTest/HTML_Files/", # This should be a URL path to your HTML file hosted on a server
                            "ProtPaint_3samp_Shiny_",GenBuild,Strand,".html?", #HTML Script
                            "sample_name1_1=",samp1_1,        #Sample 1 path - FT 1
                            "&sample_name2_1=",samp2_1,       #Sample 2 path - FT 1
                            "&sample_name3_1=",samp3_1,       #Sample 3 path - FT 1
                            "&sample_name1_2=",samp1_2,        #Sample 1 path - FT 2
                            "&sample_name2_2=",samp2_2,       #Sample 2 path - FT 2
                            "&sample_name3_2=",samp3_2,       #Sample 3 path - FT 2
                            "&position=", position,       #Position Chosen
                            "&highlight=", highlighted_position,       #Position Chosen
                            "&file_type1=",file_type1,    #File type
                            "&file_type2=",file_type2,    #File type
                            "&samp1=",samp_name1,            #Sample 1 Name
                            "&samp2=",samp_name2,            #Sample 2 Name
                            "&samp3=",samp_name3,            #Sample 3 Name
                            sep="")
            my_test <- tags$iframe(src=urllink, frameBorder="0", height=2000, width=2035)
            print(my_test)
                
        }
        
    })
    
}



















