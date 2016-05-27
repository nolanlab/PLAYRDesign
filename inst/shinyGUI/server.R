
reactiveGraph <- function (outputId)
{
      HTML(paste("<div id=\"", outputId, "\" class=\"shiny-reactiveGraph-output\"><svg /></div>", sep=""))
}

render_ui <- function(working.dir, ...){renderUI({
      fluidPage(
            fluidRow(
                  column(4,
                         fluidRow(
                               h1("PLAYRDesign"),
                               selectInput("est_file", "Select EST file", width = "100%", choices = ""),
                               selectInput("txdb_file", "Select Transcript database file", width = "100%", choices = ""),
                               selectInput("rep_db", "Select database of repetitive sequence", width = "100%", choices = ""),
                               selectInput("rna_db", "Select target transcript database", width = "100%", choices = ""),
                               selectInput("input_file", "Select input file", width = "100%", choices = c("", list.files(path = working.dir, pattern = "*.fasta$")))
                         ),
                         fluidRow(
                               column(4, numericInput("tm_min", "Tm Min", 61, min = 0, max = 100)),
                               column(4, numericInput("tm_opt", "Tm Opt", 63, min = 0, max = 100)),
                               column(4, numericInput("tm_max", "Tm Max", 65, min = 0, max = 100))
                         ),
                         fluidRow(
                               column(4,numericInput("len_min", "Oligo length Min", 20, min = 1, max = 100)),
                               column(4, numericInput("len_opt", "Oligo length Opt", 24, min = 1, max = 100)),
                               column(4, numericInput("len_max", "Oligo length Max", 26, min = 1, max = 100))
                         ),
                         fluidRow(
                               column(4,numericInput("product_min", "Product size Min", 50, min = 1, max = 1000)),
                               column(4, numericInput("product_max", "Product size Max", 60, min = 1, max = 1000))
                         ),
                         fluidRow(
                               column(4,numericInput("gc_min", "GC Min", 40, min = 1, max = 100)),
                               column(4, numericInput("gc_max", "GC Max", 70, min = 1, max = 100))
                         ),
                         fluidRow(
                               column(4, numericInput("oligos_to_report", "Number of oligos to report", 50))
                         ),
                         fluidRow(
                               actionButton("start_button", "Start analysis")
                         ),
                         fluidRow(
                               selectInput("selected_oligos", "Currently selected oligos", choices = c(""), multiple = T)
                         ),
                         fluidRow(
                               selectInput("PLAYR_system", "Select PLAYR system", choices =  unique(read.table(system.file("PLAYR_Systems.txt", package = "PLAYRDesign"), header = T, stringsAsFactors = F)$Name))
                         ),
                         fluidRow(
                               numericInput("start_playr_id", "Enter id for first oligo", 1)
                         ),
                         fluidRow(
                               actionButton("write_oligos", "Write oligos")
                         )
                         #fluidRow(
                         #      h4("Currently selected oligo"),
                         #      verbatimTextOutput("currently_selected_oligo")
                         #)
                  ),
                  column(8,
                         singleton(tags$head(tags$script(src = "d3.min.js"))),
                         singleton(tags$head(tags$script(src = "PLAYRDesign.js"))),
                         singleton(tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'PLAYRDesign.css'))),
                         reactiveGraph(outputId = "main_graph")
                  )
            )
      )
})}


read_options_file <- function(f_name)
{
      tab <- read.table(f_name, header = F, sep = "=", stringsAsFactors = F)
      ret <- as.list(tab$V2)
      names(ret) <- tab$V1
      return(ret)
}


do_est_analysis <- function(f_name, gr.est, txdb_file)
{
      
      id <- PLAYRDesign:::get_refseq_id_from_fasta(f_name)
      gr <- PLAYRDesign:::get_exons_for_transcript(id, txdb_file)
      v <- PLAYRDesign:::get_exons_skips(gr, gr.est)
      orig.seq <- PLAYRDesign:::get_seq_from_file(f_name)
      genomic.seq <- PLAYRDesign:::get_seq_from_ranges(gr)
      
      ret <- PLAYRDesign:::project_data_on_seq(orig.seq, genomic.seq, v)
      ret <- data.frame(pos = 1:length(orig.seq), exon_skip = ret)
      return(ret)
}


shinyServer(function(input, output, session)
{
      working.dir <- dirname(file.choose())
      output$PLAYRDesignUI <- render_ui(working.dir, input, output, session)
      playrdesign_opt <- read_options_file(file.path(working.dir, "playrdesign_conf.txt"))
      if(!is.null(playrdesign_opt$PLAYRDESIGN_DATA) && playrdesign_opt$PLAYRDESIGN_DATA != "")
      {
            updateSelectInput(session, 
                  "est_file", choices = c("", list.files(playrdesign_opt$PLAYRDESIGN_DATA, pattern = ".RData$")))
            updateSelectInput(session,
                  "txdb_file", choices = c("", list.files(playrdesign_opt$PLAYRDESIGN_DATA, pattern = ".sqlite$")))
      }
      updateSelectInput(session,
            "rep_db", choices = c("", list.files(playrdesign_opt$BLASTN_DB, pattern = ".fa$")))
      updateSelectInput(session,
            "rna_db", choices = c("", list.files(playrdesign_opt$BLASTN_DB, pattern = ".fa$")))
      
      
      cur_primer3_data <- NULL
      
      gr.est <- reactive({
            gr.est <- NULL
            if(!is.null(input$est_file) && input$est_file != "")
            {
                  print("Loading EST data...")
                  gr.est <- readRDS(file.path(playrdesign_opt$PLAYRDESIGN_DATA, input$est_file))
                  print("Done")
            }
            return(gr.est)
      })
      
      
      seq_data <- reactive({
            if(!is.null(input$input_file) && input$input_file != "")
            {
                  f_name <- paste(working.dir, input$input_file, sep = .Platform$file.sep)
                  print(sprintf("Running sequence analysis for %s", f_name))
                  blast_refseq <- PLAYRDesign:::run_blast_analysis_for_seq(f_name, db = input$rna_db, filter_same_gi = TRUE, playrdesign_opt)
                  blast_repbase <- PLAYRDesign:::run_blast_analysis_for_seq(f_name, db = input$rep_db, filter_same_gi = FALSE, playrdesign_opt)
                  blast <- data.frame(pos = blast_repbase$pos, blast1 = blast_refseq$percIdentity, blast2 = blast_repbase$percIdentity)
                  seq_char <- PLAYRDesign:::get_sequence_characteristics(f_name)
                  est.data <- gr.est()
                  est <- NULL
                  if(!is.null(est.data))
                        est <- do_est_analysis(f_name, est.data, file.path(playrdesign_opt$PLAYRDESIGN_DATA, input$txdb_file))
                  
                  return(list(blast = blast, seq_char = seq_char, est = est))
            }
            else
                  return(NULL)
      })      
      
      primer3_data <- reactive({
                  len <- c(input$len_min, input$len_opt, input$len_max)
                  tm <- c(input$tm_min, input$tm_opt, input$tm_max)
                  gc <- c(input$gc_min, input$gc_max)
                  product_size <- c(input$product_min, input$product_max)
                  f_name <- paste(working.dir, input$input_file, sep = .Platform$file.sep)
                  primer3 <- PLAYRDesign:::run_primer3(f_name, 
                              n = input$oligos_to_report, len = len, tm = tm, gc = gc, product_size = product_size, playrdesign_opt = playrdesign_opt)
                  return(primer3$tab_primers)
      })
      
      output$main_graph <- reactive({
            if(!is.null(input$start_button) && input$start_button != 0)
            {
                  primer3 <- primer3_data()
                  print(primer3)
                  cur_primer3_data <<- primer3
                  seq.data <- seq_data()
                  ret <- c(seq.data, list(primer3 = primer3))
                  updateSelectInput(session, "selected_oligos", selected = "", choices = "")
                  return(ret)
            }
      })
      
      output$currently_selected_oligo <- reactive({
            if(is.null(input$currently_selected_oligo))
                  return("None")
            else
                  return(input$currently_selected_oligo)
      })
      
      observe({
            if(!is.null(input$write_oligos) && input$write_oligos != 0)
            {
                  isolate({
                        tab <- cur_primer3_data
                        f_name <- paste(input$input_file, "playrdesign_out.txt", sep = ".")
                        gene_name <- gsub(".fasta", "", input$input_file)
                        PLAYRDesign:::write_selected_oligos(tab, input$selected_oligos, input$PLAYR_system, file.path(working.dir, f_name), input$start_playr_id, gene_name)
                  })
            }
      })
      
      observe({
            if(!is.null(input$click_select_oligo) && input$click_select_oligo != "")
            {
                  print(input$click_select_oligo)
                  isolate({
                        cur.sel <- input$selected_oligos
                        v <- c(cur.sel, input$click_select_oligo)
                        updateSelectInput(session, "selected_oligos", selected = v, choices = v)
                  })
            }
      })
})