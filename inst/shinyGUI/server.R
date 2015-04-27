
reactiveGraph <- function (outputId)
{
      HTML(paste("<div id=\"", outputId, "\" class=\"shiny-reactiveGraph-output\"><svg /></div>", sep=""))
}

render_ui <- function(working.dir, ...){renderUI({
      fluidPage(
            fluidRow(
                  column(4,
                         fluidRow(
                               selectInput("input_file", "Select input file", choices = c("", list.files(path = working.dir, pattern = "*.fasta$")))
                         ),
                         fluidRow(
                               column(4,numericInput("tm_min", "Tm Min", 61, min = 0, max = 100)),
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
                               selectInput("selected_oligos", "Select oligos", choices = c(""), multiple = T)
                         ),
                         fluidRow(
                               selectInput("PLAYR_system", "Select PLAYR system", choices =  unique(read.table(system.file("PLAYR_Systems.txt", package = "PLAYRDesign"), header = T, stringsAsFactors = F)$Name))
                         ),
                         fluidRow(
                               numericInput("start_playr_id", "Enter id for first oligo", 1)
                         ),
                         fluidRow(
                               actionButton("write_oligos", "Write oligos")
                         ),
                         fluidRow(
                               h4("Currently selected oligo"),
                               verbatimTextOutput("currently_selected_oligo")
                         )
                  ),
                  column(8,
                         singleton(tags$head(tags$script(src = "http://d3js.org/d3.v3.min.js"))),
                         singleton(tags$head(tags$script(src = "RLADesign.js"))),
                         singleton(tags$head(tags$link(rel = 'stylesheet', type = 'text/css', href = 'RLADesign.css'))),
                         reactiveGraph(outputId = "main_graph")
                  )
            )
      )
})}


get_data_for_graph <- function(f_name, to_report, len, tm, gc, product_size)
{
      print(f_name)
      blast_refseq <- RLADesign:::run_blast_analysis_for_seq(f_name, db = "rna_human_high_qual.fa", filter_same_gi = TRUE)
      blast_repbase <- RLADesign:::run_blast_analysis_for_seq(f_name, db = "repbase.fa", filter_same_gi = FALSE)
      blast <- data.frame(pos = blast_repbase$pos, blast1 = blast_refseq$percIdentity, blast2 = blast_repbase$percIdentity)
      primer3 <- RLADesign:::run_primer3(f_name, 
                        n = to_report, len = len, tm = tm, gc = gc, product_size = product_size)
      seq_char <- RLADesign:::get_sequence_characteristics(f_name)

      return(list(blast = blast, primer3 = primer3$tab_primers, seq_char = seq_char))
      
      
      #return(list(pos = df$pos, eVal = df$eVal, bitScore = df$bitScore, percIdentity = df$percIdentity))
}


do_est_analysis <- function(f_name, gr.est)
{
      
      id <- RLADesign:::get_refseq_id_from_fasta(f_name)
      gr <- RLADesign:::get_exons_for_transcript(id)
      v <- RLADesign:::get_exons_skips(gr, gr.est)
      orig.seq <- RLADesign:::get_seq_from_file(f_name)
      genomic.seq <- RLADesign:::get_seq_from_ranges(gr)
      
      ret <- RLADesign:::project_data_on_seq(orig.seq, genomic.seq, v)
      ret <- data.frame(pos = 1:length(orig.seq), exon_skip = ret)
      return(ret)
}


shinyServer(function(input, output, session)
{
      working.dir <- dirname(file.choose())
      output$RLADesignUI <- render_ui(working.dir, input, output, session)
      
      print("Loading EST data...")
      gr.est <- readRDS(system.file("spliced_est_hg19.RData", package = "RLADesign"))
      print("Done")
      
      cur_primer3_data <- NULL
      
      seq_data <- reactive({
            if(!is.null(input$input_file) && input$input_file != "")
            {
                  f_name <- paste(working.dir, input$input_file, sep = .Platform$file.sep)
                  print(sprintf("Running sequence analysis for %s", f_name))
                  blast_refseq <- PLAYRDesign:::run_blast_analysis_for_seq(f_name, db = "rna_human_high_qual.fa", filter_same_gi = TRUE)
                  blast_repbase <- PLAYRDesign:::run_blast_analysis_for_seq(f_name, db = "repbase.fa", filter_same_gi = FALSE)
                  blast <- data.frame(pos = blast_repbase$pos, blast1 = blast_refseq$percIdentity, blast2 = blast_repbase$percIdentity)
                  seq_char <- PLAYRDesign:::get_sequence_characteristics(f_name)
                  est <- do_est_analysis(f_name, gr.est)
                  return(list(blast = blast, seq_char = seq_char, est = est))
            }
            else
                  return(NULL)
      })      
      
      primer3_data <- reactive({
            if(!is.null(input$start_button) && input$start_button != 0)
                  isolate({
                        len <- c(input$len_min, input$len_opt, input$len_max)
                        tm <- c(input$tm_min, input$tm_opt, input$tm_max)
                        gc <- c(input$gc_min, input$gc_max)
                        product_size <- c(input$product_min, input$product_max)
                        f_name <- paste(working.dir, input$input_file, sep = .Platform$file.sep)
                        primer3 <- PLAYRDesign:::run_primer3(f_name, 
                                                           n = input$oligos_to_report, len = len, tm = tm, gc = gc, product_size = product_size)
                        return(primer3$tab_primers)
                  })
      })
      
      output$main_graph <- reactive({
            temp <- primer3_data()
            print(temp)
            if(!is.null(temp))
            {
                  seq.data <- seq_data()
                  if(!is.null(cur_primer3_data))
                  {
                        a <- max(cur_primer3_data$pair_row_id)
                        temp$pair_row_id <- temp$pair_row_id + a + 1
                        a <- max(cur_primer3_data$id)
                        temp$id <- temp$id + a + 1
                        a <- max(cur_primer3_data$unique_id)
                        temp$unique_id <- temp$unique_id + a + 1
                        temp <- rbind(cur_primer3_data, temp)
                  }
                  ret <- c(seq.data, list(primer3 = temp))
                  cur_primer3_data <<- temp
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
                        f_name <- paste(input$input_file, "playr_design_out.txt", sep = ".")
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