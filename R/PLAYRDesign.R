#' @import plyr
#' @import zoo
#' @import reshape2
#' @import shiny
#' @import BSgenome
#' @import GenomicRanges
#' @import GenomicFeatures
#' @import AnnotationDbi
#' @import AnnotationFuncs
#' @import IRanges
#' @import BSgenome.Hsapiens.UCSC.hg19
#' @import org.Hs.eg.db
#' @import XML
#' @import Biostrings



#' @export
PLAYRDesign.run <- function()
{
      runApp(appDir = file.path(system.file(package = "PLAYRDesign"), "shinyGUI"), launch.browser=T)
}


my_load <- function(f_name)
{
      con <- file(f_name, "rb")
      retval <- unserialize(con)
      close(con)
      return(retval)
}



load_est <- function(f_name)
{
      tab <- my_load(f_name)
      
      ret <- GRanges(seqnames = Rle(tab$tName), ranges = IRanges(start = tab$tStart, end = tab$tEnd), 
                    strand = Rle(tab$strand))
      
      mcols(ret) <- tab[, c("blockSizes", "tStarts")]
      names(ret) <- tab$qName
      return(ret)
      
      
      
      
}

load_hexeventdb <- function(f_name) 
{
      ret <- read.table(f_name, header = T, sep = "\t", quote = "", stringsAsFactors = F, comment.char = "")
      ret <- ret[, c("chromo", "strand", "start", "end", "count", "constitLevel", "OnlyESTexons", "genename")]
}




sort_primer_pairs <- function(tab)
{
      limits <- ddply(tab, ~id, function(x) {return(data.frame(start = min(x$start), end = max(x$start + x$len)))})
      ir <- IRanges(limits$start, limits$end)
      limits <- data.frame(id = limits$id, pair_row_id = disjointBins(ir) - 1)
      res <- merge(tab, limits)
      return(res)
}

select_primer3_output_lines <- function(s)
{
      v <- sapply(s, function(x) {return(x[[1]])})
      v <- strsplit(v, "_")
      i <- sapply(v, function(x) {length(x) > 1 && x[[2]] %in% c("PAIR", "LEFT", "RIGHT") && grepl("^[[:digit:]]*$", x[3])})
      return(s[i])
}

revcomp <- function(s)
{
      x <- DNAStringSet(s)
      return(as.character(reverseComplement(x)))
      
}


write_selected_oligos <- function(tab, selected_oligos, playr_system, f_name, start_id, gene_name)
{
      playr.systems <- read.table(system.file("PLAYR_Systems.txt", package = "PLAYRDesign"), 
                                header = T, sep = "\t", quote = "", stringsAsFactors = F)
      playr.systems <- playr.systems[playr.systems$Name == playr_system,]
      
      ret <- NULL      
      id <- start_id
      for(i in selected_oligos)
      {
            temp <- strsplit(i, "_")[[1]]
            m <- NULL
            if(length(temp) == 1)
            {
                  m <- tab[tab$id == i,]
            }
            else
            {
                  m <- tab[tab$unique_id %in% temp, ]
                  if(m[1, "start"] < m[2, "start"])
                        v <- c("LEFT", "RIGHT")
                  else
                        v <- c("RIGHT", "LEFT")
                  m$type <- v
            }
            R <- m$type == "RIGHT"
            L <- m$type == "LEFT"
            
            ret <- rbind(ret, data.frame(Name = sprintf("%s_%d_%s", gene_name, id, playr_system), 
                                         Sequence = paste(revcomp(m[L, "sequence"]), playr.systems[playr.systems$Type == "PLAYR2", "Sequence"], sep = ""),  m[L, c("tm", "start", "len")]))
            id <- id + 1
            
            #When parsing primer3 we reverse complement the right primer. Therefore it needs
            #to be reverse complemented here too.
            ret <- rbind(ret, data.frame(Name = sprintf("%s_%d_%s", gene_name, id, playr_system), 
                                         Sequence = paste(revcomp(m[R, "sequence"]), playr.systems[playr.systems$Type == "PLAYR1", "Sequence"], sep = ""), m[R, c("tm", "start", "len")]))
            id <- id + 1      
      }
      print(ret)
      ret <- data.frame(ret)
      names(ret)[1:2] <- c("Name", "Sequence")
      write.table(ret, f_name, row.names = F, quote = F, sep = "\t")
}

write_selected_oligos_old <- function(tab, playr_system, f_name, start_id, gene_name)
{
      playr.systems <- read.table(system.file("PLAYR_Systems.txt", package = "PLAYRDesign"), 
                                header = T, sep = "\t", quote = "", stringsAsFactors = F)
      playr.systems <- playr.systems[playr.systems$Name == playr_system,]
      
      x <- DNAStringSet(tab[tab$type == "LEFT", "sequence"])
      x <- as.character(reverseComplement(x))
      tab[tab$type == "LEFT", "sequence"] <- x
      ret <- NULL      
      id <- start_id
      for(i in unique(tab$id))
      {
            m <- tab[tab$id == i,]
            ret <- rbind(ret, c(sprintf("%s_%d_%s", gene_name, id, playr_system), 
                                paste(m[m$type == "RIGHT", "sequence"], playr.systems[playr.systems$Type == "PLAYR1", "Sequence"], sep = "")))
            id <- id + 1
            ret <- rbind(ret, c(sprintf("%s_%d_%s", gene_name, id, playr_system), 
                                paste(m[m$type == "LEFT", "sequence"], playr.systems[playr.systems$Type == "PLAYR2", "Sequence"], sep = "")))
            id <- id + 1
      }
      ret <- data.frame(ret)
      names(ret) <- c("Name", "Sequence")
      write.table(ret, f_name, col.names = T, row.names = F, quote = F, sep = "\t")
}



parse_primer3_output <- function(f_name)
{
      v <- readLines(f_name)
      s <- strsplit(v, "=")
      s <- select_primer3_output_lines(s)
      s <- t(sapply(s, c, simplify = T))
      
      temp <- strsplit(s[,1], "_")
      
      ret <- sapply(temp, function(x) 
            {
                  pair <- x[[3]]
                  type <- x[[2]]
                  key <- ""
                  if(length(x) > 3)
                        key <- paste(x[4:length(x)], collapse = "_")
                  else
                        key <- "POS"
                  return(c(pair, type, key))
            
      })
      ret <- t(ret)
      ret <- data.frame(ret, s[,2], stringsAsFactors = F)
      names(ret) <- c("id", "type", "key", "val")
      tab_pairs <- ret[ret$type == "PAIR",]
      tab_primers <- ret[ret$type != "PAIR",]
      ret <- list(tab_primers = tab_primers, tab_pairs = tab_pairs)
      ret <- lapply(ret, function(x) {dcast(x, id+type~key, value.var = "val")})
      ret <- lapply(ret, data_frame_factor_to_char)
      
      temp <- ret$tab_primers
      v <- strsplit(temp$POS, ",")
      v <- t(sapply(v,c))
      colnames(v) <- c("START", "LEN")
      temp$POS <- NULL
      temp <- data.frame(temp, v, stringsAsFactors = F)
      ret$tab_primers <- temp
      
      ret <- lapply(ret, function(x)
            {
                  i <- !(names(x) %in% c("type", "SEQUENCE"))
                  x[i] <- lapply(x[i], as.numeric)
                  x <- data.frame(x)
                  names(x) <- tolower(names(x))
                  return(x)
            })
      ret$tab_primers$start <- ret$tab_primers$start + 1
      #Primer3 reports the right primer as the reverse complement, so the first base is actually the last on the mRNA
      ret$tab_primers <- reverse_right_primers(ret$tab_primers)
      
      ret$tab_primers <- sort_primer_pairs(ret$tab_primers)
      ret$tab_primers <- cbind(ret$tab_primers, unique_id = 1:nrow(ret$tab_primers))
      return(ret)
}

reverse_right_primers <- function(tab)
{
      x <- tab$type == "RIGHT"
      tab[x, "sequence"] <- revcomp(tab[x, "sequence"])
      tab[x, "start"] <- tab[x, "start"] - tab[x, "len"] + 1
      
      return(tab)
}


get_sequence_characteristics <- function(f_name)
{
      s <- readDNAStringSet(f_name)
      s <- as.character(s, use.names = F)
      s <- strsplit(s, "")[[1]]
      gc <- rollapply(s, get_gc, width = 20)
      tm <- rollapply(s, get_tm, width = 20)
      return(data.frame(pos = 1:length(tm), gc, tm))
}

get_gc <- function(s)
{
      tot <- sum(s %in% c("G", "C"))
      return(tot / length(s))
}

get_tm <- function(s)
{
      tot <- sum(s %in% c("G", "C"))
      return(64.9 + (41 *(tot - 16.4)) / length(s))
      
}


      
      
run_primer3 <- function(f_name, n, len, tm, gc, product_size)
{
      opt <- list()
      v <- c("MIN", "OPT", "MAX")
      s <- paste("PRIMER_", v, "_SIZE=", len, sep = "")
      s <- c(s, paste("PRIMER_", v, "_TM=", tm, sep = ""))
      s <- c(s, sprintf("PRIMER_MIN_GC=%f", gc[1]), sprintf("PRIMER_MAX_GC=%f", gc[2]))
      #s <- c(s, sprintf("PRIMER_OPT_GC_PERCENT=%f", gc[2]))
      s <- c(s, sprintf("PRIMER_PRODUCT_SIZE_RANGE=%s", paste(product_size, collapse = "-")))
      
      seq <- readDNAStringSet(f_name)
      seq <- as.character(seq, use.names = F)
      s <- c(s, sprintf("SEQUENCE_TEMPLATE=%s", seq))
      s <- c(s, sprintf("PRIMER_NUM_RETURN=%d", n))
      
      template <- list.files(path = file.path(system.file(package = "PLAYRDesign")), pattern = "primer3_settings_template.txt", full.names = T)
      template <- readLines(template)
      cat(s, template, file = "PLAYRDesign_primer3_input.txt", sep = "\n")      
      system("primer3_core -output=PLAYRDesign_primer3_output.txt PLAYRDesign_primer3_input.txt")
      return(parse_primer3_output("PLAYRDesign_primer3_output.txt"))
}


refseq_to_symbol <- function(a)
{
      a <- gsub("\\.[0-9]*$", "", a)
      tt.a <- AnnotationFuncs::translate(a, from = org.Hs.egREFSEQ2EG, to = org.Hs.egSYMBOL)
      excl <- a[!(a %in% names(tt.a))]
      if(length(excl) > 0)
      {
            print(sprintf("Can't convert %s", paste(excl, collapse = ",")))
            temp <- as.list(excl)
            names(excl) <- excl
            tt.a <- c(tt.a, excl)
      }
      return(unlist(tt.a[a]))
      
}

parse_blast_result_txt <- function(f_name, filter_same_gi)
{
      #This assumes a single query sequence
      tab <- read.table(f_name, header = F, sep = "\t", quote = "", stringsAsFactors = F)
      names(tab) <- c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount", 
                      "gapOpenCount", "queryStart", "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore")
      tab <- tab[(tab$queryId != tab$subjectId),]
      if(filter_same_gi)
      {
            a <- sapply(strsplit(tab$queryId, "\\|"), function(x) {x[[4]]})
            b <- sapply(strsplit(tab$subjectId, "\\|"), function(x) {x[[4]]})
            a <- refseq_to_symbol(a)
            b <- refseq_to_symbol(b)
            tab <- tab[a != b,]
            
            
      }
      strand <- rep("+", nrow(tab))
      
      print(tab)
      #strand[tab$subjectEnd < tab$subjectStart] <- "-"
      res <- GRanges(seqnames = Rle(tab$queryId), strand = Rle(strand), 
                     ranges = IRanges(start = tab$queryStart, end = tab$queryEnd), percIdentity = tab$percIdentity, eVal = tab$eVal, bitScore = tab$bitScore,
                     gapOpenCount = tab$gapOpenCount, subjectId = tab$subjectId)
      print(res)
      return(res)
}

get_seq_from_ranges <- function(gr)
{
      gen <- BSgenome.Hsapiens.UCSC.hg19
      
      s <- getSeq(gen, names = gr)
      
      return(unlist(s))
}


get_exons_for_transcript <- function(id)
{
      txdb <- loadDb(system.file("UCSC_Refseq_transcripts.sqlite", package = "PLAYRDesign"))
      #Remove weird chromosomes
      tab <- select(txdb, keys = id, keytype="TXNAME", columns = columns(txdb))
      tab <- tab[grep("_", tab$EXONCHROM, invert = T),]
      u <- unique(tab$EXONCHROM)
      if(length(u) > 1)
      {
            print("Transcript maps to multiple chromosomes, only using first one")
            tab <- tab[tab$EXONCHROM == u[1],]
      }
      gr <- GRanges(seqnames = Rle(tab$EXONCHROM), ranges = IRanges(start = tab$EXONSTART, end = tab$EXONEND), 
                    strand = Rle(tab$EXONSTRAND))
      print(gr)
      return(gr)
}
      

aligned_sequence_to_index <- function(s)
{
      temp <- start(s) - 1
      s <- as.character(s)
      s <- strsplit(s, "")[[1]]
      
      ret <- sapply(s, function(x)
      {
            if(x != "-")
            {
                  temp <<- temp + 1
                  return(temp)
            }
            else
                  return(NA)
      })
      return(ret)
}


project_data_on_seq <- function(query, subject, v.subject)
{
      aln <- pairwiseAlignment(query, subject, type = "global")
      query.idx <- aligned_sequence_to_index(pattern(aln))
      subject.idx <- aligned_sequence_to_index(subject(aln))
      
      subject.idx <- subject.idx[!is.na(query.idx)]
      query.idx <- query.idx[!is.na(query.idx)]
      
      ret <- numeric(length(query))
      ret[query.idx] <- v.subject[subject.idx]
      
      return(ret)
}


get_refseq_id_from_fasta <- function(f_name)
{
      s <- names(readDNAStringSet(f_name))
      s <- strsplit(s, "\\|")[[1]][4]
      s <- strsplit(s, "\\.")[[1]][1]
      return(s)
}


get_seq_from_file <- function(f_name)
{
      return(unlist(readDNAStringSet(f_name)))
}

      

expand_est_alignments <- function(gr.est)
{
      l <- split(gr.est, names(gr.est))
      ret <- lapply(l, function(x)
      {
            start <- strsplit(mcols(x)$tStarts, ",")[[1]]
            size <- strsplit(mcols(x)$blockSizes, ",")[[1]]
            start <- as.numeric(start)
            size <- as.numeric(size)
            return(GRanges(seqnames = seqnames(x), ranges = IRanges(start = start, width = size), strand = strand(x)))
      })
      ret <- GRangesList(ret)
      print(ret)
      return(ret[names(gr.est)])
}


get_exons_skips <- function(gr.gene, gr.est)
{
      sel.est <- subsetByOverlaps(gr.est, gr.gene)
      sel.est.exp <- expand_est_alignments(sel.est)
      ret <-  numeric(length(gr.gene))
      names(ret) <- as.character(1:length(gr.gene))
      
      k <- table(queryHits(findOverlaps(gr.gene, sel.est)))
      ret[names(k)] <- k
      
      k <- table(queryHits(findOverlaps(gr.gene, sel.est.exp)))
      ret[names(k)] <- ret[names(k)] - k
      print("Exons sizes")
      print(width(gr.gene))
      print("Exons skips")
      print(ret)
      v <- width(gr.gene)
      v <- as.vector(Rle(ret, v))
      return(v)      
}


run_blast_analysis_for_seq <- function(f_name, db, filter_same_gi)
{
      seq_length <- nchar(readDNAStringSet(f_name))
      blast.f_name <- sprintf("%s_%s.blast_out.txt", f_name, db)
      print("Running BLAST")
      system(paste("blastn -db", db, "-query", f_name, "-task blastn -outfmt 6 >",  blast.f_name))
      print("Done")
      res <- parse_blast_result_txt(blast.f_name, filter_same_gi)
      dj <- disjoin_reduce_blast_ranges(res)
      df <- expand_ranges_to_dataframe(dj, seq_length)
      return(df)      
}

disjoin_reduce_blast_ranges <- function(gr)
{
      dj <- disjoin(gr)
      ov <- as.data.frame(findOverlaps(gr, dj))
      mc <- as.data.frame(mcols(gr))
      mc <- mc[, sapply(mc, is.numeric)]
      
      res <- ddply(ov, ~subjectHits, function(x, mc) {return(colMeans(mc[x$queryHits,]))}, mc = mc)
      mcols(dj) <- res[, names(res) != "subjectHits"]
      return(dj)
}
      
expand_ranges_to_dataframe <- function(gr, seq_length)
{
      mc <- as.data.frame(mcols(gr))
      max.end <- max(end(ranges(gr)))
      df <- data.frame(pos = 1:seq_length)
      
      for(i in 1:ncol(mc))
      {
            v <- Rle(0, seq_length)
            for(j in 1:length(ranges(gr)))
            {
                  r <- ranges(gr)[j]
                  v[r] <- mc[j, i]
            }
            df <- cbind(df, as.vector(v))
      }
      names(df) <- c("pos", names(mc))
      return(df)
}
      
      
                  
                  
                  
data_frame_factor_to_char <- function(df)
{
      i <- sapply(df, is.factor)
      df[i] <- lapply(df[i], as.character)
      return(df)
}
      
      
parse_blast_result_xml <- function(f_name)
{
      result <- xmlTreeParse(f_name, useInternalNodes=TRUE,
                              error = xmlErrorCumulator(immediate=FALSE))
      parse_single_hit <- function(x)
      {
            id <- unlist(xpathApply(x,  "./Hit_id", xmlValue))
            nums <- unlist(xpathApply(x,  "./Hit_hsps/Hsp/Hsp_num", xmlValue))
            qseq <- unlist(xpathApply(x,  "./Hit_hsps/Hsp/Hsp_qseq", xmlValue))
            hseq <- unlist(xpathApply(x,  "./Hit_hsps/Hsp/Hsp_hseq", xmlValue))
            return(data.frame(id, hsp_num = nums, qseq, hseq))
      }
      df <- xpathApply(result, "//Hit", parse_single_hit)
      df <- data.frame(do.call(rbind, df))
      df <- data_frame_factor_to_char(df)     
      res <- list()
      print(df)
      for (i in 1:nrow(df))
      {
            aln <- DNAMultipleAlignment(c(df$qseq[i], df$hseq[i]), 
                                           rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), "NormalIRanges"))
            res[i] <- aln
      }
      names(res) <- df$id
      return(res)
}



      
      
get_blast_result <-  function (x, database = "nr", hitListSize = "10", filter = "L", expect = "10", program = "blastn", equery = "",
                               attempts = 10) 
{
      baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
      query <- paste("QUERY=", as.character(x), "&DATABASE=", database, "&ENTREZ_QUERY=", equery,
                     "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
                     expect, "&PROGRAM=", program, sep = "")
      url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
      results <- tempfile()
      Sys.sleep(5)
      require(XML)
      post <- htmlTreeParse(url0, useInternalNodes = TRUE)
      x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
      rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
      rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
                             x))
      url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
                      rid)
      Sys.sleep(rtoe)
      .tryParseResult <- function(url, attempts)
            {
            for (i in 1:(attempts+1)) 
            {
                  result <- tryCatch({xmlTreeParse(url, useInternalNodes=TRUE,
                                     error = xmlErrorCumulator(immediate=FALSE))}, error=function(err) NULL)
                  if (!is.null(result)) return(result)
                  Sys.sleep(10)
            }
            stop(paste("no results after ", attempts, 
                       " attempts; please try again later", sep = ""))
      }
      result <- .tryParseResult(url1, attempts)
      return(result)
}

get_reduced_constitutive_exons <- function(gene, hexevent)
{
      tab <- hexevent[hexevent$genename == gene,]
      gr <- GRanges(seqnames = Rle(tab$chromo), ranges = IRanges(start = tab$start, end = tab$end), 
                    strand = Rle(tab$strand), constitLevel = tab$constitLevel, count = tab$count)
      red <- reduce(gr, with.mapping  = T)
      mapping <- mcols(red)$mapping
      #Do a weighted mean of constitLevel, where the weight is going to be given by count
      avg <- sapply(mapping, function(x, gr)
            {
                  m <- mcols(gr[x])
                  val <- sum(m$count * m$constitLevel) / sum(m$count) 
                  return(val)
            }, gr = gr)

      mcols(red)$avg <- avg
      return(red)
}
    
#gr <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)), ranges = IRanges(1:10, end = 7:16, names = head(letters, 10)), strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#            score = 1:10, GC = seq(1, 0, length = 10))