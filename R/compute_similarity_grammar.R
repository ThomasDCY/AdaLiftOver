
#' compute_similarity_grammar
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding list of target regions.
#'
#' @param gr_query
#' The query genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' @param gr_target_list
#' The corresponding candidate target regions of \code{\link[GenomicRanges]{GRangesList}} class.
#' The length of gr_target_list should be the same as the length of gr_query
#' @param query_genome
#' The query genome build identifier (e.g. 'mm10', 'hg38').
#' Make sure the corresponding BSgenome package is installed.
#' Or a \code{\link[BSgenome]{BSgenome}} class.
#' @param target_genome
#' The target genome build identifier (e.g. 'mm10', 'hg38').
#' Make sure the corresponding BSgenome package is installed.
#' Or a \code{\link[BSgenome]{BSgenome}} class.
#' @param motif_list
#' \code{\link[TFBSTools]{PWMatrixList}}, or \code{\link[TFBSTools]{PFMatrixList}}.
#' @param grammar_size
#' The nearby region (bp) to extract regulatory grammars. Default is 500L.
#' @param metric
#' The metric for computing similarities between indicators.
#' One of 'cosine' and 'jaccard'. Default is 'cosine'.
#' @param verbose
#' Print messages or not. Default is TRUE.
#'
#' @return
#' A \code{\link[GenomicRanges]{GRangesList}} class the same structure as gr_target_list
#' with an extra meta column named "grammar" measuring the similarity score of sequence grammars.
#'
#' @examples
#'
#'
#' data('data_example')
#' data('jaspar_pfm_list')
#' gr_target_list <- compute_similarity_grammar(gr, gr_list, 'mm10', 'hg38', jaspar_pfm_list)
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import motifmatchr
#' @import Matrix
#' @rawNamespace import(data.table, except = shift)
#' @export
compute_similarity_grammar <- function(gr_query,
                                       gr_target_list,
                                       query_genome,
                                       target_genome,
                                       motif_list,
                                       grammar_size = 500L,
                                       metric = 'cosine',
                                       verbose = TRUE) {
  stopifnot(
    class(gr_query) == 'GRanges',
    class(gr_target_list) %in% c('GRangesList', 'CompressedGRangesList'),
    length(gr_query) == length(gr_target_list),
    metric %in% c('cosine', 'jaccard')
  )

  if (verbose) {
    message('Computing grammar similarity scores.')
  }

  similarity <- compute_similarity_grammar_flat(
    gr_query = gr_query[rep(1:length(gr_query), lengths(gr_target_list))],
    gr_target = unlist(gr_target_list),
    query_genome = query_genome,
    target_genome = target_genome,
    motif_list = motif_list,
    grammar_size = grammar_size,
    metric = metric
  )

  gr_target_df <- as.data.table(gr_target_list)
  gr_target_df$group <-
    factor(gr_target_df$group, levels = 1:length(gr_target_list))
  gr_target_df$`group_name` <- NULL
  gr_target_df$grammar <- similarity

  gr_target_list <- makeGRangesListFromDataFrame(gr_target_df,
                                                 split.field = 'group',
                                                 keep.extra.columns = TRUE)
  names(gr_target_list) <- NULL

  return(gr_target_list)
}




#' compute_similarity_grammar_flat
#'
#' This function computes the sequence grammar similarities between query regions and
#' the corresponding target regions.
#'
#' @param gr_query
#' The query genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' @param gr_target
#' The target genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' The length of gr_target should be the same as the length of gr_query
#' @param query_genome
#' The query genome build identifier (e.g. 'mm10', 'hg38').
#' Make sure the corresponding BSgenome package is installed.
#' Or a \code{\link[BSgenome]{BSgenome}} class.
#' @param target_genome
#' The target genome build identifier (e.g. 'mm10', 'hg38').
#' Make sure the corresponding BSgenome package is installed.
#' Or a \code{\link[BSgenome]{BSgenome}} class.
#' @param motif_list
#' \code{\link[TFBSTools]{PWMatrixList}}, or \code{\link[TFBSTools]{PFMatrixList}}.
#' @param grammar_size
#' The nearby region (bp) to extract regulatory grammars. Default is 500L.
#' @param metric
#' The metric for computing similarities between indicators.
#' One of 'cosine' and 'jaccard'. Default is 'cosine'.
#'
#' @return
#' A numeric vector with values between 0 and 1.
#'
#' @examples
#'
#' data('data_example')
#' data('jaspar_pfm_list')
#' similarity <- compute_similarity_grammar_flat(gr[2], gr_list[[2]], 'mm10', 'hg38', jaspar_pfm_list)
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import motifmatchr
#' @import Matrix
#' @rawNamespace import(data.table, except = shift)
#' @export
compute_similarity_grammar_flat <- function(gr_query,
                                            gr_target,
                                            query_genome,
                                            target_genome,
                                            motif_list,
                                            grammar_size = 500L,
                                            metric = 'cosine') {
  stopifnot(
    # gr_query and gr_target have to be GRanges classes with the same length
    class(gr_query) == 'GRanges',
    class(gr_target) == 'GRanges',
    length(gr_query) == 1 ||
      length(gr_target) == 1 ||
      length(gr_query) == length(gr_target),
    metric %in% c('cosine', 'jaccard'),
    grammar_size >= 100
  )

  if (length(gr_query) == 0 || length(gr_target) == 0) {
    return(NULL) # has to be NULL
  }

  query_mat <-
    compute_grammar_matrix(gr_query, motif_list, query_genome, grammar_size)
  target_mat <-
    compute_grammar_matrix(gr_target, motif_list, target_genome, grammar_size)
  similarity <-
    compute_similarity_from_matrix_grammar(query_mat, target_mat, metric)
  return(similarity)
}




## A helper function that computes the sequence grammar signal indicator matrix
## (length(gr) by length(motif_list))
compute_grammar_matrix <- function(gr,
                                       motif_list,
                                       genome,
                                       grammar_size = 500L) {
  stopifnot(class(gr) == 'GRanges',
            length(gr) > 0,
            grammar_size >= 150)
  grammar_size <- as.integer(grammar_size)
  if (is.character(genome)) {
    suppressPackageStartupMessages(library(BSgenome))
    genome_obj <- getBSgenome(genome)
  } else {
    genome_obj <- genome
  }
  chr_sizes <- seqlengths(genome_obj)

  gr_safe <- gr

  for (i in seq_along(gr_safe)) {
    chr_name <- as.character(seqnames(gr_safe[i]))
    chr_len <- chr_sizes[chr_name]

    if (is.na(chr_len)) {
      warning("Chromosome ", chr_name, " not found in genome. Skipping.")
      next
    }

    region_width <- width(gr_safe[i])
    if (region_width < grammar_size) {
      half_extra <- (grammar_size - region_width) %/% 2
      new_start <- max(1, start(gr_safe[i]) - half_extra)
      new_end <- min(chr_len, end(gr_safe[i]) + half_extra)
      start(gr_safe[i]) <- new_start
      end(gr_safe[i]) <- new_end
    }

    if (start(gr_safe[i]) < 1) start(gr_safe[i]) <- 1
    if (end(gr_safe[i]) > chr_len) end(gr_safe[i]) <- chr_len
  }
  tryCatch({
    grammar_ix <- matchMotifs(motif_list, gr_safe, genome = genome, bg = 'even')
    grammar_ix <- motifMatches(grammar_ix)
    return(grammar_ix)
  }, error = function(e) {
    for (i in seq_along(gr_safe)) {
      chr_name <- as.character(seqnames(gr_safe[i]))
      chr_len <- chr_sizes[chr_name]
      cat("Region ", i, ": ", chr_name, ":", start(gr_safe[i]), "-", end(gr_safe[i]),
          " (chromosome length: ", chr_len, ")\n", sep="")
    }
    stop("Error in matchMotifs: ", e$message)
  })
}






## A helper function that computes cosine/jaccard similarities between indicator matrices rowwisely.
compute_similarity_from_matrix_grammar <- function(query_mat,
                                                   target_mat,
                                                   metric = 'cosine') {
  stopifnot(
    nrow(query_mat) == nrow(target_mat) ||
      nrow(query_mat) == 1 || nrow(target_mat) == 1,
    ncol(query_mat) == ncol(target_mat),
    metric %in% c('cosine', 'jaccard')
  )

  if (nrow(query_mat) != nrow(target_mat)) {
    if (nrow(query_mat) == 1) {
      query_mat <- query_mat[rep(1, nrow(target_mat)),]
    }

    else {
      target_mat <- target_mat[rep(1, nrow(query_mat)),]
    }
  }

  q_vec <- Matrix::rowSums(query_mat)
  t_vec <- Matrix::rowSums(target_mat)
  qt_vec <- Matrix::rowSums(query_mat & target_mat)

  if (metric == 'cosine') {
    similarity <- qt_vec / sqrt(q_vec * t_vec)
  }

  if (metric == 'jaccard') {
    similarity <- qt_vec / (q_vec + t_vec - qt_vec)
  }

  return(similarity)
}
