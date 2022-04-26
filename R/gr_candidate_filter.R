
#'
#' gr_candidate_filter
#'
#' This function takes GRangesList with 'epigenome' and 'grammar' meta columns
#' and filter the results based on these two columns.
#'
#' @param gr_candidate
#' A \code{\link[GenomicRanges]{GRangesList}} object with 'epigenome' and 'grammar' meta columns.
#' @param best_k
#' Select the best k regions for each element of gr_candidate. Default is 1L.
#' @param b_intercept
#' The regression coefficient of intercept. Default is -2.922484.
#' @param b_epigenome
#' The regression coefficient of epigenome similarity. Default is 5.032385.
#' @param b_grammar
#' The regression coefficient of grammar similarity. Default is 3.223607.
#' @param b_interaction
#' The regression coefficient of the interaction term. Default is NULL.
#' @param threshold
#' Discard the genomic regions with score below this threshold. Default is 0.
#' @param random
#' Pick best_k regions randomly for each element of gr_candidate. Default is FALSE.
#' @param verbose
#' Print messages or not. Default is TRUE.
#'
#' @return
#' A \code{\link[GenomicRanges]{GRangesList}} class the same structure as gr_candidate
#' with an extra meta column named "score" measuring the overall similarity score.
#' Any genomic regions with score below the threshold will be filtered away.
#'
#'
#' @examples
#'
#' data('data_example')
#' gr_list_filtered <- gr_candidate_filter(gr_list_score, threshold = 0.5)
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @rawNamespace import(data.table, except = shift)
#' @export
gr_candidate_filter <- function(
    gr_candidate,
    best_k = 1L,
    b_intercept = -2.922484,
    b_epigenome = 5.032385,
    b_grammar = 3.223607,
    b_interaction = NULL,
    threshold = 0,
    random = FALSE,
    verbose = TRUE
) {
    stopifnot(
        class(gr_candidate) %in% c('CompressedGRangesList', 'GrangesList'),
        best_k >= 1
    )
    best_k <- as.integer(best_k)

    gr_candidate_dt <- as.data.table(gr_candidate)
    gr_candidate_dt$`group_name` <- NULL

    stopifnot(
        'epigenome' %in% colnames(gr_candidate_dt),
        'grammar' %in% colnames(gr_candidate_dt)
    )

    if(verbose) {
        message('Filtering syntenic regions based on epigenome and grammar similarities.')
    }

    gr_candidate_dt[is.na(gr_candidate_dt)] <- 0
    score <- b_intercept + b_epigenome * gr_candidate_dt$epigenome + b_grammar * gr_candidate_dt$grammar

    if(!is.null(b_interaction)) {
        score <- score + b_interaction * gr_candidate_dt$epigenome * gr_candidate_dt$grammar
    }

    gr_candidate_dt[, score := 1/(1+exp(-score))]
    gr_candidate_dt <- gr_candidate_dt[score >= threshold]

    if(random) {
        gr_candidate_dt <- gr_candidate_dt[, .SD[sample(1:nrow(.SD), min(nrow(.SD), best_k))], by = group]
    }
    else {
        gr_candidate_dt <- gr_candidate_dt[, .SD[order(-score)][1:min(nrow(.SD), best_k)], by = group]
    }

    gr_candidate_dt$group <- factor(gr_candidate_dt$group, levels = 1:length(gr_candidate))

    gr_candidate <- makeGRangesListFromDataFrame(
        gr_candidate_dt,
        split.field = 'group',
        keep.extra.columns = TRUE
    )

    names(gr_candidate) <- NULL

    return(gr_candidate)
}
