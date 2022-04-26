
#' compute_similarity_epigenome
#'
#' This function computes the epigenome profile similarities between query regions and
#' the corresponding list of target regions.
#'
#' @param gr_query
#' The query genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' @param gr_target_list
#' The corresponding candidate target regions of \code{\link[GenomicRanges]{GRangesList}} class.
#' The length of gr_target_list should be the same as the length of gr_query
#' @param query_grlist
#' The epigenomic datasets in query genome, \code{\link[GenomicRanges]{GRangesList}} class.
#' @param target_grlist
#' The epigenomic datasets in target genome, \code{\link[GenomicRanges]{GRangesList}} class.
#' Should match with query_grlist.
#' @param weights
#' The weights of epigenomic datasets. A numeric vector with length equal to the number of epigenomic datasets.
#' Default is NULL (equal weights).
#' @param maxgap
#' The maxgap allowed for overlapping. Default is 500L.
#' @param metric
#' The metric for computing similarities between indicators.
#' One of 'cosine' and 'jaccard'. Default is 'cosine'.
#' @param verbose
#' Print messages or not. Default is TRUE.
#'
#' @return
#' A \code{\link[GenomicRanges]{GRangesList}} class the same structure as gr_target_list
#' with an extra meta column named "epigenome" measuring the similarity score of epigenome
#' datasets.
#'
#' @examples
#'
#' \dontrun{
#' data('data_example')
#' ## the epigenomic datasets in this example can be downloaded from github
#' gr_target_list <- compute_similarity_epigenome(gr, gr_list, epigenome_mm10, epigenome_hg38)
#' }
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import Matrix
#' @rawNamespace import(data.table, except = shift)
#' @export
compute_similarity_epigenome <- function(
    gr_query,
    gr_target_list,
    query_grlist,
    target_grlist,
    weights = NULL,
    maxgap = 500L,
    metric = 'cosine',
    verbose = TRUE
) {

    stopifnot(
        class(gr_query) == 'GRanges',
        class(gr_target_list) %in% c('GRangesList', 'CompressedGRangesList'),
        length(gr_query) == length(gr_target_list),
        class(query_grlist) %in% c('GRangesList', 'CompressedGRangesList'),
        class(target_grlist) %in% c('GRangesList', 'CompressedGRangesList'),
        length(query_grlist) == length(target_grlist),
        is.null(weights) || length(weights) == length(query_grlist),
        length(query_grlist) > 1,
        metric %in% c('cosine', 'jaccard')
    )

    if(verbose) {
        message('Computing epigenome similarity scores.')
    }

    similarity <- compute_similarity_epigenome_flat(
        gr_query = gr_query[rep(1:length(gr_query), lengths(gr_target_list))],
        gr_target = unlist(gr_target_list),
        query_grlist = query_grlist,
        target_grlist = target_grlist,
        weights = weights,
        maxgap = maxgap,
        metric = metric
    )

    gr_target_df <- as.data.table(gr_target_list)
    gr_target_df$group <- factor(gr_target_df$group, levels = 1:length(gr_target_list))
    gr_target_df$`group_name` <- NULL
    gr_target_df$epigenome <- similarity

    gr_target_list <- makeGRangesListFromDataFrame(
        gr_target_df,
        split.field = 'group',
        keep.extra.columns = TRUE
    )
    names(gr_target_list) <- NULL

    return(gr_target_list)
}





#' compute_similarity_epigenome_flat
#'
#' This function computes the epigenome profile similarities between query regions and
#' the corresponding target regions.
#'
#'
#' @param gr_query
#' The query genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' @param gr_target
#' The target genomic regions of \code{\link[GenomicRanges]{GRanges}} class.
#' The length of gr_target has to be equal to the length of gr_query.
#' Or either gr_query or gr_target has length 1 that can be broadcasted.
#' @param query_grlist
#' The epigenomic datasets in query genome, \code{\link[GenomicRanges]{GRangesList}} class.
#' @param target_grlist
#' The epigenomic datasets in target genome, \code{\link[GenomicRanges]{GRangesList}} class.
#' Should match with query_grlist.
#' @param weights
#' The weights of epigenomic datasets. A numeric vector with length equal to the number of epigenomic datasets.
#' Default is NULL (equal weights).
#' @param maxgap
#' The maxgap allowed for overlapping. Default is 500L.
#' @param metric
#' The metric for computing similarities between indicators.
#' One of 'cosine' and 'jaccard'. Default is 'cosine'.
#' @return
#' A numeric vector with values between 0 and 1.
#'
#' @examples
#'
#' \dontrun{
#' data('data_example')
#' ## the epigenomic datasets in this example can be downloaded from github
#' similarity <- compute_similarity_epigenome_flat(gr[2], gr_list[[2]], epigenome_mm10, epigenome_hg38)
#' }
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import Matrix
#' @export
compute_similarity_epigenome_flat <- function(
    gr_query,
    gr_target,
    query_grlist,
    target_grlist,
    weights = NULL,
    maxgap = 500L,
    metric = 'cosine'
) {

    stopifnot(
        # gr_query and gr_target have to be GRanges classes with the same length
        class(gr_query) == 'GRanges',
        class(gr_target) == 'GRanges',
        length(gr_query) == 1 || length(gr_target) == 1 || length(gr_query) == length(gr_target),
        class(query_grlist) %in% c('GRangesList', 'CompressedGRangesList'),
        class(target_grlist) %in% c('GRangesList', 'CompressedGRangesList'),
        length(query_grlist) == length(target_grlist),
        is.null(weights) || length(weights) == length(query_grlist),
        length(query_grlist) > 1,
        metric %in% c('cosine', 'jaccard')
    )

    if(length(gr_query) == 0 || length(gr_target) == 0) {
        return(NULL) # has to be NULL
    }

    query_mat <- compute_epigenome_matrix(gr_query, query_grlist, maxgap)
    target_mat <- compute_epigenome_matrix(gr_target, target_grlist, maxgap)
    similarity <- compute_similarity_from_matrix_epigenome(query_mat, target_mat, metric, weights)
    return(similarity)
}





## A helper function that computes the epigenome signal indicator matrix
## (length(gr) by length(grList))
compute_epigenome_matrix <- function(
    gr,
    grlist,
    maxgap = 500L
) {
    stopifnot(
        class(gr) == 'GRanges',
        length(gr) > 0,
        maxgap >= 0,
        class(grlist) %in% c('GRangesList', 'CompressedGRangesList')
    )

    maxgap <- as.integer(maxgap)

    overlap <- findOverlaps(gr, grlist, maxgap = maxgap)

    epigenome_ix <- matrix(FALSE, nrow = length(gr), ncol = length(grlist))
    epigenome_ix[as.matrix(overlap)] <- TRUE

    return(epigenome_ix)
}



## A helper function that computes cosine/jaccard similarities between indicator matrices rowwisely.
compute_similarity_from_matrix_epigenome <- function(
    query_mat,
    target_mat,
    metric = 'cosine',
    weights = NULL
) {

    stopifnot(
        nrow(query_mat) == nrow(target_mat) || nrow(query_mat) == 1 || nrow(target_mat) == 1,
        ncol(query_mat) == ncol(target_mat),
        metric %in% c('cosine', 'jaccard'),
        is.null(weights) || length(weights) == ncol(query_mat)
    )

    if(nrow(query_mat) != nrow(target_mat)) {
        if(nrow(query_mat) == 1) {
            query_mat <- query_mat[rep(1, nrow(target_mat)), ]
        }

        else {
            target_mat <- target_mat[rep(1, nrow(query_mat)), ]
        }
    }

    if(is.null(weights)) {
        q_vec <- Matrix::rowSums(query_mat)
        t_vec <- Matrix::rowSums(target_mat)
        qt_vec <- Matrix::rowSums(query_mat & target_mat)
    }
    else {
        q_vec <- as.vector(query_mat %*% weights)
        t_vec <- as.vector(target_mat %*% weights)
        qt_vec <- as.vector((query_mat & target_mat) %*% weights)
    }

    if(metric == 'cosine'){
        similarity <- qt_vec / sqrt(q_vec * t_vec)
    }

    if(metric == 'jaccard'){
        similarity <- qt_vec / (q_vec + t_vec - qt_vec)
    }

    return(similarity)
}




