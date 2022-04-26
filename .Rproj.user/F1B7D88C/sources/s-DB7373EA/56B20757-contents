
#' adaptive_liftover
#'
#' Generate candidate target regions
#'
#' @param gr
#' The query regions of \code{\link[GenomicRanges]{GRanges}} class.
#' @param chain
#' The chain object obtained from \code{rtracklayer::import.chain()}.
#' @param window
#' The window size (bp) around the query regions. Default is 2000L.
#' @param option
#' The option for adaptive_liftover
#' option = 'direct', apply direct liftOver to the query regions;
#' option = 'local', apply local liftOver by extending the query regions by a window size;
#' option = 'adaptive', apply direct liftOver first, for these failing to map, apply local liftOver.
#' @param gap
#' Ignore the gaps among target regions if below this value (bp). Default it 5L.
#' @param step_size
#' The step_size for generating candidate target regions. Default is 200L.
#' @param verbose
#' Print messages or not. Default is TRUE.
#'
#' @return
#' The candidate target regions of \code{\link[GenomicRanges]{GRangesList}} class with
#' length equal to the query regions.
#'
#' @examples
#'
#' \dontrun{
#' data('data_example')
#' ## the reciprocal chain file can be downloaded from UCSC
#' chain <- import.chain('mm10.hg38.rbest.chain')
#' gr_target_list <- adaptive_liftover(gr, chain, window = 2000, step_size = 10000)
#' gr_target_list <- adaptive_liftover(gr, chain, window = 2000, step_size = 200)
#' }
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import rtracklayer
#' @export
adaptive_liftover <- function(
    gr,
    chain,
    window = 2000L,
    option = 'adaptive',
    gap = 5L,
    step_size = 200L,
    verbose = TRUE
) {

    if(verbose) {
        message('Computing orthologous regions.')
    }

    gr_syntenic <- adaptive_liftover_syntenic(
        gr = gr,
        chain = chain,
        window = window,
        option = option
    )


    if(verbose) {
        message('Merging small gaps.')
    }

    stopifnot(
        gap >= 0
    )

    gap <- as.integer(gap)
    gr_syntenic <- reduce(gr_syntenic + gap)


    if(verbose) {
        message('Generating candidate regions.')
    }
    gr_syntenic <- generate_candidate_regions(
        gr_query = gr,
        gr_target_list = gr_syntenic,
        step_size = step_size
    )

    return(gr_syntenic)
}





## This helper function calls liftOver with three options
adaptive_liftover_syntenic <- function(
    gr,
    chain,
    window = 2000L,
    option = 'adaptive'
) {

    stopifnot(
        class(gr) == 'GRanges',
        option %in% c('adaptive', 'direct', 'local'),
        window >= 0
    )

    window <- as.integer(window)

    if(option == 'direct') {
        gr_syntenic <- liftOver(gr, chain)
    }

    else if(option == 'local') {
        gr_syntenic <- liftOver(gr+window, chain)
    }

    else {
        gr_syntenic <- liftOver(gr, chain)
        ind <- sapply(gr_syntenic, length) == 0
        suppressWarnings(gr_syntenic[ind] <- liftOver(gr[ind]+window, chain))
        # to suppress potential warnings from .Seqinfo.mergexy(x, y)
    }
    return(gr_syntenic)
}






## The main idea of this helper function is to generate the candidate regions
## by using the center points of all disjoint segments of target regions.
## Then move these center points towards both ends of the segments with step_size to generate anchor points.
## With these anchor points, we then resize the regions with widths equal to gr_query.
generate_candidate_regions <- function(
    gr_query,
    gr_target_list,
    step_size = 200L
) {

    stopifnot(
        class(gr_query) == 'GRanges',
        class(gr_target_list) %in% c('CompressedGRangesList', 'GrangesList'),
        length(gr_query) == length(gr_target_list),
        length(gr_query) > 0,
        step_size > 1
    )

    step_size <- as.integer(step_size)

    gr_target_df <- as.data.table(gr_target_list)

    if(nrow(gr_target_df) == 0) {
        warning('None of the input query regions can be mapped.')
        return(gr_target_list)
    }

    gr_target_df[, mid := (start + end) %/% 2]
    gr_target_df[, start := start + (mid - start) %% step_size]

    gr_candidate <- lapply(1:nrow(gr_target_df), function(i) {
        points <- seq(gr_target_df$start[i], gr_target_df$end[i], step_size)
        return(
            data.table(
                group = gr_target_df$group[i],
                seqnames = gr_target_df$seqnames[i],
                start = points,
                end = points,
                strand = '*'
            )
        )
    })
    gr_candidate <- rbindlist(gr_candidate)
    gr_candidate$group <- factor(gr_candidate$group, levels = 1:length(gr_target_list))

    ext <- (width(gr_query)[gr_candidate$group]) %/% 2
    gr_candidate$start <- gr_candidate$start - ext
    gr_candidate$end <- gr_candidate$end + ext

    gr_candidate <- makeGRangesListFromDataFrame(
        gr_candidate,
        split.field = 'group',
        keep.extra.columns = FALSE
    )
    names(gr_candidate) <- NULL

    return(gr_candidate)
}

