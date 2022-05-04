

#'
#' training_module
#'
#' This function takes GRangesList with 'epigenome' and 'grammar' meta columns
#' and evaluate the predictive power of overlapping with the true orthologous genomic regions.
#'
#' @param gr_candidate
#' A \code{\link[GenomicRanges]{GRangesList}} object with 'epigenome' and 'grammar' meta columns.
#' @param gr_true
#' A \code{\link[GenomicRanges]{GRanges}} object as the true orthologous genomic regions
#' that we aim to derive.
#' @param max_filter_proportion
#' The maximum proportion of candidate target regions that are allowed to be filtered out. Default is 0.4.
#' @param interaction
#' If we include the interaction term in the logistic regression. Default is FALSE.
#'
#' @return
#' A data.table with the following columns:
#' \tabular{ll}{
#' \code{auroc} \tab The optimal AUROC. \cr
#' \code{aupr} \tab The optimal AUPR. \cr
#' \code{b.intercept} \tab The estimated regression coefficient for the intercept. \cr
#' \code{b.epigenome} \tab The estimated regression coefficient for the epigenome signal. \cr
#' \code{b.grammar} \tab The estimated regression coefficient for the sequence grammar. \cr
#' \code{b.interaction} \tab The estimated regression coefficient for the interaction term if included. \cr
#' \code{max_filter_proportion} \tab The input parameter max_filter_proportion. \cr
#' \code{threshold} \tab The suggested score threshold. \cr
#' \code{precision} \tab The precision associated with the suggested score threshold. \cr
#' }
#'
#'
#' @examples
#'
#' data('training_module_example')
#' res <- training_module(gr_candidate, gr_true)
#'
#'
#' @author Chenyang Dong \email{cdong@stat.wisc.edu}
#' @import GenomicRanges
#' @import PRROC
#' @rawNamespace import(data.table, except = shift)
#' @export
training_module <- function(gr_candidate,
                            gr_true,
                            max_filter_proportion = 0.4,
                            interaction = FALSE) {
    stopifnot(
        class(gr_candidate) %in% c('GRangesList', 'CompressedGRangesList'),
        class(gr_true) == 'GRanges',
        max_filter_proportion >= 0 && max_filter_proportion <= 1
    )

    gr_dt <- as.data.table(gr_candidate)

    stopifnot(
        nrow(gr_dt) >= 100,
        # there are too few candidate target regions
        'epigenome' %in% colnames(gr_dt),
        'grammar' %in% colnames(gr_dt)
    )

    gr_dt$`group_name` <- NULL
    gr_dt$group <-
        factor(gr_dt$group, levels = 1:length(gr_candidate))

    gr <- makeGRangesFromDataFrame(gr_dt)
    gr_dt$label <- gr %over% gr_true

    gr_dt[is.na(gr_dt)] <- 0

    if (interaction) {
        model <-
            glm(
                label ~ epigenome + grammar + epigenome * grammar,
                family = 'binomial',
                data = gr_dt
            )
        coeff <- as.vector(summary(model)$coefficients[, 1])

        gr_dt$logit <- predict(model, type = 'response')
        roc <-
            roc.curve(scores.class0 = gr_dt[label == TRUE]$logit,
                      scores.class1 = gr_dt[label == FALSE]$logit)
        pr <-
            pr.curve(scores.class0 = gr_dt[label == TRUE]$logit,
                     scores.class1 = gr_dt[label == FALSE]$logit)

        gr_dt_max <-
            gr_dt[, .(logit = max(logit), label = any(label)), by = group]
        gr_dt_max <- gr_dt_max[order(-logit)]

        gr_dt_max$precision <-
            cumsum(gr_dt_max$label) / (1:nrow(gr_dt_max))
        gr_dt_max <- gr_dt_max[order(logit)]

        max_filter_number <-
            max(1, as.integer(nrow(gr_dt_max) * max_filter_proportion))
        ind <- which.max(gr_dt_max$precision[1:max_filter_number])

        threshold <- gr_dt_max$logit[ind]
        precision <- gr_dt_max$precision[ind]

        return(
            data.table(
                auroc = roc$auc,
                aupr = pr$auc.integral,
                b.intercept = coeff[1],
                b.epigenome = coeff[2],
                b.grammar = coeff[3],
                b.interaction = coeff[4],
                max_filter_proportion = max_filter_proportion,
                threshold = threshold,
                precision = precision
            )
        )
    }
    else {
        model <-
            glm(label ~ epigenome + grammar,
                family = 'binomial',
                data = gr_dt)
        coeff <- as.vector(summary(model)$coefficients[, 1])

        gr_dt$logit <- predict(model, type = 'response')
        roc <-
            roc.curve(scores.class0 = gr_dt[label == TRUE]$logit,
                      scores.class1 = gr_dt[label == FALSE]$logit)
        pr <-
            pr.curve(scores.class0 = gr_dt[label == TRUE]$logit,
                     scores.class1 = gr_dt[label == FALSE]$logit)

        gr_dt_max <-
            gr_dt[, .(logit = max(logit), label = any(label)), by = group]
        gr_dt_max <- gr_dt_max[order(-logit)]

        gr_dt_max$precision <-
            cumsum(gr_dt_max$label) / (1:nrow(gr_dt_max))
        gr_dt_max <- gr_dt_max[order(logit)]

        max_filter_number <-
            max(1, as.integer(nrow(gr_dt_max) * max_filter_proportion))
        ind <- which.max(gr_dt_max$precision[1:max_filter_number])

        threshold <- gr_dt_max$logit[ind]
        precision <- gr_dt_max$precision[ind]

        return(
            data.table(
                auroc = roc$auc,
                aupr = pr$auc.integral,
                b.intercept = coeff[1],
                b.epigenome = coeff[2],
                b.grammar = coeff[3],
                max_filter_proportion = max_filter_proportion,
                threshold = threshold,
                precision = precision
            )
        )
    }
}
