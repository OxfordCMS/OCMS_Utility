#' alpha_diversity
#'
#' Calculates alpha diversity metrics from count table: Shannon's D,
#' Simpson's index, evenness, richness. Calculations are based on
#' Jost et al 2006 https://doi.org/10.1111/j.2006.0030-1299.14714.x
#'
#' @param count_df Count table with samples in columns, taxa in rows
#'
#' @returns data frame of alpha diversity metrics for each sample
#'
#' @import vegan
#' @import dplyr
#' @import tibble
#' @export
#'
#' @examples
#' # get example data
#' data(asv_example)
#'
#' # rownames have to be features
#' asv_counts <- data.frame(asv_example[2:ncol(asv_example)], row.names=asv_example$sequence)
#' alpha_div <- alpha_diversity(asv_counts)
#'

alpha_diversity <- function(count_df) {
  alpha <- vegan::diversity(count_df, index='shannon', MARGIN=2)
  alpha <- as.data.frame(alpha)
  alpha <- tibble::rownames_to_column(alpha, 'sample_id')
  colnames(alpha)[2] <- 'shannonH'
  alpha$shannonD <- exp(alpha$shannonH)

  # richness
  richness <- as.data.frame(specnumber(t(count_df))) %>%
    tibble::rownames_to_column('sample_id') %>%
    rename(richness="specnumber(t(count_df))")

  alpha <- alpha %>%
    left_join(richness, 'sample_id') %>%
    mutate(evenness = shannonH / log(richness))

  return(alpha)
}