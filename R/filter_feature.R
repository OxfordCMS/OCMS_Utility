#' filter_feature
#'
#' Filter features from a count table based on abundance and prevalence criteria, 
#' with adjustable behavior for ASV or metagenomics data.
#'
#' @param count_df dataframe. Count table with samples in columns and features in rows. Feature ID in rownames.
#' @param tax_df dataframe. FeatureID must match rownames of \code{count_df}.
#'     Required if \code{data_type = "asv"}. Contains columns \code{'featureID','Kingdom','Phylum',
#'     'Class','Order','Family','Genus','Species'}.
#' @param filter_method character. One of \code{"abs_count", "percent_sample", "percent_dataset"}.
#'     Determines how the abundance cutoff is applied.
#' @param cutoff numeric. The abundance threshold:
#'     - If \code{filter_method = 'abs_count'}, it is the minimum read count.
#'     - If \code{filter_method = 'percent_sample' or 'percent_dataset'}, it is the minimum relative abundance (0–100).
#' @param prev_cutoff integer. The minimum number of samples in which a feature must meet the \code{cutoff} to be retained.
#' @param data_type character. One of \code{"asv","metagenomics"}.
#'     - If \code{"asv"}, the filtering includes taxonomic data and outputs a filtered taxonomy table.
#'     - If \code{"metagenomics"}, taxonomy data is not used.
#'
#' @return
#'     A list containing:
#'     \itemize{
#'       \item \code{taxonomy}: Filtered taxonomy table (only if \code{data_type = "asv"}).
#'       \item \code{filtered}: Filtered count table.
#'       \item \code{p_agg}: Plot of features kept/removed based on aggregated abundance and prevalence.
#'       \item \code{p_exp}: Expanded plot showing abundance in each sample.
#'       \item \code{feature_keep}: Vector of features retained after filtering.
#'       \item \code{feature_remove}: Vector of features removed during filtering.
#'       \item \code{msg}: Summary message describing filtering results.
#'     }
#'
#' @details
#' Filtering is based on feature abundance and sample prevalence. 
#' \code{filter_method} determines how \code{cutoff} is interpreted:
#' \itemize{
#'   \item \code{"abs_count"}: \code{cutoff} is an absolute read count.
#'   \item \code{"percent_sample"}: \code{cutoff} is a percentage of each sample's total reads.
#'   \item \code{"percent_dataset"}: \code{cutoff} is a percentage of the total reads in the entire dataset.
#' }
#'
#' \code{prev_cutoff} sets the minimum number of samples that must contain the feature at or above the \code{cutoff}.
#'
#' For \code{data_type = "asv"}, a taxonomy table is expected and included in the output.
#' For \code{data_type = "metagenomics"}, no taxonomy table is needed or returned.
#'
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @export
#'
#' @examples
#' # Example for ASV data:
#' # data(dss_example)
#' # tax_df <- dss_example$merged_taxonomy
#' # count_df <- dss_example$merged_abundance_id %>%
#' #   column_to_rownames('featureID')
#' # filtered_ls <- filter_feature(count_df, tax_df, 'percent_sample', 0.001, 2, data_type = 'asv')
#'
#' # Example for metagenomics data:
#' # Suppose counts_table is a count matrix of metagenomics features:
#' # filtered_results <- filter_feature(count_df = counts_table, filter_method = 'percent_sample', 
#' #                                    cutoff = 0.01, prev_cutoff = 2, data_type = 'metagenomics')
filter_feature <- function(count_df, tax_df = NULL,
                           filter_method = 'abs_count',
                           cutoff = 1, prev_cutoff = 2,
                           data_type = 'asv') {

  # input checks ---------------------------------------------------------------
  if(!data_type %in% c('asv','metagenomics')) {
    stop("data_type must be 'asv' or 'metagenomics'")
  }

  if(class(count_df) != 'data.frame') {
    stop("count_df must be a dataframe")
  }

  if(!filter_method %in% c('abs_count','percent_sample','percent_dataset')) {
    stop("filter_method must be 'abs_count','percent_sample', or 'percent_dataset'")
  }

  if(prev_cutoff > ncol(count_df)) {
    stop("prev_cutoff cannot exceed the number of samples")
  }

  if(filter_method == 'abs_count') {
    if(cutoff %% 1 != 0) {
      stop("cutoff must be an integer when filter_method = 'abs_count'")
    }
    if(cutoff > max(count_df)) {
      stop("cutoff cannot be greater than the maximum number of reads when filter_method = 'abs_count'")
    }
  } else {
    if(cutoff < 0 | cutoff > 100) {
      stop("cutoff must be between 0-100 when filter_method is 'percent_sample' or 'percent_dataset'")
    }
  }

  if(data_type == 'asv') {
    if(is.null(tax_df)) {
      stop("tax_df is required when data_type = 'asv'")
    }
    if(!identical(rownames(count_df), tax_df$featureID)) {
      stop("Order of features in count_df and tax_df must be identical")
    }
  }

  # calculate relative abundances -----------------------------------------------
  relab_by_sample <- apply(count_df, 2, function(x) x / sum(x))
  relab_by_data <- count_df / sum(count_df)
  prev_data <- rowSums(count_df != 0)

  # perform filtering ----------------------------------------------------------
  filtered_count <- switch(
    filter_method,
    abs_count = count_df[rowSums(count_df >= cutoff) >= prev_cutoff,],
    percent_sample = relab_by_sample[rowSums(relab_by_sample >= cutoff) >= prev_cutoff,],
    percent_dataset = relab_by_data[rowSums(relab_by_data >= cutoff) >= prev_cutoff,]
  )

  if(is.null(filtered_count)) {
    stop("No features remain after filtering. Please select different filter cutoff.")
  }

  # remove samples with no reads -----------------------------------------------
  sampleID <- colnames(filtered_count)
  empty_sample_ind <- which(colSums(filtered_count)==0)
  if(length(empty_sample_ind) > 0) {
    filtered_count <- filtered_count[, sampleID[-empty_sample_ind]]
    message(sprintf("%s samples contained 0 reads after filtering and have been removed: %s",
                    length(empty_sample_ind), paste(sampleID[empty_sample_ind], collapse=", ")))
  }

  # remove empty features ------------------------------------------------------
  featID <- rownames(filtered_count)
  empty_feat_ind <- which(rowSums(filtered_count)==0)
  if(length(empty_feat_ind) > 0) {
    filtered_count <- filtered_count[featID[-empty_feat_ind],]
    message(sprintf("%s features contained 0 reads after filtering and have been removed: %s",
                    length(empty_feat_ind), paste(featID[empty_feat_ind], collapse=", ")))
  }

  # record kept and removed features -------------------------------------------
  feature_keep <- rownames(filtered_count)
  feature_remove <- setdiff(rownames(count_df), rownames(filtered_count))

  # filter taxonomy if asv -----------------------------------------------------
  if(data_type == 'asv') {
    tax_filtered <- tax_df %>%
      filter(featureID %in% feature_keep)
  } else {
    tax_filtered <- NULL
  }

  # preparing plot data --------------------------------------------------------
  if(filter_method == 'abs_count') {
    pdata_val <- count_df
    ylab <- 'Read count'
  } else if(filter_method == 'percent_dataset') {
    pdata_val <- relab_by_data
    ylab <- 'Relative abundance\n(% of total reads)'
  } else {
    pdata_val <- relab_by_sample
    ylab <- 'Relative abundance\n(% of sample)'
  }

  pdata <- pdata_val %>%
    as.data.frame() %>%
    tibble::rownames_to_column('featureID') %>%
    pivot_longer(cols = -featureID, names_to = 'sampleID', values_to = 'value') %>%
    left_join(enframe(prev_data, 'featureID','prevalence'), 'featureID') %>%
    mutate(colour = ifelse(featureID %in% feature_keep, 'to keep', 'to remove'))

  pdata_agg <- pdata %>%
    group_by(featureID) %>%
    mutate(y = mean(value)) %>%
    distinct(prevalence, y, colour)

  scaleFUN <- function(x) sprintf("%.0f", x)

  p_agg <- ggplot(pdata_agg, aes(x = prevalence, y = y, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab(sprintf("Mean %s", ylab)) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, ncol(count_df))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  p_exp <- ggplot(pdata, aes(x = prevalence, y = value, colour = colour)) +
    geom_point(alpha = 0.5, size = 2) +
    xlab('Feature prevalence (# of samples)') +
    ylab(ylab) +
    scale_x_continuous(labels = scaleFUN, limits = c(0, length(sampleID))) +
    scale_y_continuous(trans = 'log10') +
    theme_bw(12) +
    theme(legend.position = 'bottom')

  # summary message ------------------------------------------------------------
  msg1 <- length(feature_keep)
  msg2 <- nrow(count_df)
  msg3 <- round(msg1 / msg2 * 100, 2)
  msg4 <- cutoff
  msg5 <- prev_cutoff
  msg6 <- ncol(count_df)
  msg7 <- round(msg5 / msg6 * 100, 1)

  if(filter_method == 'abs_count') {
    msg <- sprintf('Kept %s/%s (%s%%) features with read counts >= %s in >= %s/%s (%s%%) samples', 
                   msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_dataset') {
    msg <- sprintf('Kept %s/%s (%s%%) features with >= %s%% of dataset reads in >= %s/%s (%s%%) samples', 
                   msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  } else if(filter_method == 'percent_sample') {
    msg <- sprintf('Kept %s/%s (%s%%) features with >= %s%% of sample reads in >= %s/%s (%s%%) samples', 
                   msg1, msg2, msg3, msg4, msg5, msg6, msg7)
  }

  message(msg)

  # return ---------------------------------------------------------------------
  res <- list(
    'filtered' = filtered_count,
    'p_agg' = p_agg,
    'p_exp' = p_exp,
    'feature_remove' = feature_remove,
    'feature_keep' = feature_keep,
    'msg' = msg
  )
  
  if(data_type == 'asv') {
    res$taxonomy <- tax_filtered
  }

  return(res)
}
