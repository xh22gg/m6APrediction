#' Encode DNA sequence strings into a positional factor table
#'
#' Splits each DNA sequence (for example a 5-mer) into single bases and returns
#' a data frame with one factor column per position, using the fixed levels
#' A/T/C/G. This is the encoding expected by the fitted random forest model.
#'
#' @param dna_strings A character vector of DNA sequences of equal length.
#' @return A data.frame with one factor column per position, named nt_pos1 up
#'   to nt_posN.
#' @examples
#' dna_encoding(c("GGACA", "AGACT"))
dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A sites for many genomic positions
#'
#' Uses a fitted random forest model to predict the m6A modification probability
#' for each row of a feature table, then labels each site as Positive or
#' Negative using a probability threshold.
#'
#' @param ml_fit A fitted randomForest model (an RDS object).
#' @param feature_df A data.frame with the required feature columns:
#'   gc_content, RNA_type, RNA_region, exon_length, distance_to_junction,
#'   evolutionary_conservation and DNA_5mer.
#' @param positive_threshold Probability above which a site is labelled
#'   Positive (default 0.5). Values equal to the threshold are Negative.
#' @return The input data.frame with two added columns: predicted_m6A_prob
#'   (numeric) and predicted_m6A_status ("Positive"/"Negative").
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' example_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, head(example_df, 5))
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df)
}

#' Predict the m6A status of a single genomic position
#'
#' Convenience wrapper that predicts the m6A probability and status for one site
#' from its individual feature values.
#'
#' @param ml_fit A fitted randomForest model (an RDS object).
#' @param gc_content Numeric GC content (0 to 1).
#' @param RNA_type One of "mRNA", "lincRNA", "lncRNA", "pseudogene".
#' @param RNA_region One of "CDS", "intron", "3'UTR", "5'UTR".
#' @param exon_length Numeric, log2-transformed exon length.
#' @param distance_to_junction Numeric, log2-transformed distance to junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer A five-base context with A in the central position (for
#'   example "GGACA").
#' @param positive_threshold Probability threshold (default 0.5).
#' @return A named character vector with elements predicted_m6A_prob and
#'   predicted_m6A_status.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                               package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA", RNA_region = "CDS",
#'                   exon_length = 12, distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA")
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector)
}
