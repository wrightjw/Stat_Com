# score_se: Compute the Squared Error score
# Input:
#   pred : data.frame with (at least) a column "mu"
#   y : data vector
# Output:
#   a vector of scores

score_se <- function(pred, y) {
  (y - pred$mu)^2
}

# score_ds: Compute the Dawid-Sebastiani score
# Input:
#   pred : data.frame with (at least) columns "mu" and "sigma"
#   y : data vector
# Output:
#   a vector of scores

score_ds <- function(pred, y) {
  ((y - pred$mu) / pred$sigma)^2 + 2 * log(pred$sigma)
}

# mean_score: Compute average scores
# Input:
#   data: A data frame containing columns with computed scores
#         Scores in columns named SE, SD, and Brier are handled.
#   by.ID: If TRUE, the mean scores are computed within unique IDs
# Output:
#   a data.frame with averaged scores
# Usage example:
#   testscores <- data.frame(ID = rep(c("A", "B"), c(3, 2)),
#                            SE = 1:5,
#                            DS = rnorm(5))
#   mean_score(testscores, by.ID = FALSE)
#   mean_score(testscores, by.ID = TRUE)

suppressPackageStartupMessages(library(tidyverse)) # Required library
mean_score <- function(data, by.ID = FALSE) {
  if (by.ID) {
    out <- data %>%
      group_by(ID) %>%
      summarise(SE = ifelse(is.null(data$SE), NA, mean(SE)),
                DS = ifelse(is.null(data$DS), NA, mean(DS)),
                Brier = ifelse(is.null(data$Brier), NA, mean(Brier))) %>%
      ungroup() %>%
      select(ID, SE, DS, Brier)
  } else {
    out <- data %>%
      summarise(SE = ifelse(is.null(data$SE), NA, mean(SE)),
                DS = ifelse(is.null(data$DS), NA, mean(DS)),
                Brier = ifelse(is.null(data$Brier), NA, mean(Brier))) %>%
      select(SE, DS, Brier)
  }
  as.data.frame(out)
}
