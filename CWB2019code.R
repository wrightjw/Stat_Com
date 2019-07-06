suppressPackageStartupMessages(library(tidyverse)) # Required library

# logit: Compute a logit transformation
# Input:
#   x : scalar or vector of values between 0 and 1
# Output:
#   The logit transformation of x

logit <- function(x) {
  log(x) - log(1-x)
}

# ilogit: Compute the inverse logit transformation
# Input:
#   x : scalar or vector
# Output:
#   The inverse logit transformation of x; values between 0 and 1

ilogit <- function(x) {
  1 / (1 + exp(-x))
}


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

# score_brier: Compute the Brier score
#   Note: this implementation uses the prediction probabilities
#         and event indicator values directly, unlike the CWA version!
# Input:
#   prob : vector with event prediction probabilities
#   event : data vector of event indicators
# Output:
#   a vector of scores

score_brier <- function(prob, event) {
  (event - prob)^2
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


# Implementation of an answer from Coursework A:
prediction <- function(newdata, fit) {
  pred <- predict(fit, newdata = newdata, se.fit = TRUE)
  # Compute an estimate of the prediction standard deviation:
  pred_sigma <- sqrt(pred$se.fit^2 + sigma(fit)^2)
  data.frame(mu = pred$fit,
             sigma = pred_sigma,
             prob = pnorm(0, mean = pred$fit, sd = pred_sigma))
}

# Extra helper function to avoid duplicating code in Coursework A.
# Input:
#   newdata: A data.frame with data to be used for model assessment
#   pred: output from the prediction() function
#   response_name: A ""-quoted string with the name of the response
#                  variable in the data objects
# Output:
#   A data.frame with computed scores
collect_scores <- function(newdata, pred, response_name) {
  # Deviation for the coursework: score_brier is defined differently, see above.
  score <- data.frame(
    SE = score_se(pred, newdata[, response_name]),
    DS = score_ds(pred, newdata[, response_name]),
    Brier = score_brier(pred$prob, (newdata[, response_name] <= "0")))
  if (is.null(newdata$ID)) {
    score
  } else {
    cbind(data.frame(ID = newdata$ID), score)
  }
}


# Estimate a model on training data and return validation scores
# Input:
#   train: A data.frame with data to be used for model estimation
#   valid: A data.frame with data to be used for model assessment
#   formula: The formula defining the model estimable with lm()
# Output:
#   The prediction scores (SE and DS for y predictions,
#   and Brier for y <= 0 events) for the validation data.
train_and_validate <- function(train, valid, formula) {
  fit <- lm(formula, data = train)
  pred <- prediction(newdata = valid, fit = fit)
  # This seems to be the most reliable approach to extract the response name:
  response_name <- all.vars(update(formula, . ~ 1))
  collect_scores(valid, pred, response_name)
}

# From lab 6 solutions ----

# See the lab06 Q7 solution for information about this function
make_formula <- function(max_order) {
  form <-
    paste0("Value ~ Longitude + Latitude",
           " + Elevation + I(DecYear-Year)")
  if (max_order > 0) {
    form <-
      paste0(form,
             paste0(" + I(cos(2*pi*DecYear*", seq_len(max_order), "))",
                    collapse = ""),
             paste0(" + I(sin(2*pi*DecYear*", seq_len(max_order), "))",
                    collapse = ""))
  }
  as.formula(form)
}

cvk_define_splits <- function(N, K) {
  if (K > N) {
    stop("K-fold cross validation cannot have K > N")
  }
  sub_size <- ceiling((1:K) * N / K) - ceiling((0:(K-1)) * N / K)
  sample(rep(1:K, times = sub_size), size = N, replace = FALSE)
}


cvk_split <- function(data, splits, k) {
  list(train = data[splits != k, , drop = FALSE],
       valid = data[splits == k, , drop = FALSE])
}


cvk_do_all <- function(data, splits, formula) {
  K <- max(splits)
  S_hat <- data.frame(SE = numeric(K),
                      DS = numeric(K),
                      Brier = numeric(K))
  for (k in 1:K) {
    data_k <- cvk_split(data, splits, k)
    scores_k <- train_and_validate(train = data_k$train,
                                   valid = data_k$valid,
                                   formula = formula)
    S_hat[k,] <- mean_score(scores_k, by.ID = FALSE)[c("SE", "DS", "Brier")]
  }
  S_hat
}

# read_TMIN_data: Read the Scottish TMIN data
# Input:
#   path : Optional path to the data file location
# Output:
#   A list with elements "train" and "test" where
#   the "obs" data are in "train" and the "test" data in "test".

read_TMIN_data <- function(path = NULL) {
  file_obs <- "TMINallobs.csv"
  file_test <-"TMINalltest.csv"
  if (!is.null(path)) {
    file_obs <- file.path(path, file_obs)
    file_test <- file.path(path, file_test)
  }
  list(train = read.csv(file = file_obs,
                        header = TRUE,
                        stringsAsFactors = FALSE),
       test = read.csv(file = file_test,
                        header = TRUE,
                        stringsAsFactors = FALSE))
}


# boot_resample: Draw a random sample of the same size as the input
# Input:
#   data: A data.frame
# Output:
#   A data.frame.

boot_resample <- function(data) {
  data[sample(nrow(data), size = nrow(data), replace = TRUE), , drop = FALSE]
}


# data_list_resample: Resample with replacement across data.frame:s in a list
# Input:
#   data_list : list of data.frame:s
# Output:
#   A list of data.frame:s of the same size and the same names
#   as in the input list.

data_list_resample <- function(data_list) {
  resampled <- boot_resample(do.call(rbind, data_list))
  N <- vapply(data_list, function(x) nrow(x), 1L)
  idx_offset <- c(0, cumsum(N))
  new_data <-
    lapply(seq_along(data_list),
           function(x) resampled[idx_offset[x] + seq_len(N[x]), ,
                                 drop = FALSE])
  names(new_data) <- names(data_list)
  new_data
}

# cwb4_scores: Function for Coursework B Question 4.

cwb4_scores <- function(data_list, K) {
  # Step 1:
  the_splits <- cvk_define_splits(nrow(data_list$train), K)
  # Step 2:
  scores_cvk <- cvk_do_all(data_list$train,
                           the_splits,
                           make_formula(2))
  # Step 3:
  scores_test <-
    mean_score(
      train_and_validate(train = data_list$train,
                         valid = data_list$test,
                         formula = make_formula(2)),
      by.ID = FALSE)[c("SE", "DS", "Brier")]
  # Step 4:
  list(cvk_mean = colMeans(scores_cvk),
       cvk_std_err = apply(scores_cvk, 2, sd) / sqrt(K),
       test = scores_test)
}
