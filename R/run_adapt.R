#' Conducts simple adaptive survey with standard Bernoulli Thompson sampling
#'
#' Allows users to iteratively update treatment assignment probability in an adaptive survey.
#' Returns all necessary statistics including success rate and posterior probabilities.
#' This function is to be used iteratively in an adaptive design, as the elements in a returned list will be iteratively updated over periods.
#'
#' @param prior_survey_data (default = \code{NULL}) data.frame containing survey results from prior period. Must contain the following two columns: \code{arm} that contains treatment arm numbers 1 through k, and \code{Y} that contains binary survey outcome (0 or 1). Ignored in the first period (i.e., when \code{current_period == 1}).
#' @param prior_adapt_matrix (default = \code{NULL}) result of \code{run_adapt()} from the previous batch. Ignored in the first period (i.e., when \code{current_period == 1}).
#' @param current_period (default = \code{NULL}) numeric value indicating the number of current batch to be implemented (i.e., total number of batches previously done + 1). For example, if a user has conducted a survey for the first batch and is looking to update treatment probability for the second batch, the value should be 2.
#' @param arms (default = \code{NULL}) Total number of arms in the survey.
#' @param periods (default = \code{NULL}) Total number of periods to be implemented in the survey.
#' @param floor (default = \code{NULL}) A numeric value between 0 and 1 indicating the minimum probability of treatment assignment imposed on each arm.
#' @param seed (default = 1004) Integer to specify seeds.
#' @importFrom stats rbeta
#' @return \code{run_adapt} returns a list containing the following elements:
#'  \itemize{
#'    \item \code{m_posterior}: matrix containing posterior probability for each arm (columns) for each period (rows).
#'    \item \code{m_trials}: matrix containing cumulative number of observations assigned to each arm (columns) for each period (rows).
#'    \item \code{m_success}: matrix containing cumulative number of successes (i.e., when outcome equals 1) for to each arm (columns) for each period (rows).
#'    \item \code{prior_survey_data}: appended survey data with batch number, treatment assignment probability, and inverse-probability weight appended.
#'  }
#' @references Lee and Green. (2023+). Discovering optimal ballot wording using adaptive survey design.
#' @export

run_adapt <- function(prior_survey_data = NULL, prior_adapt_matrix = NULL, current_period = NULL, arms = NULL, periods = NULL, floor = NULL, seed = 1004){

  # First period: set up matrices and assign equal treatment probability
  if (current_period == 1){
    if (is.null(arms)){
      stop('arms must be specified.')
    }
    if (is.null(periods)){
      stop('periods must be specified')
    }
    m_trials <- m_success <- matrix(nrow = periods, ncol = arms)
    m_posterior <- matrix(nrow = periods + 1, ncol = arms)
    m_posterior[current_period,] <- 1 / arms
    df <- NULL
  }

  # Subsequent period: update matrices and update treatment assignment probability with Thompson sampling
  if (current_period > 1){
    m_trials         <- prior_adapt_matrix[['m_trials']]
    m_success   <- prior_adapt_matrix[['m_success']]
    m_posterior <- prior_adapt_matrix[['m_posterior']]
    df          <- prior_adapt_matrix[['prior_survey_data']]

    if (is.null(arms))    arms    <- ncol(m_trials)
    if (is.null(periods)) periods <- nrow(m_trials)

    success <- trials <- NULL
    for (i in 1:arms) trials[i]  <- sum(prior_survey_data$arm == i)
    for (i in 1:arms) success[i] <- sum(prior_survey_data$arm == i & prior_survey_data$Y == 1)

    if (current_period == 2){
      m_trials[current_period - 1,]       <- trials
      m_success[current_period - 1,] <- success
    }
    if (current_period > 2){
      m_trials[current_period - 1,] <- m_trials[current_period - 2, ] + trials
      m_success[current_period - 1,] <- m_success[current_period - 2, ] + success
    }

    # Thompson sampling
    set.seed(seed)
    draws    <- replicate(10000, rbeta(arms, m_success[current_period - 1,] + 1, m_trials[current_period - 1,] - m_success[current_period - 1,] + 1))
    argmax   <- apply(draws, 2, which.max)
    probs_ts <- table(cut(argmax, 0:arms)) / 10000

    # Setting floors if specified
    if (!is.null(floor)){
      new <- pmax(probs_ts, floor)
      total_slack <- sum(new) - 1
      indiv_slack <- new - floor
      probs_ts <- new - (total_slack / sum(indiv_slack)) * indiv_slack
    }
    m_posterior[current_period, ] <- probs_ts

    # Appending Dataset
    ## batch number
    prior_survey_data$period <- current_period - 1
    ## treatment probability
    prior_survey_data$probs <- NA
    for (i in 1:arms){
      prior_survey_data$probs[prior_survey_data$arm == i]  <- m_posterior[current_period - 1, i]
    }
    ## inverse probability weight
    prior_survey_data$ipw <- 1 / prior_survey_data$probs

    df <- rbind(df, prior_survey_data)
  }

  cat(paste('Treatment assignment probabilities for period', current_period, 'are as follows:\n',
               paste(m_posterior[current_period,], collapse = ', '), '\n'))


  return(list(m_posterior = m_posterior, m_trials = m_trials, m_success = m_success, prior_survey_data = df))
}
