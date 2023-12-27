#' Simulate adaptive or static experiments while varying several parameters and produce key simulated statistics
#'
#' Simulate experiments varying arms, periods, sample size, and design and outputs simulated results including selection of best arms, RMSE, bias, coverage.
#'
#' @param probs (default = \code{NULL}) distribution of true mean outcomes for each arm
#' @param static (default = \code{FALSE}) whether the simulated experiment should be adaptive or static
#' @param periods (default = \code{10}) total number of batches in a single experiment (relevant for adaptive design only)
#' @param n (default = \code{1000}) total sample size
#' @param n_first (default = \code{NULL}) total sample size to be allotted to the first batch. If \code{NULL}, \code{n} is evenly distributed across batches
#' @param floor_rate (default = \code{0.01}) minimum sample probability to implement to ensure all arms receive non-zero sample in each batch
#' @param iter (default = \code{1000}) number of simulations to run
#' @import estimatr
#' @import tidyr
#' @import dplyr
#' @import randomizr
#' @import bandit
#' @importFrom reshape2 melt
#' @importFrom stats rbinom coef
#' @return \code{simulate} returns a list containing the following elements:
#'  \itemize{
#'    \item \code{outmat}: simulation summary statistics containing proportion of iterations that each arm is identified as the best arm, bias, RMSE, and coverage.
#'    \item \code{modmat}: IPW estimates (result of \code{estimatr::lm_robust()}) for each arm for all iterations.
#'    \item \code{postmat}: data.frame containing posterior probabilities for each arm (columns) for each period (rows) for all iterations (rows).
#'  }
#' @references Lee and Green. (2023+). Discovering optimal ballot wording using adaptive survey design.
#' @export

simulate <- function(probs = NA, static = FALSE, periods = 10, n = 1000, n_first = NA, iter = 1000, floor_rate = 0.01){
  # Housekeeping
  if (is.na(n_first) & n%%periods != 0){
    stop('`n_total` must be divisible by `periods`.')
  }
  if (!is.na(n_first) & ((n-n_first)/(periods-1))%%periods != 0){
    stop('`n_total` must be divisible by `periods`.')
  }

  K <- length(probs)
  outmat <- matrix(rep(NA, 5*iter), ncol = 5)
  d_fit  <- postmat <- NULL
  colnames(outmat) <- c('correct', 'posterior_best', 'bias', 'rmse_best', 'cov_best')

  cat("\n")
  cat("-----------------------------\n")
  cat(paste(ifelse(static == TRUE, "Static", "Adaptive"), "simulated experiment will be performed with following parameters:\n"))
  cat(paste(" Number of Arms:", K, "\n"))
  cat(paste(" Number of Periods:", periods, "\n"))
  cat(paste(" Total Sample Size:", n, "\n"))
  cat(paste(" Sample Size for First Period:", ifelse(is.na(n_first), n/periods, n_first), "\n"))
  cat(paste(" Sample Size for Remaining Periods:", ifelse(is.na(n_first), n/periods, (n-n_first)/(periods-1)), "\n"))
  cat("-----------------------------\n")
  cat("\n")

  cat(' Iteration: ')
  for(i in 1:iter){

    vals <- sim_out(probs = probs, arms = K, static = static, periods = periods, n_total = n, n_first = n_first, floor = floor_rate)
    true_val  <- max(probs)
    true_vali <- paste0('zvals', which.max(probs)) # index of max probs
    postmat   <- bind_rows(postmat,
                           bind_rows(data.frame(matrix(data = rep(1/K, K), ncol = K)) %>% mutate(period = 1, iter = i),
                                     data.frame(vals$ppmat) %>% mutate(period = 2:(1+nrow(vals$ppmat)), iter = i)))
    lm_r <- ipw(vals, static)

    # rmse
    outmat[i,'rmse_best'] <- sqrt((coef(lm_r)[true_vali]-true_val)^2)
    # coverage
    outmat[i, 'cov_best'] <- 1 * ((lm_r[['conf.low']][true_vali] < true_val) & (lm_r[['conf.high']][true_vali] > true_val))

    d_fit <- bind_rows(d_fit, tidy(lm_r) %>%
                         complete(term = paste0('zvals', 1:K)) %>%
                         mutate(ord = as.numeric(gsub('zvals', '', term))) %>%
                         arrange(ord) %>%
                         bind_cols(posterior = vals$ppmat[nrow(vals$ppmat),]) %>%
                         mutate(iter = i) %>%
                         select(-ord))

    # selected correct arm
    outmat[i,'correct'] <- (which.max(vals$ppmat[nrow(vals$ppmat),])==which.max(probs))*1
    # posterior probability of best arm
    outmat[i,'posterior_best'] <- vals$ppmat[nrow(vals$ppmat), which.max(probs)]
    # bias
    outmat[i,'bias'] <- coef(lm_r)[true_vali]-true_val
    cat(i, '..')
  }
  cat('\n')

  # summary
  est <- d_fit %>%
    left_join(data.frame(term = paste0('zvals', 1:K), true = probs)) %>%
    group_by(iter) %>%
    mutate(best = max(posterior) == posterior) %>%
    group_by(term) %>%
    summarize(true = mean(true, na.rm = T),
              best = mean(best, na.rm = T),
              bias = mean(true - estimate, na.rm = T),
              rmse = sqrt(mean((estimate - true)^2, na.rm = T)),
              coverage = mean(conf.low < true & conf.high > true, na.rm = T)) %>%
    arrange(-true) %>%
    mutate_if(is.numeric, ~sprintf('%0.3f', .x)) %>%
    mutate(term = gsub('zvals', 'Arm ', term))

  return(list(outmat = est, modmat = d_fit, postmat = postmat))
}
