#' Estimate inverse probability weighted mean outcome for each arm
#'
#' Estimates mean outcome using inverse probability weight and creates a coefficient plot based on the result
#'
#' @param adapt_matrix (default = \code{NULL}) result of \code{run_adapt()}
#' @import estimatr
#' @import ggplot2
#' @importFrom generics tidy
#' @importFrom forcats fct_reorder
#' @return \code{adapt_est} returns a list containing the following elements:
#'  \itemize{
#'    \item \code{est}: IPW estimates in an object of \code{estimatr::lm_robust()}.
#'    \item \code{est_plot}: coefficient plot created with \code{ggplot2::ggplot()}.
#'  }
#' @references Lee and Green. (2023+). Discovering optimal ballot wording using adaptive survey design.
#' @export

adapt_est <- function(adapt_matrix = NULL){

  adapt_matrix$prior_survey_data$arm <- factor(adapt_matrix$prior_survey_data$arm)
  mod <- lm_robust(Y ~ arm - 1, data = adapt_matrix$prior_survey_data, weights = adapt_matrix$prior_survey_data$ipw)

  data  <- tidy(mod)
  arm_n <- table(adapt_matrix$prior_survey_data$arm)
  data$arm_n <- NA
  for (i in 1:length(arm_n)){
    data$arm_n[data$term == paste0('arm', i)] <- arm_n[i]
  }
  data$term <- gsub('arm', 'Arm ', data$term)
  data$term <- fct_reorder(factor(data$term), data$estimate)
  data$label<- paste0(sprintf("%.3f", data$estimate), ' (', sprintf("%.3f", data$std.error), ') [', data$arm_n, ']')
  data$label<- ifelse(data$term == levels(data$term)[length(arm_n)], paste('Est (S.E.) [N]:', data$label), data$label)

  p <- ggplot(data, aes(x = data$estimate, y = data$term, xmin = data$conf.low, xmax = data$conf.high)) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point() +
    geom_errorbarh(height = 0) +
    xlab('Mean outcome') +
    ylab('Treatment Arm') +
    theme_bw() +
    geom_text(aes(label = data$label), nudge_y = .3, nudge_x = 0, size = 3, show.legend = FALSE) +
    theme(strip.background = element_blank(),
          axis.title.y = element_blank())

  return(list(est = mod, est_plot = p))
}
