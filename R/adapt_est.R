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

  moddata  <- tidy(mod)
  arm_n <- table(adapt_matrix$prior_survey_data$arm)
  moddata$arm_n <- NA
  for (i in 1:length(arm_n)){
    moddata$arm_n[moddata$term == paste0('arm', i)] <- arm_n[i]
  }
  moddata$term <- gsub('arm', 'Arm ', moddata$term)
  moddata$term <- fct_reorder(factor(moddata$term), moddata$estimate)
  moddata$label<- paste0(sprintf("%.3f", moddata$estimate), ' (', sprintf("%.3f", moddata$std.error), ') [', moddata$arm_n, ']')
  moddata$label<- ifelse(moddata$term == levels(moddata$term)[length(arm_n)], paste('Est (S.E.) [N]:', moddata$label), moddata$label)

  p <- ggplot(data = moddata, aes(x = moddata$estimate, y = moddata$term, xmin = moddata$conf.low, xmax = moddata$conf.high)) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    geom_point() +
    geom_errorbarh(height = 0) +
    xlab('Mean outcome') +
    ylab('Treatment Arm') +
    theme_bw() +
    geom_text(aes(label = moddata$label), nudge_y = .3, nudge_x = 0, size = 3, show.legend = FALSE) +
    theme(strip.background = element_blank(),
          axis.title.y = element_blank())

  return(list(est = mod, est_plot = p))
}
