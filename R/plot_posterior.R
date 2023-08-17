#' Creates a plot for posterior probability
#'
#' Creates a ggplot showing over time development of the posterior probability that each arm is best from an adaptive experiment.
#'
#' @param adapt_matrix (default = \code{NULL}) result of \code{run_adapt()}
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @return \code{plot_posterior} returns ggplot-based plot
#' @references Lee and Green. (2023+). Discovering optimal ballot wording using adaptive survey design.
#' @export

plot_posterior <- function(adapt_matrix = NULL){

  arms <- ncol(adapt_matrix$m_posterior)
  period <- nrow(adapt_matrix$m_posterior) - 1
  post <- NULL
  for (j in 1:arms){
    post <- rbind(post, data.frame(period = 1:(period+1), arm = rep(paste('Arm', j), period+1), posterior = adapt_matrix$m_posterior[, j]))
  }
  post <- post[post$period <= period,]
  post$p_text <- factor(post$arm, levels = post[post$period == period,][order(post$posterior[post$period == period], decreasing = TRUE), 'arm'])
  post_text  <- post[post$period == period,]

  ggplot(data = post, aes(x = post$period, y = post$posterior, color = post$p_text)) +
    geom_line() +
    coord_cartesian(xlim = c(0.5, (period + 2)),  ylim = c(0, 1), clip = 'off') +
    scale_x_continuous(breaks = c(seq(1, period, 1))) +
    scale_y_continuous(breaks = c(seq(0, 1, 0.1))) +
    scale_colour_manual(
      name = '',
      breaks = levels(post$p_text),
      labels = levels(post$p_text),
      values = paste0('gray', seq(0, 99, floor(99/arms)))) +
    ylab('Posterior probability of being the best arm') + xlab('Batch number') +
    geom_text_repel(data = post_text, aes(x = post_text$period, y = post_text$posterior, color = post_text$p_text, label = post_text$p_text),
                    nudge_x = 10, hjust = 1, segment.size = .2,
                    seed = 343, direction = 'y', size = 3.5) +
    theme_bw() + theme(legend.position = 'none',
                       panel.grid.minor = element_blank(),
                       plot.caption = element_text(hjust = 0)) -> p_posterior

  return(p_posterior)
}
