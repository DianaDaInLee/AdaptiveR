#' Creates a plot for posterior probability
#'
#' Creates a ggplot showing over time development of the posterior probability that each arm is best from an adaptive experiment.
#'
#' @param adapt_matrix (default = \code{NULL}) result of \code{run_adapt()}
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @return \code{plot_posterior} returns ggplot-based plot
#' @references Lee and Green. (2023+). Discovering optimal ballot wording using adaptive survey design.
#' @export

plot_posterior <- function(adapt_matrix = NULL){

  post <- as.data.frame(adapt_matrix$m_posterior)
  post$period <- 1:nrow(post)
  post <- pivot_longer(data = post, cols = V1:V11)
  post$name <- gsub('V', 'Arm ', post$name)
  return(post)
  # fct_ord <- as.vector(post[post$period == nrow(post),] %>% arrange(-value) %>% select(name))
  # post$p.cat = factor(post$name, levels = fct_ord$name)
  #
  # post_t <- post[post$period == max(post$period),] %>%
  #   group_by(period, name, p.cat) %>%
  #   summarize(value = mean(value))
  #
  # ggplot(data = post, aes(x = period, y = value, group = name, color = p.cat, linetype = p.cat)) +
  #   geom_line() +
  #   coord_cartesian(xlim = c(0.5, (10 + 4)),  ylim = c(0, 0.4), clip = 'off') +
  #   scale_x_continuous(breaks = c(seq(1, 10, 1))) +
  #   scale_y_continuous(breaks = c(seq(0, 1, 0.1))) +
  #   scale_colour_manual(
  #     name = 'Position in last period',
  #     breaks = levels(post$p.cat),
  #     labels = levels(post$p.cat),
  #     # values = RColorBrewer::brewer.pal(11, 'Paired')) +
  #     values = c('red', 'blue', 'forestgreen', 'orange', paste0('gray', seq(0, 99, floor(90/7))))) +
  #   scale_linetype_manual(name = '',
  #                         values = c(rep('solid', 4), rep('longdash', 3), rep('dotdash', 2), rep('twodash', 2))) +
  #   ylab('Posterior probability of being the best arm') + xlab('Batch number') +
  #   geom_text_repel(data = post_t, aes(label = p.cat), nudge_x = 10, hjust = 1, segment.size = .2,
  #                   seed = 343, direction = 'y', size = 3) +
  #   theme_bw() + theme(legend.position = 'none',
  #                      panel.grid.minor = element_blank(),
  #                      plot.caption = element_text(hjust = 0)) -> p_posterior

  # return(plot_posterior)
}
