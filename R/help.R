sim_out <- function(n_total = NA, periods = NA, arms = NA, probs, static = FALSE, n_first = NA, floor = NA){

  # account for first batch size if not all equal
  if(is.na(n_first)){
    n_first <- n <- n_total/periods
  } else {
    n <- (n_total-n_first)/(periods-1)
  }

  xmat <- nmat <- ppmat <- matrix(NA, ncol = arms, nrow = periods)

  i <- 1
  nmat[i,] <- table(complete_ra(N = n_first, prob_each = rep(1/arms, arms)))

  # outcomes: sample from a binomial distribution according to each arm's true p
  xmat[i,] <- mapply(rbinom, n = 1, size = nmat[i,], prob = probs)
  # use x, n, to get posterior probability
  tsprob <- best_binomial_bandit_sim(xmat[i,], nmat[i,])
  # apply floor, if appropriate
  if (!is.na(floor)){
    ppmat[i,] <- impose_floor(tsprob, floor)
  } else{
    ppmat[i,] <- tsprob
  }

  # For subsequent arms, sampling is proportionate to posterior probability (thompson sampling)
  if(periods>i){
    repeat{
      i <-  i+1
      if(static == TRUE){
        newvals <- table(complete_ra(N = n, prob_each = rep(1/arms, arms)))
      } else {
        newvals <- table(complete_ra(N = n, prob_each = ppmat[(i-1),]))
      }
      nmat[i,] <- nmat[(i-1),] + newvals # cumulative
      xmat[i,] <- xmat[(i-1),] + mapply(rbinom, n = 1, size = newvals,
                                        prob = probs)
      ppmat[i,] <- best_binomial_bandit_sim(xmat[i,], nmat[i,])

      if (i == periods) break
    }}
  out <- list(nmat = nmat, xmat = xmat, ppmat = ppmat)
  return(out)
}

ipw <- function(vals, static = FALSE, se_type='HC2'){

  bbmat <- melt(rbind(vals$nmat[1,],
                      diff(vals$nmat)))
  bbmat$success <- melt(rbind(vals$xmat[1,],
                              diff(vals$xmat)))$value

  yvals <- as.vector(unlist(apply(bbmat, 1, function(x)
    c(rep(0, x['value']-x['success']), rep(1, x['success']) ))))

  zvals <- as.factor(unlist(apply(bbmat, 1, function(x)
    c(rep(x['Var2'], x['value'])))))

  if(static == TRUE){ # static, no control

    wvals <- rep(1/length(unique(zvals)), length(yvals))
    lmfit <- lm_robust(yvals ~ -1 + zvals, weights = 1/wvals, se_type = se_type)

  }
  else if(static == FALSE){ # adaptive, no control

    bbmat$weight <- melt(rbind(rep(1/ncol(vals$ppmat), ncol(vals$ppmat)),
                               vals$ppmat[1:( nrow(vals$ppmat)-1),]))$value
    wvals <- unlist(apply(bbmat, 1, function(x)
      c(rep(x['weight'], x['value']))))

    lmfit <- lm_robust(yvals ~ -1 + zvals, weights = 1/wvals, se_type = se_type)
  }
  return(lmfit)
}

impose_floor <- function(a, amin){
  new <- pmax(a, amin)
  total_slack <- sum(new) - 1
  individual_slack <- new - amin
  c <- total_slack / sum(individual_slack)
  new - c * individual_slack
}
