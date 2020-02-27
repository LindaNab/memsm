# makes grid to be used in the sensitivity analysis
make_grid <- function(min_p0, max_p0, min_p1, max_p1, 
                     pi_star_1, pi_star_0, ell, gamma_star){
  size <- 20
  size_p0 <- size
  size_p1 <- size
  # is a range selected or one point?
  if (min_p0 == max_p0 & min_p1 != max_p1){
    size_df <- size
    size_p0 <- 1
  }
  else if (min_p0 != max_p0 & min_p1 == max_p1){
    size_df <- size
    size_p1 <- 1
  }
  else if (min_p0 == max_p0 & min_p1 == max_p1){
    size_df <- 1
    size_p0 <- 1
    size_p1 <- 1
  }
  else {size_df <- size * size}
  grid <- data.frame("p0_unif" = numeric(size_df), 
                     "p1_unif" = numeric(size_df), 
                     "p0_trap" = numeric(size_df), 
                     "p1_trap" = numeric(size_df), 
                     "p0_trian" = numeric(size_df), 
                     "p1_trian" = numeric(size_df), 
                     "pi_star_0" = numeric(size_df), 
                     "pi_star_1" = numeric(size_df),
                     "ell" = numeric(size_df), 
                     "gamma_star" = numeric(size_df))
  grid$pi_star_0 <- pi_star_0
  grid$pi_star_1 <- pi_star_1
  grid$ell <- ell
  grid$gamma_star <- gamma_star
  # uniform
  p0_sample <- seq(from = min_p0, to = max_p0, length.out = size_p0)
  p1_sample <-  seq(from = min_p1, to = max_p1, length.out = size_p1)
  grid_unif <- expand.grid("p_0" = p0_sample, "p_1" = p1_sample)
  grid$p0_unif <- grid_unif$p_0
  grid$p1_unif <- grid_unif$p_1
  # trapezidiol
  u <- seq(from = 0, to = 1, length.out = size)
  mod_lowp0 <- ( min_p0 + ( (max_p0 - min_p0) * 1 / 3 ) ) * 100
  mod_upp0 <- ( min_p0 + ( (max_p0 - min_p0) * 2 / 3 ) ) * 100
  mod_lowp1 <- ( min_p1 + ( (max_p1 - min_p1) * 1 / 3 ) ) * 100
  mod_upp1 <- ( min_p1 + ( (max_p1 - min_p1) * 2 / 3 ) ) * 100
  if (min_p0 != max_p0){
    trap0 <- sapply(u, sample_trap, min_p0*100, mod_lowp0, mod_upp0, max_p0*100)
  }
  else trap0 <- min_p0 * 100
  if (min_p1 != max_p1){
    trap1 <- sapply(u, sample_trap, min_p1*100, mod_lowp1, mod_upp1, max_p1*100)
  }
  else trap1 <- min_p1 * 100
  grid_trap <- expand.grid("p_0" = trap0/100, "p_1" = trap1/100)
  grid$p0_trap <- grid_trap$p_0
  grid$p1_trap <- grid_trap$p_1
  # triangular
  mod_lowp0 <- ( min_p0 + ( (max_p0 - min_p0) * 1 / 2 ) ) * 100
  mod_upp0 <- mod_lowp0
  mod_lowp1 <- ( min_p1 + ( (max_p1 - min_p1) * 1 / 2 ) ) * 100
  mod_upp1 <- mod_lowp1
  if (min_p0 != max_p0){
    trian0 <- sapply(u, sample_trap, min_p0*100, mod_lowp0, mod_upp0, max_p0*100)
  }
  else trian0 <- min_p0 * 100
  if (min_p1 != max_p1){
    trian1 <- sapply(u, sample_trap, min_p1*100, mod_lowp1, mod_upp1, max_p1*100)
  }
  else trian1 <- min_p1 * 100
  grid_trian <- expand.grid("p_0" = trian0/100, "p_1" = trian1/100)
  grid$p0_trian <- grid_trian$p_0
  grid$p1_trian <- grid_trian$p_1
  grid
}

sample_trap <- function(u, min, mod_low, mod_up, max){
  s <- (min + mod_low + u * (max + mod_up - min - mod_low)) / 2
  if (s > mod_up){
    trap <- max - sqrt(2 * (max - mod_up) * (s - mod_up))
  } else if (s < mod_low){
    trap <- min + sqrt((mod_low - min) * u * (max + mod_up - min - mod_low))
  }
  else trap <- s
  return(trap)
}
