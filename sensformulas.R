# Sensitivity analysis----------------------------------------------------------
# Parameters in data
# in a sens analyse, you must estimate the following parameters in your data
# ell = P(L*)
# pi_star_1 = P(A|L*=1)
# pi_star_0 = P(A|L*=0)
# omega = P(A)
# gamma_star = coefficient from the linear regression modelling 
# E[Y|A,L,AL] = alpha_star + beta_star * A + gamma_star * L* + delta * A * L
# Further, one guesses values for p_0 and p_1
# From these parameters, we calculate lambda and pi_0 and pi_1:
calc_lambda <- function(ell, p_0, p_1){
  if (p_0 != p_1){
    (ell - p_0) / (p_1 - p_0)}
  else p_0
} 
calc_pi_1 <- function(pi_star_0, pi_star_1, ell, p_0, p_1){
  if (p_0 == 1 & p_1 != 0){
    pi_1 = pi_star_0
  } else if (p_1 == 0 & p_0 != 0){
    term1 <- pi_star_0 * (1 - ell) - pi_star_1 * (1 - p_0) * 
      (1 - calc_lambda(ell, p_0, p_1))
    term2 <- calc_lambda(ell, p_0, p_1)
    pi_1 <- term1 / term2
  } else {
    term1 <- (pi_star_1 * ell - pi_star_0 * (1 - ell) * (p_0 / (1 - p_0))) /
      (p_1 * calc_lambda(ell, p_0, p_1))
    term2 <- ((1 - p_0) * p_1) / ((1 - p_0) * p_1 - (1 - p_1) * p_0)
    pi_1 <- term1 * term2
  }
  return(pi_1)
}
calc_pi_0 <- function(pi_star_0, pi_star_1, ell, p_0, p_1){
  if (p_0 == 1 & p_1 != 0){
    term1 <- pi_star_0 * p_1 * lambda - pi_star_1 * ell
    term2 <- 1 - calc_lambda(ell, p_0, p_1)
    pi_0 <- term1 / term2
  }
  else if (p_1 == 0 & p_0 != 0){
    pi_0 = pi_star_1
  }
  else {
    pi_0 <- (pi_star_0 * (1 - ell) - calc_pi_1(pi_star_0, pi_star_1, ell, p_0, p_1) * 
       (1 - p_1) * calc_lambda(ell, p_0, p_1)) / 
       ((1 - p_0) * (1 - calc_lambda(ell, p_0, p_1)))
  }
  return(pi_0)
}
# Using these, we calculate phi_{al*} = P(L=1|A=a, L*=l*)
calc_phi <- function(pi_star_0, pi_star_1, ell, p_0, p_1, a, l_star){
  pi_0 <- calc_pi_0(pi_star_0, pi_star_1, ell, p_0, p_1)
  pi_1 <- calc_pi_1(pi_star_0, pi_star_1, ell, p_0, p_1)
  lambda <- calc_lambda(ell, p_0, p_1)
  #select correct pi_star
  pi_star <- function(pi_star_0, pi_star_1, l_star){
    if(l_star == 0) out <- pi_star_0
    else out <- pi_star_1
    return(out)}
  out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
    ((1-pi_star(pi_star_0, pi_star_1, l_star))^(1-a) * 
       pi_star(pi_star_0, pi_star_1, l_star)^a*(1-ell)^(1-l_star)*ell^l_star)
  return(out)
}
# And gamma = gamma_star / (phi_01 - phi00)
calc_gamma <- function(gamma_star, pi_star_0, pi_star_1, ell, p_0, p_1){
  gamma_star/(calc_phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 1) - 
                calc_phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 0))
}
# sensitivity analysis msm 
sensanalyse <- function(df, namep_0 = "p_0", namep_1 = "p_1"){
  # input
  pi_star_0 <- df["pi_star_0"]
  pi_star_1 <- df["pi_star_1"]
  ell <- df["ell"]
  gamma_star <- df["gamma_star"]
  p_0 <- df[namep_0]
  p_1 <- df[namep_1]
  # calculate parameters
  lambda <- calc_lambda(ell, p_0, p_1)
  pi_1 <- calc_pi_1(pi_star_0, pi_star_1, ell, p_0, p_1)
  pi_0 <- calc_pi_0(pi_star_0, pi_star_1, ell, p_0, p_1)
  phi_01 <- calc_phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 1)
  phi_00 <- calc_phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 0)
  gamma <- calc_gamma(gamma_star, pi_star_0, pi_star_1, ell, p_0, p_1)
  # calculate bias
  if (p_0 != p_1){
    bias_msm <- calc_bias_msm(p_1, p_0, pi_0, pi_1, gamma, lambda)
    bias_cm <- calc_bias_cm(p_1, p_0, pi_0, pi_1, gamma, lambda)
  }
  else{
    bias_msm <- NA
    bias_cm <- NA
  }
  return(c(bias_msm = unname(bias_msm), bias_cm = unname(bias_cm)))
}