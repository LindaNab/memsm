# bias in ATE estimated in marginal structural model using ipw
calc_bias_msm <- function(p_1, p_0, pi_0, pi_1, gamma, lambda){
  # P(L_star=1)
  ell <- p_0 * (1 - lambda) + p_1 * lambda
  # pi_star = P(A|L_star=l_star)
  pi_star <- function(l_star){
    temp0 <- pi_0 * ((1-p_0)^(1-l_star)*p_0^l_star*(1-lambda)) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    temp1 <- pi_1 * ((1-p_1)^(1-l_star)*p_1^l_star*lambda) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    return(temp0 + temp1)
  }
  #phi = P(L|A=a,L_star=l_star)
  phi <- function(a, l_star){
    out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
      ((1-pi_star(l_star))^(1-a)*pi_star(l_star)^a*(1-ell)^(1-l_star)*ell^l_star)
    return(out)
  }
  bias_msm <- gamma * (phi(1,0) - phi(0,0)) * (1 - ell) + 
    gamma * (phi(1,1) - phi(0,1)) * ell
  return(bias_msm)}
# bias in ATE estimated in a conditional model
calc_bias_cm <- function(p_1, p_0, pi_0, pi_1, gamma, lambda){
  # omega = P(A)
  omega <- pi_0 * (1-lambda) + pi_1 * lambda
  # ell = P(L^*)
  ell <- p_0 * (1-lambda) + p_1 * lambda
  # pi_star = P(A|L_star=l_star)
  pi_star <- function(l_star){
    term0 <- pi_0 * ((1-p_0)^(1-l_star)*p_0^l_star*(1-lambda)) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    term1 <- pi_1 * ((1-p_1)^(1-l_star)*p_1^l_star*lambda) / 
      ((1-ell)^(1-l_star)*ell^l_star)
    return(term0 + term1)
  }
  # phi = P(L|A=a,L_star=l_star)
  phi <- function(a, l_star){
    out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
     ((1-pi_star(l_star))^(1-a)*pi_star(l_star)^a*(1-ell)^(1-l_star)*ell^l_star)
    return(out)
  }
  #bias
  temp <- ( pi_star(1) - pi_star(1) ^ 2 ) / 
    ( pi_star(1) * ell + pi_star(0) * ( 1 - ell ) - 
        pi_star(1) ^ 2 * ell - pi_star(0) ^ 2 * ( 1 - ell) )
  bias_cm <- gamma * (phi(1,0) - phi(0,0)) * (1 - ell * temp) + 
    gamma * (phi(1,1) - phi(0,1)) * ell * temp
  return(bias_cm)
}