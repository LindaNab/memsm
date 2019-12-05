#generate data deterministically
nobs <- 1e3
alpha <- 1
beta  <- 1
gamma <- 2
lambda <- 0.6 #P(L)
p_0 <- 0.10 #P(L*|L=0)
p_1 <- 0.9 #P(L*|L=1)
pi_0 <- 0.5 #P(A|L=0)
pi_1 <- 0 #P(A|L=1)
L <- c(rep(0, nobs*(1-lambda)), rep(1, nobs*lambda))
mean(L)
A <- numeric(nobs)
#n_al
n00 <- (1-pi_0)*(1-lambda)*nobs
n01 <- (1-pi_1)*lambda*nobs
n10 <- pi_0*(1-lambda)*nobs
n11 <- pi_1*lambda*nobs
A[L==0] <- c(rep(0, n00), rep(1, n10))
A[L==1] <- c(rep(0, n01), rep(1, n11))
mean(A[L==0])
mean(A[L==1])
#calculate omega := E[A] = P(A|L=0)*P(L=0) + P(A|L=1)*P(L=1)
omega <- pi_0 * (1-lambda) + pi_1 * lambda
mean(A) 
#calculate ell := E[L^*] = P(L*=1|L=0)*P(L=0) + P(L*=1|L=1)*P(L=1)
ell <- p_0 * (1-lambda) + p_1 * lambda
Ls <- numeric(nobs)
#n_lstar_a_l
n000 <- (1-p_0)*(1-pi_0)*(1-lambda)*nobs #180
n100 <- p_0*(1-pi_0)*(1-lambda)*nobs #20
n010 <- (1-p_0)*pi_0*(1-lambda)*nobs #180
n110 <- p_0*pi_0*(1-lambda)*nobs #20
n001 <- (1-p_1)*(1-pi_1)*(lambda)*nobs #42
n101 <- p_1*(1-pi_1)*(lambda)*nobs #378
n011 <- (1-p_1)*pi_1*(lambda)*nobs #18
n111 <- p_1*pi_1*(lambda)*nobs #162
n000; n100; n010; n110; n001; n101; n011; n111
Ls[A==0 & L==0] <- c(rep(0, each = n000), rep(1, times = n100))
Ls[A==1 & L==0] <- c(rep(0, times = n010), rep(1, times = n110))
Ls[A==0 & L==1] <- c(rep(0, times = 3), rep(1, times = 27)) #weird behaviour of rep, fill out times manually
Ls[A==1 & L==1] <- c(rep(0, times = 57), rep(1, times = 513)) #weird behaviour of rep, fill outtimes manually
mean(Ls)
mean(Ls[L==0])
mean(Ls[L==1])

#outcome model
Y <- alpha + beta * A + gamma * L
data <- data.frame(Y, A, L, Ls)

#pi_star = P(A|L_star=l_star)
pi_star <- function(l_star){
  term0 <- pi_0 * ((1-p_0)^(1-l_star)*p_0^l_star*(1-lambda)) / ((1-ell)^(1-l_star)*ell^l_star)
  term1 <- pi_1 * ((1-p_1)^(1-l_star)*p_1^l_star*lambda) / ((1-ell)^(1-l_star)*ell^l_star)
  return(term0 + term1)
}
pi_star_0 <- pi_star(0)
mean(A[Ls==0])
pi_star_1 <- pi_star(1)
mean(A[Ls==1])

#phi = P(L|A=a,L_star=l_star)
phi <- function(a, l_star){
  out <- (lambda*(1-pi_1)^(1-a)*pi_1^a*(1-p_1)^(1-l_star)*p_1^l_star) / 
    ((1-pi_star(l_star))^(1-a)*pi_star(l_star)^a*(1-ell)^(1-l_star)*ell^l_star)
  return(out)
}
phi_00 <- phi(0,0)
mean(L[A==0 & Ls ==0])
phi_10 <- phi(1,0)
mean(L[A==1 & Ls ==0])
phi_01 <- phi(0,1)
mean(L[A==0 & Ls ==1])
phi_11 <- phi(1,1)
mean(L[A==1 & Ls ==1])

bias_cm <- lm(Y ~ A + Ls, data = data)$coef[2] - beta
temp <- (pi_star_1*ell*(1-omega)-pi_star_1*(pi_star_1 - pi_star_0)*ell*(1-ell)) / (omega*(1-omega)-(pi_star_1-pi_star_0)^2*ell*(1-ell))
bias_condmod <- gamma * (phi(1,0) - phi(0,0)) * (1 - temp) + gamma * (phi(1,1) - phi(0,1)) * temp

require(ipw)
require(survey)
temp <- ipwpoint(exposure = A, 
                 family = "binomial", 
                 link = "logit", 
                 denominator = ~ Ls,
                 numerator = ~ 1,
                 data = data)
msm <- svyglm(Y ~ A, 
              design = svydesign(~ 1, 
                                 weights = ~ temp$ipw.weights, 
                                 data = data))

bias_msm <- msm$coefficients[2] - beta

bias_margstrmod <- gamma * (phi(1,0) - phi(0,0)) * (1 - ell) + gamma * (phi(1,1) - phi(0,1)) * ell


#pi_1
(pi_star_1 * prev_Ls - pi_0 * p_0 * (1 - prev_L)) / (p_1 * prev_L)
#pi_0
(pi_star_0 * (1 - prev_Ls) - pi_1 * (1 - p_1) * prev_L) / ((1 - p_0) * (1 - prev_L))
#pi_1
(pi_star_1 * prev_Ls - (pi_star_0 * (1 - prev_Ls) - pi_1 * (1 - p_1) * prev_L) / ((1 - p_0) * (1 - prev_L)) * p_0 * (1 - prev_L)) / (p_1 * prev_L)
#splits termen
(pi_star_1 * prev_Ls - pi_star_0 * (1 - prev_Ls) * (p_0 / (1 - p_0))) / (p_1 * prev_L) 
(pi_1 * (1 - p_1) * prev_L * (p_0/ (1 - p_0))) / (p_1 * prev_L)
#term 1 different
(pi_1 * (1 - p_1) * p_0) / ((1 - p_0) * p_1) 


#this is code used to calc bias in cm
Q <- A * Ls
lm(Q ~ A + Ls)$coef[2]
u1 <- (var(Ls)*cov(A, Q) - cov(A, Ls)*cov(Ls, Q)) / (var(Ls)*var(A) - cov(A, Ls)^2)

#calculate l := E[L^*] = P(L*=1|L=0)*P(L=0) + P(L*=1|L=1)*P(L=1)
l <- p_0 * (1 - lambda) + p_1 * lambda
mean(Ls)
#calculate omega := E[A] = P(A|L=0)*P(L=0) + P(A|L=1)*P(L=1)
omega <- pi_0 * (1-lambda) + pi_1 * (1-lambda)
mean(A)
#var(Ls) = l*(1-l)^2 + (1-l)*l^2
varLs <- l*(1-l)
var(Ls)*(nobs-1)/nobs
#var(A) = omega(1-omega)
varA <- omega*(1-omega)
var(A)*(nobs-1)/nobs

#P(L=ind_l|L^*=ls) (needed to calc pi_star_lstar = sum_l(P(A=1|L=l) X P(L=l|L*=l*))
#so this function is to calculate the last bit of the formula)
#note that l is here 'ell' or P(L*)
P_L_Ls <- function(p_0, p_1, lambda, l, ind_l, ls){
  out <- ( ((1-p_0)^(1-ls)*p_0^(ls))^(1-ind_l) * ((1-p_1)^(1-ls)*p_1^(ls))^ind_l * (1-lambda)^(1-ind_l) * lambda^ind_l ) / ((1-l)^(1-ls)*l^ls)
  return(out) 
}
#pi*_l^* = P(A=1|L^*=l^*)
pi_s_ls <- function(p_0, p_1, lambda, l, pi_0, pi_1, ls){
  pi_0 * P_L_Ls(p_0, p_1, lambda, l, 0, ls) + pi_1 * P_L_Ls(p_0, p_1, lambda, l, 1, ls) 
}
pi_s_0 <- pi_s_ls(p_0, p_1, lambda, l, pi_0, pi_1, 0) #dit klopt niet :(
mean(A[Ls==0])
pi_s_1 <- pi_s_ls(p_0, p_1, lambda, l, pi_0, pi_1, 1)
mean(A[Ls==1])


#Cov(A*L^*, A)
lm(Q ~ A)
mean(Ls[A==1])
mean(A[Ls==1])*mean(Ls)/mean(A)
pi_s_1 * l / omega
covQA <- pi_s_1 * l / omega * varA
cov(Q, A)

#Cov(A*L^*, L^*)
covQLs <- pi_s_ls(p_0, p_1, lambda, l, pi_0, pi_1, 1) * varLs
cov(Q, Ls)

#Cov(A, L^*)
covALs <- (pi_s_ls(p_0, p_1, lambda, l, pi_0, pi_1, 1) - 
             pi_s_ls(p_0, p_1, lambda, l, pi_0, pi_1, 0)) * varLs
cov(A, Ls)
