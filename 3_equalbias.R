# Term between 'curly brackets' in bias in conditional model
Brackets <- function(pi1star, pi0star, omega, ell){
  ( pi1star * (1 - omega) - pi1star * (pi1star - pi0star) * (1 - ell) ) / 
    ( omega * (1 - omega) - (pi1star - pi0star) ^ 2 * ell * (1 - ell)) 
}
# Term between 'curly brackets' rewritten
Brackets2 <- function(pi1star, pi0star, ell){
  ( pi1star - pi1star ^ 2 ) / 
    ( pi1star * ell + pi0star * ( 1 - ell ) - pi1star ^ 2 * ell - pi0star ^ 2 * ( 1 - ell) )
}
pi1star <- 0.6
pi0star <- 0.3
ell <- 0.2
Brackets2(1 - pi0star, pi0star, ell)
Brackets2(pi1star, 1 - pi1star, ell)
omega <- pi1star * ell + pi0star * (1 - ell)
Brackets(pi1star, pi0star, omega, ell)
Brackets2(pi1star, pi0star, ell)

# Maximum for pi1star
plot(seq(0, 1, length.out = 100), Brackets2(seq(0, 1, length.out = 100), pi0star, ell))
abline(h = 1)
points(pi0star, 1, col = "red")
optimize(Brackets2, c(0, 1), pi0star = pi0star, ell = ell, maximum= TRUE)
# Maximum is always larger than 1, and equals 1 if pi

# Minimum for pi0star
plot(seq(0, 1, length.out = 100), Brackets2(pi1star, seq(0, 1, length.out = 100), ell))
abline(h = 1)
optimize(Brackets2, c(0, 1), pi1star = pi1star, ell = ell, maximum= F)

abc <- function(a, b, c){
  D <- b^2 - 4 * a * c
  if ( D > 0 )
  outPos <- ( - b + sqrt(D) ) / (2 * a)
  outNeg <- ( - b - sqrt(D) ) / (2 * a)
  c(outPos, outNeg)
}

