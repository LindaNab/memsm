# ------------------------------------------------------------------------------
# This code creates Figure 2 + 3 in 'Quantitative bias analysis for a misclassif
# ied confounder: a comparison between marginal structural models and conditiona 
# l models' --------------------------------------------------------------------
# Author: Linda Nab, lindanab4@gmail.com ---------------------------------------
# Date of creation: 20190512 ---------------------------------------------------
# Date last change: 20200220 ---------------------------------------------------
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Set directories, source functions and load libraries
# ------------------------------------------------------------------------------
figures_dir <- "./figures" # directory where figures will be saved to
source("./biasformulas.R")
wrap_bias_msm <- function(df){
  bias_msm(p_1 = df['sens'], 
           p_0 = (1 - df['spec']),
           pi_0 = df['pi_0'], 
           pi_1 = df['pi_1'],
           gamma = df['gamma'],
           lambda = df['lambda'])
}
wrap_bias_cm <- function(df){
  bias_cm(p_1 = df['sens'], 
           p_0 = (1 - df['spec']),
           pi_0 = df['pi_0'], 
           pi_1 = df['pi_1'],
           gamma = df['gamma'],
           lambda = df['lambda'])
}

# ------------------------------------------------------------------------------
# 1. Set parameters for the figures
# ------------------------------------------------------------------------------
# Figure 2A consist of 2 plots. The x axis of the plot varies pi_1 from 0 to 
# 1, and the y axis represents bias in the ATE. 2A1) gamma = 2, pi_0 = 0.8; 
# 2A2) gamma = -2, pi_0 = 0.5
fig2A1 <- data.frame(spec = 0.95,
                     sens = 0.9, 
                     pi_0 = 0.8, 
                     pi_1 = seq(from = 0, to = 1, length.out = 100),
                     gamma = 2,
                     lambda = 0.5)
fig2A2 <- fig2A1
fig2A2$pi_0 <- 0.5
fig2A2$gamma <- -2

# Figure 2B consist of 1 plot with 2 lines. The x axis of the plot varies gamma 
# from -5 to 5, and the y axis represents bias in the ATE. 2B1) pi_1 = 0.9 
# pi_0 = 0.5; 2B2) pi_1 = 0.5, pi_0 = 0.9
fig2B1 <- data.frame(spec = 0.95,
                     sens = 0.9, 
                     pi_0 = 0.45, 
                     pi_1 = 0.9,
                     gamma = seq(from = -5, to = 5, length.out = 100),
                     lambda = 0.5)
fig2B2 <- fig2B1
fig2B2$pi_0 = 0.9
fig2B2$pi_1 = 0.5

# Figure 3 consist of 1 plot with 2 lines. The x axis of the plot varies lambda 
# from 0 to 1, and the y axis represents bias in the ATE. 3) pi_0 = 0.5 
# pi_1 = 0.75; 2B2) pi_0 = 0.25, pi_1 = 0.75
fig31 <- data.frame(spec = 0.95,
                   sens = 0.9, 
                   pi_0 = 0.5, 
                   pi_1 = 0.75,
                   gamma = 5,
                   lambda = seq(from = 0, to = 1, length.out = 100))
fig32 <- fig31
fig32$pi_0 = 0.25


png(paste0(figures_dir,"/fig2A1.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mgp = c(0, 0.25, 0), mar = c(4, 4, 4, 4))
plot(0, type = "n", xaxt = "n", yaxt = "n", frame.plot = F, ann = F,
     xlim = c(0, 1), ylim = c(-1.5, 1.5))
rect(xleft = c(0, 0), ybottom = c(-1.5, -1.5), 
     xright = c(1, 1), ytop = c(0, 0),
     col = "lightgrey", border = NA)
rect(xleft = c(0, 0), ybottom = c(-1.5, -1.5), 
     xright = c(0.8, 0.8), ytop = c(1.5, 1.5),
     density = 10, col = "grey", border = NA)
lines(c(0.8, 0.8), c(-1.5, 1.5), col = "grey")
lines(fig2A1$pi_1, apply(fig2A1, 1, wrap_bias_msm), 
      col = "black", lty = 1, lwd = 1.5)
lines(fig2A1$pi_1, apply(fig2A1, 1, wrap_bias_cm), 
      col = "black", lty = 2, lwd = 1.5)
axis(1, at = c(0, 0.8, 1), 
     label = c(0, expression(paste(pi[0], " = 0.8")), 1), 
     tck = -0.01, cex.axis = 0.75)
text(0.8, -1.75, expression(paste(pi[0], " = 0.8")))
mtext(expression(pi[1]), side = 1, line = 1.5, cex = 1.25)
axis(4, at = c(-1.5, 0, 1.5), tck = -0.01, cex.axis = 0.75)
mtext("bias in average treatment effect", side = 4, line = 1.5)
mtext("Confounder and outcome positively associated", side = 3, line = 2)
#text(0.8, 1.55, expression(paste(pi[1], " = ", pi[0])))
legend(-0.2, 2.15, xpd = T, 
       legend = c("Bias < 0", "confounder and exposure neg ass"),
                  #expression(paste(pi[1], " < ", pi[0]))), 
       fill = c("lightgrey", "grey"),
       border = c("lightgrey", "grey"),
       density = c(NA, 20),
       horiz = T,
       bty = "n",
       cex = 0.75,
       x.intersp = 0.25,
       text.width=c(0,0.25))
legend(0, 2, xpd = T, 
       legend = c("MSM-IPW", "conditional model"),
       lty = c(1, 2),
       bty = "n",
       horiz = T,
       cex = 0.75,
       x.intersp = 0.25,
       seg.len = 1.5,
       text.width = c(0, 0.3))
dev.off()

png(paste0(figures_dir,"/fig2A2.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mgp = c(0, 0.25, 0), mar = c(4, 4, 4, 4))
plot(0, type = "n", xaxt = "n", yaxt = "n", frame.plot = F, ann = F,
     xlim = c(0, 1), ylim = c(-1.5, 1.5))
rect(xleft = c(0, 0), ybottom = c(-1.5, -1.5), 
     xright = c(1, 1), ytop = c(0, 0),
     col = "lightgrey", border = NA)
rect(xleft = c(0, 0), ybottom = c(-1.5, -1.5), 
     xright = c(0.5, 0.5), ytop = c(1.5, 1.5),
     density = 10, col = "grey", border = NA)
lines(c(0.5, 0.5), c(-1.5, 1.5), col = "grey")
lines(fig2A2$pi_1, apply(fig2A2, 1, wrap_bias_msm), 
      col = "black", lty = 1, lwd = 1.5)
lines(fig2A2$pi_1, apply(fig2A2, 1, wrap_bias_cm), 
      col = "black", lty = 2, lwd = 1.5)
axis(1, at = c(0, 0.5, 1), 
     label = c(0, expression(paste(pi[0], " = 0.5")), 1), 
     tck = -0.01, cex.axis = 0.75)
mtext(expression(pi[1]), side = 1, line = 1.5, cex = 1.25)
axis(2, at = c(-1.5, 0, 1.5), tck = -0.01, cex.axis = 0.75)
mtext("bias in average treatment effect", side = 2, line = 1.5)
mtext("Confounder and outcome negatively associated", side = 3, line = 2)
#text(0.8, 1.55, expression(paste(pi[1], " = ", pi[0])))
legend(-0.2, 2.15, xpd = T, 
       legend = c("Bias < 0", "confounder and exposure neg ass"),
       #expression(paste(pi[1], " < ", pi[0]))), 
       fill = c("lightgrey", "grey"),
       border = c("lightgrey", "grey"),
       density = c(NA, 20),
       horiz = T,
       bty = "n",
       cex = 0.75,
       x.intersp = 0.25,
       text.width=c(0,0.25))
legend(0, 2, xpd = T, 
       legend = c("MSM-IPW", "conditional model"),
       lty = c(1, 2),
       bty = "n",
       horiz = T,
       cex = 0.75,
       x.intersp = 0.25,
       seg.len = 1.5,
       text.width = c(0, 0.3))
dev.off()

png(paste0(figures_dir,"/fig2B.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mgp = c(0, 0.25, 0), mar = c(4, 4, 4, 4))
plot(0, type = "n", xaxt = "n", yaxt = "n", frame.plot = F, ann = F,
     xlim = c(-5, 5), ylim = c(-1.5, 1.5))
rect(xleft = c(-5, -5), ybottom = c(-1.5, -1.5), 
     xright = c(5, 5), ytop = c(0, 0),
     col = "lightgrey", border = NA)
lines(fig2B1$gamma, apply(fig2B1, 1, wrap_bias_msm), 
      col = "black", lty = 1, lwd = 1.5)
lines(fig2B1$gamma, apply(fig2B1, 1, wrap_bias_cm), 
      col = "black", lty = 2, lwd = 1.5)
lines(fig2B2$gamma, apply(fig2B2, 1, wrap_bias_msm), 
      col = "grey", lty = 1, lwd = 1.5)
lines(fig2B2$gamma, apply(fig2B2, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
axis(1, at = c(-5, 0, 5), 
     label = c("", 0, ""), 
     tck = -0.01, cex.axis = 0.75)
axis(1, lwd = 0, lwd.tick = 0, at = c(-2.5, 2.5), 
     label = c("negative", "positive"), cex.axis = 0.75)
mtext("Association between confounder and outcome", 
      side = 1, line = 1.5)
axis(4, at = c(-1.5, 0, 1.5), tck = -0.01, cex.axis = 0.75)
mtext("bias in average treatment effect", side = 4, line = 1.5)
legend(-8, 2.5, xpd = T, 
       legend = c("Confounder and exposure neg ass", 
                  "Confounder and exposure pos ass"),
       lty = c(1, 1),
       col = c("grey", "black"),
       bty = "n",
       horiz = F,
       cex = 1,
       x.intersp = 0.25,
       seg.len = 1.5,
       text.width = c(0, 0.3))
legend(-2, 2, xpd = T, 
       legend = c("Bias < 0"),
       fill = c("lightgrey"),
       border = c("lightgrey"),
       horiz = T,
       bty = "n",
       cex = 0.75,
       x.intersp = 0.25)
legend(-5.5, 1.85, xpd = T, 
       legend = c("MSM-IPW", "conditional model"),
       lty = c(1, 2),
       bty = "n",
       horiz = T,
       cex = 0.75,
       x.intersp = 0.25,
       seg.len = 1.5,
       text.width = c(0, 2.5))
dev.off()

png(paste0(figures_dir,"/fig3.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mgp = c(0, 0.25, 0), mar = c(4, 4, 4, 4))
plot(0, type = "n", xaxt = "n", yaxt = "n", frame.plot = F, ann = F,
     xlim = c(0, 1), ylim = c(-1.5, 1.5))
rect(xleft = c(0, 0), ybottom = c(-1.5, -1.5), 
     xright = c(1, 1), ytop = c(0, 0),
     col = "lightgrey", border = NA)
lines(fig31$lambda, apply(fig31, 1, wrap_bias_msm), 
      col = "black", lty = 1, lwd = 1.5)
lines(fig31$lambda, apply(fig31, 1, wrap_bias_cm), 
      col = "black", lty = 2, lwd = 1.5)
lines(fig32$lambda, apply(fig32, 1, wrap_bias_msm), 
      col = "grey", lty = 1, lwd = 1.5)
lines(fig32$lambda, apply(fig32, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
axis(1, at = c(0, 1), 
     tck = -0.01, cex.axis = 0.75)
mtext("Prevalence of confounding variable", side = 1, line = 1.5)
axis(2, at = c(-1.5, 1.5), tck = -0.01, cex.axis = 0.75)
mtext("bias in average treatment effect", side = 2, line = 1.5)
mtext(paste("Confounder positively associated", "\n",  "with exposure and outcome"), 
      side = 3, line = 2)
legend(0.3, 2, xpd = T, 
       legend = c("Bias < 0"),
       fill = c("lightgrey"),
       border = c("lightgrey"),
       horiz = T,
       bty = "n",
       cex = 0.75,
       x.intersp = 0.25)
legend(-0.1, 1.85, xpd = T, 
       legend = c("MSM-IPW", "conditional model"),
       lty = c(1, 2),
       bty = "n",
       horiz = T,
       cex = 0.75,
       x.intersp = 0.25,
       seg.len = 1.5,
       text.width = c(0, 0.3))
dev.off()

# ------------------------------------------------------------------------------
# 1. Set parameters
# ------------------------------------------------------------------------------
# Figure 2 consist of 4 plots. The x axis of the plot varies the spec from 0 to 
# 1, and the y axis represents bias in the ATE. In the first two plots, pi_1 > 
# pi_0 (scen A, pi_1 = 0.75). In the last two plots, pi_1 < pi_0 (scen B, pi_1 = 
# 0.5). In plot A1 + B1, sens = 0.95, in plot A2 + B2, sens = 0.80. Plot A1-B2 
# consist of four lines: 1) gamma = 2, pi_0 = 0.25; 2) gamma = 2, pi_0 = 0.5; 3) 
# gamma = -2, pi_0 = 0.25; 4) gamma = -2, pi_0 = 0.5.
fig2A1_1 <- data.frame(spec = seq(from = 0, to = 1, length.out = 100),
                       sens = 0.95, 
                       pi_0 = 0.1, 
                       pi_1 = 0.75,
                       gamma = 2,
                       lambda = 0.5)
fig2A1_2 <- fig2A1_1
fig2A1_2$pi_0 <- 0.5
fig2A1_3 <- fig2A1_1
fig2A1_3$gamma <- -2
fig2A1_4 <- fig2A1_2
fig2A1_4$gamma <- -2

fig2A2_1 <- fig2A1_1
fig2A2_1$sens = 0.80
fig2A2_2 <- fig2A1_2
fig2A2_2$sens = 0.80
fig2A2_3 <- fig2A1_3
fig2A2_3$sens = 0.80
fig2A2_4 <- fig2A1_4
fig2A2_4$sens = 0.80

png(paste0(figures_dir,"/fig2A1.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 1, 1, 4))
plot(fig2A1_1$spec, apply(fig2A1_1, 1, wrap_bias_msm), type = 'l', pch = 20, 
     col = "black", lty = 1, xaxt = "n", yaxt = "n", lwd = 1.5, 
     frame.plot = F, ann = F,
     xlim = c(0, 1), ylim = c(-1.5, 1.5))
lines(fig2A1_1$spec, apply(fig2A1_1, 1, wrap_bias_cm), 
      col = "grey", lwd = 1.5)
lines(fig2A1_2$spec, apply(fig2A1_2, 1, wrap_bias_msm), 
      col = "black", lty = 2,  lwd = 1.5)
lines(fig2A1_2$spec, apply(fig2A1_2, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
lines(fig2A1_3$spec, apply(fig2A1_3, 1, wrap_bias_msm), 
      col = "black", lty = 1,  lwd = 1.5)
lines(fig2A1_3$spec, apply(fig2A1_3, 1, wrap_bias_cm), 
      col = "grey", lty = 1, lwd = 1.5)
lines(fig2A1_4$spec, apply(fig2A1_4, 1, wrap_bias_msm), 
      col = "black", lty = 2,  lwd = 1.5)
lines(fig2A1_4$spec, apply(fig2A1_4, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
axis(1, at = c(0, 1))
mtext('specificity', side = 1, line = 1.5)
axis(4, at = c(-1.5, 0, 1.5))
mtext("bias in average treatment effect", side = 4, line = 2.5)
text(fig2A1_1$spec[15], apply(fig2A1_1, 1, wrap_bias_msm)[1] - 0.11, 
     expression(paste(gamma, " = 2, ", pi[0], " = 0.1")), cex = 0.75)
text(fig2A1_2$spec[15], apply(fig2A1_2, 1, wrap_bias_msm)[1] - 0.1, 
     expression(paste(gamma, " = 2, ", pi[0], " = 0.5")), cex = 0.75)
text(fig2A1_3$spec[15], apply(fig2A1_3, 1, wrap_bias_msm)[1] + 0.1, 
     expression(paste(gamma, " = -2, ", pi[0], " = 0.1")), cex = 0.75)
text(fig2A1_4$spec[15], apply(fig2A1_4, 1, wrap_bias_msm)[1] + 0.11, 
     expression(paste(gamma, " = -2, ", pi[0], " = 0.5")), cex = 0.75)
dev.off()

png(paste0(figures_dir,"/fig2A2.png"), 
    width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 4, 1, 1))
plot(fig2A2_1$spec, apply(fig2A2_1, 1, wrap_bias_msm), type = 'l', pch = 20, 
     col = "black", lty = 1, xaxt = "n", yaxt = "n", lwd = 1.5, 
     frame.plot = F, ann = F,
     xlim = c(0, 1), ylim = c(-1.5, 1.5))
lines(fig2A2_1$spec, apply(fig2A2_1, 1, wrap_bias_cm), 
      col = "grey", lwd = 1.5)
lines(fig2A2_2$spec, apply(fig2A2_2, 1, wrap_bias_msm), 
      col = "black", lty = 2,  lwd = 1.5)
lines(fig2A2_2$spec, apply(fig2A2_2, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
lines(fig2A2_3$spec, apply(fig2A2_3, 1, wrap_bias_msm), 
      col = "black", lty = 1,  lwd = 1.5)
lines(fig2A2_3$spec, apply(fig2A2_3, 1, wrap_bias_cm), 
      col = "grey", lty = 1, lwd = 1.5)
lines(fig2A2_4$spec, apply(fig2A2_4, 1, wrap_bias_msm), 
      col = "black", lty = 2,  lwd = 1.5)
lines(fig2A2_4$spec, apply(fig2A2_4, 1, wrap_bias_cm), 
      col = "grey", lty = 2, lwd = 1.5)
axis(1, at = c(0, 1))
mtext('specificity', side = 1, line = 1.5)
axis(2, at = c(-1.5, 0, 1.5))
mtext("bias in average treatment effect", side = 2, line = 2.5)
text(fig2A2_1$spec[85], apply(fig2A2_1, 1, wrap_bias_cm)[100] - 0.1, 
     expression(paste(gamma, " = 2, ", pi[0], " = 0.1")), cex = 0.75)
text(fig2A2_2$spec[85], apply(fig2A2_2, 1, wrap_bias_msm)[100] - 0.1, 
     expression(paste(gamma, " = 2, ", pi[0], " = 0.5")), cex = 0.75)
text(fig2A2_3$spec[85], apply(fig2A2_3, 1, wrap_bias_cm)[100] + 0.1, 
     expression(paste(gamma, " = -2, ", pi[0], " = 0.1")), cex = 0.75)
text(fig2A2_4$spec[85], apply(fig2A2_4, 1, wrap_bias_msm)[100] + 0.1, 
     expression(paste(gamma, " = -2, ", pi[0], " = 0.5")), cex = 0.75)
dev.off()

#Scenarios------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#Plot a. Varying pi_1 and the pi_1/pi_0 ratio and gamma > 0-----------------------------------------------
#---------------------------------------------------------------------------------------------------------
#a1: beta = 1, gamma = 2, lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = varying (from 0 to 1)---
a1 <- matrix(nrow = 100, ncol = 8)
a1[,1] <- 1; a1[,2] <- 2; a1[,3] <- 0.5; a1[,4] <- 0.05; a1[,5] <- 0.9
a1[,6] <- 0.5; a1[,7] <- seq(from = 0, to = 1, length.out = 100); a1[,8] <- a1[,7] / a1[,6] #pi_1/pi_0
#a2: beta = 1, gamma = 2, lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.9, pi_1 = varying (from 0 to 1)---
a2 <- matrix(nrow = 100, ncol = 8)
a2[,1] <- 1; a2[,2] <- 2; a2[,3] <- 0.5; a2[,4] <- 0.05; a2[,5] <- 0.9
a2[,6] <- 0.9; a2[,7] <- seq(from = 0, to = 1, length.out = 100); a2[,8] <- a2[,7] / a2[,6] #pi_1/pi_0
#plot a1 + a2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/pi1pi0_1.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, rep(1, 2), 4))
plot(a1[,8], apply(a1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, 
     frame.plot = F, ann = F,
     xlim = c(0, 2), ylim = c(-1.5, 1.5))
points(a1[1,8], bias_cm(a1[1,]), pch = 21, col = "grey", bg = "white")
points(a1[100,8], bias_cm(a1[100,]), pch = 21, col = "grey", bg = "white")
lines(a1[,8], apply(a1, 1, bias_msm), col = "grey", lwd = 1.5)
points(a1[1,8], bias_msm(a1[1,]), pch = 21, col = "grey", bg = "white")
points(a1[100,8], bias_msm(a1[100,]), pch = 21, col = "grey", bg = "white")
lines(a2[,8], apply(a2, 1, bias_cm), lty = 2, lwd = 1.5)
points(a2[1,8], bias_cm(a2[1,]), pch = 21, bg = "white")
points(a2[100,8], bias_cm(a2[100,]), pch = 21, bg = "white")
lines(a2[,8], apply(a2, 1, bias_msm), lwd = 1.5)
points(a2[1,8], bias_msm(a2[1,]), pch = 21, bg = "white")
points(a2[100,8], bias_msm(a2[100,]), pch = 21, bg = "white")
axis(1, at = c(0, 1, 2))
mtext(expression(pi[1]:pi[0]), side = 1, line = 2.5)
axis(4, at = c(-1.5, 0, 1.5))
mtext("bias in average treatment effect", side = 4, line = 2.5)
legend("topleft", legend = c(expression(paste(pi[0], " = 0.50 ")), expression(paste(pi[0], " = 0.90 "))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n',cex = 0.75)
legend("bottomright", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot b. Varying pi_1 and the pi_1/pi_0 ratio and gamma < 0-----------------------------------------------
#---------------------------------------------------------------------------------------------------------
#b1: beta = 1, gamma = -2, lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = varying (from 0 to 1)--
b1 <- matrix(nrow = 100, ncol = 8)
b1[,1] <- 1; b1[,2] <- -2; b1[,3] <- 0.5; b1[,4] <- 0.05; b1[,5] <- 0.9
b1[,6] <- 0.5; b1[,7] <- seq(from = 0, to = 1, length.out = 100); b1[,8] <- b1[,7] / b1[,6] #pi_1/pi_0
#b2: beta = 1, gamma = -2, lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.9, pi_1 = varying (from 0 to 1)--
b2 <- matrix(nrow = 100, ncol = 8)
b2[,1] <- 1; b2[,2] <- -2; b2[,3] <- 0.5; b2[,4] <- 0.05; b2[,5] <- 0.9; 
b2[,6] <- 0.9; b2[,7] <- seq(from = 0, to = 1, length.out = 100); b2[,8] <- b2[,7] / b2[,6] #pi_1/pi_0
#plot b1 + b2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/pi1pi0_2.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 4, rep(1, 2)))
plot(b1[,8], apply(b1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, 
     frame.plot = F, ann = F, 
     xlim = c(0, 2), ylim = c(-1.5, 1.5))
points(b1[1,8], bias_cm(b1[1,]), pch = 21, col = "grey", bg = "white")
points(b1[100,8], bias_cm(b1[100,]), pch = 21, col = "grey", bg = "white")
lines(b1[,8], apply(b1, 1, bias_msm), col = "grey", lwd = 1.5)
points(b1[1,8], bias_msm(b1[1,]), pch = 21, col = "grey", bg = "white")
points(b1[100,8], bias_msm(b1[100,]), pch = 21, col = "grey", bg = "white")
lines(b2[,8], apply(b2, 1, bias_cm), lty = 2, lwd = 1.5)
points(b2[1,8], bias_cm(b2[1,]), pch = 21, bg = "white")
points(b2[100,8], bias_cm(b2[100,]), pch = 21, bg = "white")
lines(b2[,8], apply(b2, 1, bias_msm), lwd = 1.5)
points(b2[1,8], bias_msm(b2[1,]), pch = 21, bg = "white")
points(b2[100,8], bias_msm(b2[100,]), pch = 21, bg = "white")
axis(1, at = c(0, 1, 2))
mtext(expression(pi[1]:pi[0]), side = 1, line = 2.5)
axis(2, at = c(-1.5, 0, 1.5))
mtext("bias in average treatment effect", side = 2, line = 2.5)
legend("topright", legend = c(expression(paste(pi[0], " = 0.50 ")), expression(paste(pi[0], " = 0.90 "))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("bottomleft", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot c. Varying gamma from -5 to 5 and pi_1/pi_0>1-------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#c1: beta = 1, gamma = varying(from -5 to 5), lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.45, pi_1 = 0.9
c1 <- matrix(nrow = 100, ncol = 7)
c1[,1] <- 1; c1[,2] <- seq(from = -5, to = 5, length.out = 100); c1[,3] <- 0.5; c1[,4] <- 0.05
c1[,5] <- 0.9; c1[,6] <- 0.45; c1[,7] <- 0.9
#c2: beta = 1, gamma = varying(from -5 to 5), lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.20, pi_1 = 0.9
c2 <- matrix(nrow = 100, ncol = 7)
c2[,1] <- 1; c2[,2] <- seq(from = -5, to = 5, length.out = 100); c2[,3] <- 0.5; c2[,4] <- 0.05
c2[,5] <- 0.9; c2[,6] <- 0.2; c2[,7] <- 0.9
#plot c1 + c2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/gamma_1.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, rep(1, 2), 4))
plot(c1[,2], apply(c1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(-5, 5), ylim = c(-2, 2))
axis(1, at = c(-5, 0, 5))
axis(4, at = c(-2, 0, 2))
mtext(expression(gamma), side = 1, line = 2.5)
mtext("bias in average treatment effect", side = 4, line = 2.5)
lines(c1[,2], apply(c1, 1, bias_msm), col = "grey", lwd = 1.5)
lines(c2[,2], apply(c2, 1, bias_cm), lty = 2, lwd = 1.5)
lines(c2[,2], apply(c2, 1, bias_msm), lwd = 1.5)
legend("topleft", legend = c(expression(paste(pi[0], " = 0.45, ", pi[1], " = 0.90")), 
                             expression(paste(pi[0], " = 0.20, ", pi[1], " = 0.90"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("bottomright", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot d. Varying gamma from -5 to 5 and pi_1/pi_0 < 1-----------------------------------------------------
#---------------------------------------------------------------------------------------------------------
#d1: beta = 1, gamma = varying(from -5 to 5), lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.9, pi_1 = 0.45
d1 <- matrix(nrow = 100, ncol = 8)
d1[,1] <- 1; d1[,2] <- seq(from = -5, to = 5, length.out = 100); d1[,3] <- 0.5; d1[,4] <- 0.05
d1[,5] <- 0.9 ;d1[,6] <- 0.9; d1[,7] <- 0.45
#d2: beta = 1, gamma = varying (from -5 to 5), lambda = 0.5, p_0 = 0.05, p_1 = 0.9, pi_0 = 0.9, pi_1 = 0.2
d2 <- matrix(nrow = 100, ncol = 8)
d2[,1] <- 1; d2[,2] <- seq(from = -5, to = 5, length.out = 100); d2[,3] <- 0.5; d2[,4] <- 0.05
d2[,5] <- 0.9 ;d2[,6] <- 0.9; d2[,7] <- 0.2 
#plot d1 + d2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/gamma_2.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 4, rep(1, 2)))
plot(d1[,2], apply(d1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(-5, 5), ylim = c(-2, 2))
axis(1, at = c(-5, 0, 5))
axis(2, at = c(-2, 0, 2))
mtext(expression(gamma), side = 1, line = 2.5)
mtext("bias in average treatment effect", side = 2, line = 2.5)
lines(d1[,2], apply(d1, 1, bias_msm), col = "grey", lwd = 1.5)
lines(d2[,2], apply(d2, 1, bias_cm), lty = 2, lwd = 1.5)
lines(d2[,2], apply(d2, 1, bias_msm), lwd = 1.5)
legend("topright", legend = c(expression(paste(pi[0], " = 0.90, ", pi[1], " = 0.45")), 
                              expression(paste(pi[0], " = 0.90, ", pi[1], " = 0.20"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("bottomleft", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot e. Varying lambda from -5 to 5, pi_1/pi_0 > 1 and gamma > 0-----------------------------------------
#---------------------------------------------------------------------------------------------------------
#e1: beta = 1, gamma = 2, lambda = (0,1), p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = 0.75, pi_1/pi_0 = 1.5
e1 <- matrix(nrow = 100, ncol = 7)
e1[,1] <- 1; e1[,2] <- 2; e1[,3] <- seq(0, 1, length.out = 100); e1[,4] <- 0.05; e1[,5] <- 0.9 
e1[,6] <- 0.5; e1[,7] <- 0.75 
#e2: beta = 1, gamma = 2, lambda = (0,1), p_0 = 0.05, p_1 = 0.9, pi_0 = 0.25, pi_1 = 0.75, pi_1/pi_0 = 3
e2 <- matrix(nrow = 100, ncol = 7)
e2[,1] <- 1; e2[,2] <- 2; e2[,3] <- seq(0, 1, length.out = 100); e2[,4] <- 0.05; e2[,5] <- 0.9
e2[,6] <- 0.25; e2[,7] <- 0.75  
#plot e1 + e2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/lambda_1.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 1, 1, 4))
plot(e1[,3], apply(e1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(0, 1), ylim = c(0, 0.5))
lines(e1[,3], apply(e1, 1, bias_msm), col = "grey", lwd = 1.5)
lines(e2[,3], apply(e2, 1, bias_cm), lty = 2, lwd = 1.5)
lines(e2[,3], apply(e2, 1, bias_msm), lwd = 1.5)
#points(e1[1,3], bias_cm(e1[1,]), pch = 21, bg = "white")
#points(e1[100,3], bias_cm(e1[100,]), pch = 21, bg = "white")
axis(1, at = c(0, 1))
mtext(expression(lambda), side = 1, line = 2.5)
axis(4, at = c(0, 0.5))
mtext("bias in average treatment effect", side = 4, line = 2.5)
legend("topleft", legend = c(expression(paste(pi[0], " = 0.50, ", pi[1], " = 0.75")), 
                             expression(paste(pi[0], " = 0.25, ", pi[1], " = 0.75"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("topright", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot f. Varying lambda from -5 to 5, pi_1/pi_0 > 1 and gamma < 0-----------------------------------------
#---------------------------------------------------------------------------------------------------------
#f1: beta = 1, gamma = -2, lambda = (0,1), p_0 = 0.05, p_1 = 0.9, pi_0 = 0.5, pi_1 = 0.75, pi_1/pi_0 = 1.5
f1 <- matrix(nrow = 100, ncol = 7)
f1[,1] <- 1; f1[,2] <- -2; f1[,3] <- seq(0, 1, length.out = 100); f1[,4] <- 0.05; f1[,5] <- 0.9
f1[,6] <- 0.5; f1[,7] <- 0.75 
#f2: beta = 1, gamma = -2, lambda = (0,1), p_0 = 0.05, p_1 = 0.9, pi_0 = 0.25, pi_1 = 0.75, pi_1/pi_0 = 3
f2 <- matrix(nrow = 100, ncol = 7)
f2[,1] <- 1; f2[,2] <- -2; f2[,3] <- seq(0, 1, length.out = 100); f2[,4] <- 0.05; f2[,5] <- 0.9
f2[,6] <- 0.25; f2[,7] <- 0.75
#plot f1 + f2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/lambda_2.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 4, 1, 1))
plot(f1[,3], apply(f1, 1, bias_cm), type = 'l', pch = 20, 
     col = "grey", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(0, 1), ylim = c(-0.5, 0))
lines(f1[,3], apply(f1, 1, bias_msm), col = "grey", lwd = 1.5)
lines(f2[,3], apply(f2, 1, bias_cm), lty = 2, lwd = 1.5)
lines(f2[,3], apply(f2, 1, bias_msm), lwd = 1.5)
#points(f1[1,3], bias_cm(m7[1,]), pch = 21, bg = "white")
#points(f1[100,3], bias_cm(m7[100,]), pch = 21, bg = "white")
axis(1, at = c(0, 1))
mtext(expression(lambda), side = 1, line = 2.5)
axis(2, at = c(-0.5, 0))
mtext("bias in average treatment effect", side = 2, line = 2.5)
legend("bottomright", legend = c(expression(paste(pi[0], " = 0.50, ", pi[1], " = 0.75")), 
                             expression(paste(pi[0], " = 0.25, ", pi[1], " = 0.75"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("bottomleft", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot g. Varying p_0 from 0 to 1, pi_1/pi_0 > 1 and gamma > 0---------------------------------------------
#---------------------------------------------------------------------------------------------------------
#g1: beta = 1, gamma = 2, lambda = 0.5, p_0 = c(0,1), p_1 = 0.95, pi_0 = 0.5, pi_1 = 0.75, pi_1/pi_0 = 1.5
g1 <- matrix(nrow = 100, ncol = 7)
g1[,1] <- 1; g1[,2] <- 2; g1[,3] <- 0.5; g1[,4] <- seq(0, 1, length.out = 100); g1[,5] <- 0.95
g1[,6] <- 0.5; g1[,7] <- 0.75 #;g1[,8] <- g1[,4]/g1[,5]
#g2: beta = 1, gamma = 2, lambda = 0.8, p_0 = c(0,1), p_1 = 0.80, pi_0 = 0.25, pi_1 = 0.75, pi_1/pi_0 = 3
g2 <- matrix(nrow = 100, ncol = 7)
g2[,1] <- 1; g2[,2] <- 2; g2[,3] <- 0.5; g2[,4] <- seq(0, 1, length.out = 100); g2[,5] <- 0.8; 
g2[,6] <- 0.25; g2[,7] <- 0.75 #;g2[,8] <- g2[,4]/g2[,5]
#plot g1 + g2---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------
png(paste0(figures_dir,"/p_0.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 1, 1, 4))
plot(g1[,4], apply(g1, 1, bias_cm), type = 'l', pch = 20, 
     col = "black", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(0, 1), ylim = c(0, 1))
lines(g1[,4], apply(g1, 1, bias_msm), col = "black", lwd = 1.5)
lines(g2[,4], apply(g2, 1, bias_cm), col = "grey", lty = 2, lwd = 1.5)
lines(g2[,4], apply(g2, 1, bias_msm), col = "grey", lwd = 1.5)
axis(1, at = c(0, 1))
mtext("1 - specificity", side = 1, line = 2.5)
axis(4, at = c(0, 1))
mtext("bias in average treatment effect", side = 4, line = 2.5)
legend("bottomright", legend = c(expression(paste("sens = 0.95, ", pi[0], " = 0.50, ", pi[1], " = 0.75")), 
                                 expression(paste("sens = 0.80, ", pi[0], " = 0.25, ", pi[1], " = 0.75"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("topleft", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()

#Plot h. Varying p_1 from 0 to 1, pi_1/pi_0 > 1 and gamma > 0---------------------------------------------
#---------------------------------------------------------------------------------------------------------
#h1: beta = 1, gamma = 2, lambda = 0.5, p_0 = 0.05, p_1 = c(0,1), pi_0 = 0.5, pi_1 = 0.75, pi_1/pi_0 = 1.5
h1 <- matrix(nrow = 100, ncol = 7)
h1[,1] <- 1; h1[,2] <- 2; h1[,3] <- 0.5; h1[,4] <- 0.05; 
h1[,5] <- seq(0, 1, length.out = 100); h1[,6] <- 0.5; h1[,7] <- 0.75 #;h1[,8] <- h1[,4]/h1[,5]
#h2: beta = 1, gamma = 2, lambda = 0.5, p_0 = 0.2, p_1 = c(0,1), pi_0 = 0.25, pi_1 = 0.75, pi_1/pi_0 = 3
h2 <- matrix(nrow = 100, ncol = 7)
h2[,1] <- 1; h2[,2] <- 2; h2[,3] <- 0.5; h2[,4] <- 0.2; h2[,5] <- seq(0, 1, length.out = 100)
h2[,6] <- 0.25; h2[,7] <- 0.75 #;h2[,8] <- h2[,4]/h2[,5]
png(paste0(figures_dir,"/p_1.png"), width = 4, height = 4, units = 'in', res = 100)
par(mar = c(3.5, 4, 1, 1))
plot(h1[,5], apply(h1, 1, bias_cm), type = 'l', pch = 20, 
     col = "black", lty = 2, xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, 
     xlim = c(0, 1), ylim = c(0, 1))
lines(h1[,5], apply(h1, 1, bias_msm), col = "black", lwd = 1.5)
lines(h2[,5], apply(h2, 1, bias_cm), col = "grey", lty = 2, lwd = 1.5)
lines(h2[,5], apply(h2, 1, bias_msm), col = "grey", lwd = 1.5)
axis(1, at = c(0, 1))
mtext("sensitivity", side = 1, line = 2.5)
axis(2, at = c(0, 1))
mtext("bias in average treatment effect", side = 2, line = 2.5)
legend("bottomleft", legend = c(expression(paste("1-spec = 0.05, ", pi[0], " = 0.50, ", pi[1], " = 0.75")), 
                                 expression(paste("1-spec = 0.20, ", pi[0], " = 0.25, ", pi[1], " = 0.75"))),
       seg.len = 1, col = c("grey", "black"), lty = 1, lwd = 1.5, bty = 'n', cex = 0.75)
legend("topright", legend = c("cm ", "msm "), bty = 'n',
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, cex = 0.75)
dev.off()
