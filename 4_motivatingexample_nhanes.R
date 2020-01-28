#----------------------------------------------------------------------------------------------
# Motivating example ME_MSM using NHANES data
#----------------------------------------------------------------------------------------------
#require(RNHANES)
require(tidyverse)
require(ipw)
require(survey)
nhanes_dir <- "./nhanes"
figures_dir <- "./figures"

# Load NHANES Data-----------------------------------------------------------------------------
# 2013-2014------------------------------------------------------------------------------------
# Data <- nhanes_load_data(c("BMX_H","BPX_H","WHQ_H"), "2013-2014", demographics = F) 
# #DEMO: DEMOGRAPHICS, BPQ: BLOOD PRES QUESTIONNAIRE,
# NROW(Data$BMX_H) #9813
# NROW(Data$BPX_H) #9813
# NROW(Data$WHQ_H) #6464
# setdiff(Data$BMX_H$SEQN, Data$BPX_H$SEQN) #difference is 0
# setdiff(Data$WHQ_H$SEQN, Data$BMX_H$SEQN) #everyone who is in WHQ_H but not in BMX_H
# NROW(setdiff(Data$WHQ_H$SEQN, Data$BMX_H$SEQN)) #198 objects are in WHQ_H but not in BMX_H
# #these are people that are >16 years, but were not part of the examined sample
# #Hence, there are 6464 - 198 = 6266 objects in the combined data 
# DataF <- Reduce(function(x, y) merge(x, y, by="SEQN"), Data)
# DataF <- DataF[, !duplicated(colnames(DataF))]
# NROW(DataF) #6266 objects
# saveRDS(DataF, file = paste0(nhanes_dir, "/NHANES.rds"))
DataF <- readRDS(file = paste0(nhanes_dir, "/NHANES.rds"))
# 2011-2012------------------------------------------------------------------------------------
# Data2 <- nhanes_load_data(c("BMX","BPX","WHQ"), "2011-2012", demographics = F)
# NROW(Data2$BMX) #9338
# NROW(Data2$BPX) #9338
# NROW(Data2$WHQ) #6175
# setdiff(Data2$BMX$SEQN, Data2$BPX$SEQN) #difference is 0
# setdiff(Data2$WHQ$SEQN, Data2$BMX$SEQN) #everyone who is in WHQ_H but not in BMX_H
# NROW(setdiff(Data2$WHQ$SEQN, Data2$BMX$SEQN)) #256 objects are in WHQ but not in BMX
# #these are people that are >16 years, but were not part of the examined sample
# #Hence, there are 6175 - 256 = 5919 objects in the combined data 
# DataF2 <- Reduce(function(x, y) merge(x, y,by="SEQN"), Data2)
# DataF2 <- DataF2[, !duplicated(colnames(DataF2))]
# NROW(DataF2) #5919 objects
# saveRDS(DataF2, paste0(nhanes_dir, "/NHANES2.rds"))
DataF2 <- readRDS(file = paste0(nhanes_dir, "/NHANES2.rds"))
#9813 + 9338 = 19151 objects are physically examined

# Recode and select variables NHANES data------------------------------------------------------
vars <- c("SEQN"="id", "BMXWT"="weight","BMXHT"="height","BMXBMI"="bmi","BPXSY1"="sbp1",
          "BPXDI1"="dbp1","BPXSY2"="sbp2","BPXDI2"="dbp2","BPXSY3"="sbp3","BPXDI3"="dbp3",
          "WHD020"="weight_self", "WHD010"="height_self")
DataS1 <- DataF %>% select(names(vars)) %>% plyr::rename(vars)
DataS2 <- DataF2 %>% select(names(vars)) %>% plyr::rename(vars)
# NROW(DataS1) #6266 participants
# NROW(DataS2) #5919 participants

# Create one data file from DataS1 and DataS2 (DataS) -----------------------------------------
DataS <- rbind(DataS1, DataS2)
# NROW(DataS) #12185
# Data
# DataS$weight  # weight in KG by trained health technicians
# DataS$height  # height in cm by trained health technicians
# DataS$bmi     # bmi by trained health technicians
# DataS$sbp1    # blood pressure systolisch in mm Hg first reading trained BP examiner 
                # (after 5 min resting)
DataS$dbp1[DataS$dbp1==0] <- NA # blood pressure diastolisch in mm Hg first reading trained BP 
                                # examiner (after 5 min resting)
# DataS$sbp2    # blood pressure systolic second reading
DataS$dbp2[DataS$dbp2==0] <- NA   # blood pressure diastolic second reading
# DataS$sbp3    # blood pressure systolic third reading
DataS$dbp3[DataS$dbp3==0] <- NA    # blood pressure diastolic third reading
DataS$weight_self[DataS$weight_self==9999|DataS$weight_self==7777] <- NA
DataS$weight_selfKG <- DataS$weight_self*0.45359237 # self reported weight in KG
DataS$height_self[DataS$height_self==9999] <- NA
DataS$height_selfCM <- DataS$height_self*2.54 #self reported height in cm
DataS$bmi_self <- round(DataS$weight_selfKG / (DataS$height_selfCM/100)^2,1)
#DataS$bmi_2 <- DataS$weight / (DataS$height/100)^2
DataS <- cbind(DataS, "sbpm" = rowMeans(with(DataS, cbind(sbp1, sbp2, sbp3)), na.rm=TRUE)) 
#NaN for the ones that have not sbp1 and sbp2 and spb3
# NROW(DataS[is.nan(DataS$sbpm),]) #494 objects do not have any sbp measure
# NROW(DataS[!is.nan(DataS$sbpm),]) #11691 objects (12185 - 494)

# Load RXQ_DRUG to get the drug id's-----------------------------------------------------------
# see https://wwwn.cdc.gov/Nchs/Nhanes/1999-2000/RXQ_DRUG.htm----------------------------------
# med <- nhanes_load_data(file_name = "RXQ_DRUG", year = "1999-2000")
# saveRDS(med, file = paste0(nhanes_dir, "/med.rds"))
med <- readRDS(file = paste0(nhanes_dir, "/med.rds"))
# Load medication data
# DRUG ID can be linked to DRUG ID in RXQ_RX data (actual use in participants)
# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/RXQ_RX_H.htm
# Data_med1 <- nhanes_load_data(file_name = "RXQ_RX", year = "2013-2014")
# Data_med2 <- nhanes_load_data(file_name = "RXQ_RX", year = "2011-2012")
# saveRDS(Data_med1, file = "Data_med1.rds")
# saveRDS(Data_med2, file="Data_med2.rds")
Data_med1 <- readRDS(file = paste0(nhanes_dir, "/Data_med1.rds"))
Data_med2 <- readRDS(file = paste0(nhanes_dir, "/Data_med2.rds"))
# NROW(intersect(DataF$SEQN, Data_med1$SEQN))#all 6266 ind of DataF are in the RXQ_RX file
# NROW(intersect(DataF2$SEQN, Data_med2$SEQN))#all 5919 ind of DataF2 are in the RXQ_RX file
# NROW(intersect(Data$BMX_H$SEQN, Data_med1$SEQN)) #all 9813 ind (total pop)
# NROW(intersect(Data2$BMX$SEQN, Data_med2$SEQN)) #all 9938 ind (total pop)

# Individuals that take diuretics (id_diu) or beta blockers (id_bb) or both (id_bb) 2013-2014--
id_diu1 <- unique(Data_med1[Data_med1$RXDDRGID %in% 
                              med[med$RXDDCN1B == "DIURETICS",]$RXDDRGID,"SEQN"])
id_bb1  <- unique(Data_med1[Data_med1$RXDDRGID %in% 
                              med[med$RXDDCN1C == "BETA BLOCKERS, CARDIOSELECTIVE" | 
                                  med$RXDDCN1C == "BETA BLOCKERS, NON-CARDIOSELECTIVE",
                                  ]$RXDDRGID, "SEQN"])
id_both1 <- intersect(id_diu1, id_bb1)
# NROW(id_diu1) # 496 participants take diuretics
# NROW(id_bb1) # 655 participants take beta blocker
# NROW(id_both1) # 225 particpants take both
id_diu1 <- id_diu1[id_diu1 %in% id_both1 == F]
id_bb1 <- id_bb1[id_bb1 %in% id_both1 == F]
# NROW(id_diu1) #271 take only a diuretic
# NROW(id_bb1) #430 take only beta blockers
# intersect(id_diu1, id_bb1) #intersect is 0
# Individuals that take diuretics (id_diu) or beta blockers (id_bb) or both (id_bb) 2011-2012--
id_diu2 <- unique(Data_med2[Data_med2$RXDDRGID %in% 
                              med[med$RXDDCN1B == "DIURETICS",]$RXDDRGID, "SEQN"])
id_bb2  <- unique(Data_med2[Data_med2$RXDDRGID %in% 
                              med[med$RXDDCN1C == "BETA BLOCKERS, CARDIOSELECTIVE" | 
                                  med$RXDDCN1C == "BETA BLOCKERS, NON-CARDIOSELECTIVE",
                                  ]$RXDDRGID, "SEQN"])
id_both2 <- intersect(id_diu2, id_bb2)
# NROW(id_diu2) # 511 participants take diuretics
# NROW(id_bb2) # 592 participants take beta blocker
# NROW(id_both2) # 198 particpants take both
id_diu2 <- id_diu2[id_diu2 %in% id_both2 == F]
id_bb2 <- id_bb2[id_bb2 %in% id_both2 == F]
# NROW(id_diu2) #313 take only a diuretic
# NROW(id_bb2) #394 take only beta blockers
# intersect(id_diu2, id_bb2) #intersect is 0
# Select individuals 2011-2014 that take diuretics or beta blockers (DataS_s)------------------ 
id_diu <- c(id_diu1, id_diu2)
id_bb <- c(id_bb1, id_bb2)
# NROW(id_diu) #584 diuretics in 2011-2014 cycle
# NROW(id_bb) #824 diuretics in 2011-2014 cycle
# NROW(id_diu) / (NROW(Data$BMX_H) + NROW(Data2$BMX)) #0.03 in total pop 2011-2014
# NROW(id_bb) / (NROW(Data$BMX_H) + NROW(Data2$BMX)) #0.04 in total pop 2011-2014
# NROW(DataS[DataS$id %in% c(id_diu),]) #551 objects are in the self-rep sample and take diuretics
# NROW(DataS[DataS$id %in% c(id_bb),]) #787 objects are in the self-rep sample and take diuretics
# select individuals that are in self-rep sample
id_diu <- DataS[DataS$id %in% id_diu, "id"]
id_bb <- DataS[DataS$id %in% id_bb, "id"]

# Select subsample of DataS with individuals that are in the self-rep sample and took meds-----
# (DataS_s)
DataS_s <- DataS[DataS$id %in% c(id_diu, id_bb),]
# NROW(DataS_s) #1338 objects
DataS_s$bpdrug <- ifelse(DataS_s$id %in% id_diu, 1, 0) #diuretic = 1, beta blocker = 0

# Select complete data (complete data analysis) (DataS_s_c)------------------------------------
DataS_s_c <- DataS_s[!is.nan(DataS_s$sbpm) & !is.na(DataS_s$bmi) & !is.na(DataS_s$bmi_self),]
# NROW(DataS_s_c) #1227
mean(DataS_s_c[DataS_s_c$bpdrug == 1, "bmi_self"], na.rm = T) #31.5
mean(DataS_s_c[DataS_s_c$bpdrug == 0, "bmi_self"], na.rm = T) #29.2
mean(DataS_s_c[DataS_s_c$bpdrug == 1, "bmi"], na.rm = T) #32.5
mean(DataS_s_c[DataS_s_c$bpdrug == 0, "bmi"], na.rm = T) #30.2

# Create bmi categories (bmi_c and bmi_self_c) ------------------------------------------------
DataS_s_c$bmi_c <- ifelse(DataS_s_c$bmi < 25, 0, 1)
DataS_s_c$bmi_self_c <- ifelse(DataS_s_c$bmi_self < 25, 0, 1)

# Spec and sens in data -----------------------------------------------------------------------
table(DataS_s_c$bmi_c, DataS_s_c$bmi_self_c, dnn = c("bmi", "bmi_self")) #misclassification
#P(bmi_self = overweight/obese | bmi = overweight/obese) = 925 / (55+925) = 0.94 (spec)
#P(bmi_self = overweight/obese | bmi = underweight/normal) = 20 / (227 + 20) = 0.08 (1-sens)

# Analyse data using using self-reported bmi categories ---------------------------------------
fit_unadj <- lm(sbpm ~ bpdrug, data = DataS_s_c) #-4.03
summary(fit_unadj)
#ci 
lower_unadj <- fit_unadj$coefficients[2] - qnorm(0.975) * summary(fit_unadj)$coef[2,2]
upper_unadj <- fit_unadj$coefficients[2] + qnorm(0.975) * summary(fit_unadj)$coef[2,2]
# fit_adj <- lm(sbpm ~ bpdrug + bmi_c, data = DataS_s_c) #-3.60
# summary(fit_adj)
# conditional model----------------------------------------------------------------------------
fit_adj_me <- lm(sbpm ~ bpdrug + bmi_self_c, data = DataS_s_c) #-3.48
summary(fit_adj_me)
lower_adj_me <- fit_adj_me$coefficients[2] - qnorm(0.975) * summary(fit_adj_me)$coef[2,2]
upper_adj_me <- fit_adj_me$coefficients[2] + qnorm(0.975) * summary(fit_adj_me)$coef[2,2]
# msm------------------------------------------------------------------------------------------
temp <- ipwpoint(exposure = bpdrug,
                 family = "binomial",
                 link = "logit",
                 denominator = ~ bmi_c,
                 numerator = ~ 1,
                 data = DataS_s_c)
fit_msm <- svyglm(sbpm ~ bpdrug,
                  design = svydesign(~ 1,
                                     weights = ~ temp$ipw.weights,
                                     data = DataS_s_c))
summary(fit_msm) #-3.59
lower_msm <- fit_msm$coefficients[2] + qnorm(0.975) * summary(fit_msm)$coef[2,2]
upper_msm <- fit_msm$coefficients[2] - qnorm(0.975) * summary(fit_msm)$coef[2,2]
temp <- ipwpoint(exposure = bpdrug, 
                 family = "binomial", 
                 link = "logit", 
                 denominator = ~ bmi_self_c,
                 numerator = ~ 1,
                 data = DataS_s_c)
fit_msm_me <- svyglm(sbpm ~ bpdrug, 
              design = svydesign(~ 1, 
                                 weights = ~ temp$ipw.weights, 
                                 data = DataS_s_c))
summary(fit_msm_me) #-3.52
lower_msm_me <- fit_msm_me$coefficients[2] - qnorm(0.975) * summary(fit_msm_me)$coef[2,2]
upper_msm_me <- fit_msm_me$coefficients[2] + qnorm(0.975) * summary(fit_msm_me)$coef[2,2]

# Sensitivity analysis------------------------------------------------------------------------
# Parameters in data--------------------------------------------------------------------------
# in a sens analyse, you'll know:
ell_data <- mean(DataS_s_c$bmi_self_c) #P(L*=1) = 0.77 (l)
pi_star_1_data <- mean(DataS_s_c$bpdrug[DataS_s_c$bmi_self_c == 1]) #P(A=1|L*=1) = 0.44 (pi_star_1)
pi_star_0_data <- mean(DataS_s_c$bpdrug[DataS_s_c$bmi_self_c == 0]) #P(A=1|L*=0) = 0.32 (pi_star_0)
omega_data <- mean(DataS_s_c$bpdrug) #P(A)=0.41
# you'll assume:
# We know what p_1, p_0 and lambda = P(L=1) is in our data
# mean(DataS_s_c$bmi_self_c[DataS_s_c$bmi_c == 0]) #p_0=0.08
# mean(DataS_s_c$bmi_self_c[DataS_s_c$bmi_c == 1]) #p_1=0.94
# mean(DataS_s_c$bmi_c) #P(L=1)=0.80
# P(L) is,
p_0 <- 0.08
p_1 <- 0.94
#from this, we can calculate lambda and pi_0 and pi_1:
calc_lambda <- function(ell, p_0, p_1){
  (ell - p_0) / (p_1 - p_0)
} 
calc_lambda(ell_data, p_0, p_1) 
# mean(DataS_s_c$bpdrug[DataS_s_c$bmi_c == 1]) #the 'real' pi_1=0.44
# mean(DataS_s_c$bpdrug[DataS_s_c$bmi_c == 0]) #the 'real' pi_0=0.32
#calc pi_1
calc_pi_1 <- function(pi_star_0, pi_star_1, ell, p_0, p_1){
  term1 <- (pi_star_1 * ell - pi_star_0 * (1 - ell) * (p_0 / (1 - p_0))) /
    (p_1 * calc_lambda(ell, p_0, p_1))
  term2 <- ((1 - p_0) * p_1) / ((1 - p_0) * p_1 - (1 - p_1) * p_0)
  return(term1 * term2)
}
calc_pi_1(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1) #0.45
#calc pi_0
calc_pi_0 <- function(pi_star_0, pi_star_1, ell, p_0, p_1){
  (pi_star_0 * (1 - ell) - calc_pi_1(pi_star_0, pi_star_1, ell, p_0, p_1) * 
     (1 - p_1) * calc_lambda(ell, p_0, p_1)) / 
    ((1 - p_0) * (1 - calc_lambda(ell, p_0, p_1)))
}
calc_pi_0(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1) #0.28
#calc phi
# mean(DataS_s_c$bmi_c[DataS_s_c$bpdrug == 0 & DataS_s_c$bmi_self_c == 0]) #0.20
# mean(DataS_s_c$bmi_c[DataS_s_c$bpdrug == 1 & DataS_s_c$bmi_self_c == 0]) #0.17
# mean(DataS_s_c$bmi_c[DataS_s_c$bpdrug == 0 & DataS_s_c$bmi_self_c == 1]) #0.97
# mean(DataS_s_c$bmi_c[DataS_s_c$bpdrug == 1 & DataS_s_c$bmi_self_c == 1]) #0.99
phi <- function(pi_star_0, pi_star_1, ell, p_0, p_1, a, l_star){
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
phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 0, 0) #0.17
phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 1, 0) #0.30
phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 0, 1) #0.97
phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 1, 1) #0.99
#calc gamma
fit_cm_int <- lm(sbpm ~ bpdrug * bmi_self_c, data = DataS_s_c)
gamma <- fit_cm_int$coef[3]/(phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 0, 1) - 
                               phi(pi_star_0_data, pi_star_1_data, ell_data, p_0, p_1, 0, 0))
beta_star <- fit_msm_me$coef[2]
DataS_s_c$q <- with(DataS_s_c, bpdrug * bmi_self_c)
beta_star_cm <- fit_adj_me$coef[2]
u1 <- lm(q ~ bpdrug + bmi_self_c, data = DataS_s_c)$coef[2]

# sensitivity analysis msm 
sens_msm <- function(df){
  pi_star_0 <- df["pi_star_0"]
  pi_star_1 <- df["pi_star_1"]
  ell <- df["ell"]
  p_0 <- df["p_0"]
  p_1 <- df["p_1"]
  gamma <- df["gamma"]
  beta_star <- df["beta_star"]
  bias_msm <- gamma * (1 - ell) * 
    (phi(pi_star_0, pi_star_1, ell, p_0, p_1, 1, 0) - 
       phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 0)) +
    gamma * ell * (phi(pi_star_0, pi_star_1, ell, p_0, p_1, 1, 1) - 
                     phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 1))
  beta_msm <- beta_star - bias_msm
  return(c(bias = unname(bias_msm), beta = unname(beta_msm)))}
df <- data.frame("pi_star_0" = pi_star_0_data, 
                      "pi_star_1" = pi_star_1_data, 
                      "ell" = ell_data, 
                      "p_0" = p_0, 
                      "p_1" = p_1, 
                      "gamma" = gamma, 
                      "beta_star" = beta_star)
sens_msm(df)

# sensitivity analysis conditional model
sens_cm <- function(df){
  pi_star_0 <- df["pi_star_0"]
  pi_star_1 <- df["pi_star_1"]
  ell <- df["ell"]
  p_0 <- df["p_0"]
  p_1 <- df["p_1"]
  gamma <- df["gamma"]
  beta_star <- df["beta_star"]
  u1 <- df["u1"]
  bias_cm <- gamma * (phi(pi_star_0, pi_star_1, ell, p_0, p_1, 1, 0) - 
                        phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 0)) * 
    (1 - u1) + gamma * (phi(pi_star_0, pi_star_1, ell, p_0, p_1, 1, 1) - 
                          phi(pi_star_0, pi_star_1, ell, p_0, p_1, 0, 1)) * u1
  beta_cm <- beta_star_cm - bias_cm
  return(c(bias = unname(bias_cm), beta = unname(beta_cm)))}
df$u1 <- u1
sens_cm(df)

# sample p_0 and p_1 from a beta distribution
p_0_sample <- 1 - rbeta(n = 50, 15, 2)
p_1_sample <- rbeta(n = 50, 15, 2)
parGridBeta <- expand.grid("p_0" = p_0_sample, "p_1" = p_1_sample)
parGridBeta$pi_star_0 <- pi_star_0_data
parGridBeta$pi_star_1 <- pi_star_1_data
parGridBeta$ell <- ell_data
parGridBeta$gamma <- gamma
parGridBeta$beta_star <- beta_star
parGridBeta$u1 <- u1
sensResCmBeta <- apply(parGridBeta, 1, sens_cm)
sensResMsmBeta <- apply(parGridBeta, 1, sens_msm)
parGridBeta$bias_cm <- sensResCmBeta[1,]
parGridBeta$beta_cm <- sensResCmBeta[2,]
parGridBeta$bias_msm <- sensResMsmBeta[1,]
parGridBeta$beta_msm <- sensResMsmBeta[2,]

# sample p_0 and p_1 from a uniform distribution
p_0_sample <- seq(from = 0.02, to = 0.10, length.out = 50)
p_1_sample <- seq(from = 0.90, to = 0.98, length.out = 50)
#p_0_sample <- runif(n = 50, min = 0.02, max = 0.10)
#p_1_sample <- runif(n = 50, min = 0.90, max = 0.98)
parGridUnif <- expand.grid("p_0" = p_0_sample, "p_1" = p_1_sample)
parGridUnif$pi_star_0 <- pi_star_0_data
parGridUnif$pi_star_1 <- pi_star_1_data
parGridUnif$ell <- ell_data
parGridUnif$gamma <- gamma
parGridUnif$beta_star <- beta_star
parGridUnif$u1 <- u1
sensResCmUnif <- apply(parGridUnif, 1, sens_cm)
sensResMsmUnif <- apply(parGridUnif, 1, sens_msm)
parGridUnif$bias_cm <- sensResCmUnif[1,]
parGridUnif$beta_cm <- sensResCmUnif[2,]
parGridUnif$bias_msm <- sensResMsmUnif[1,]
parGridUnif$beta_msm <- sensResMsmUnif[2,]

# summary statistics of sensitivity analysis
mean(parGridUnif$bias_msm)
median(parGridUnif$bias_msm)
var(parGridUnif$bias_msm)
quantile(parGridUnif$bias_msm)
mean(parGridBeta$bias_msm)
median(parGridBeta$bias_msm)
var(parGridBeta$bias_msm)
quantile(parGridBeta$bias_msm)

# create plot 
png(paste0(figures_dir, "/sensanalyse2.png"), width = 6, height = 4, 
    units = 'in', res = 300)
par(mar = c(4, 3, rep(1, 2)))
plot(density(parGridUnif$bias_msm), xaxt = "n", yaxt = "n", lwd = 1.5, frame.plot = F,
     ann = F, zero.line = F,
     xlim = c(-2.5, 0.5), ylim = c(0, 4))
lines(density(parGridUnif$bias_cm), lty = 2)
lines(density(parGridBeta$bias_msm), lty = 1, col = "grey")
lines(density(parGridBeta$bias_cm), lty = 2, col = "grey")
abline(v = 0, lty = 3)
axis(1, at = c(-2.5, 0.0, 0.5))
axis(2, labels = F, at = c(0, 4))
mtext("bias in average treatment effect", side = 1, line = 2.5)
mtext("density", side = 2, line = 1.5)
legend("topright", legend = c("cm", "msm"),
       seg.len = 1, col = "black", lty = c(2, 1), lwd = 1.5, bty = 'n')
legend("topleft", legend = c("uniform", "beta"),
       seg.len = 1, col = c("black", "grey"), lty = 1, lwd = 1.5, bty = 'n')
dev.off()
