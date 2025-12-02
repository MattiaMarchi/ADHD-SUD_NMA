# R code for replicating the dataset and the analyses in the paper #
# Pharmacological treatments for ADHD and comorbid substance use disorder: a systematic review and network meta-analysis #
# R code by Mattia Marchi (mattiamarchimd@gmail.com) 
# December 2, 2025

###--------------------------------------------------------------------------------------###
###-------------------------Pairwise Meta-Analysis ADHD-SUD------------------------------###
###--------------------------------------------------------------------------------------###
#Load required packages
library(meta)
library(tidyverse)
library(metafor)
library(netmeta)
#----------------------------1. MA ADHD treatment response--------------------------------#
adhd_r <- structure(list(ID = c("Levin et al. 2006", "Levin et al. 2007", "Levin et al. 2015", "Levin et al. 2024", "Riggs et al. 2004", "Riggs et al. 2011", "Schubiner et al. 2002", "Thurstone et al. 2010"),
                         t1 = c("sustained-release MPH", "sustained-release MPH", "extended-release mixed amphetamine salts (60mg)", "extended-release mixed amphetamine salts (80mg)", "Pemoline", "Osmotic-Release Methylphenidate", "methylphenidate", "Atomoxetine"),
                         t2 = c("PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO"),
                         contrast = c("sustained-release MPH- PBO", "sustained-release MPH- PBO", "extended-release mixed amphetamine salts (60mg) - PBO", "extended-release mixed amphetamine salts (80mg) - PBO", "Pemoline-PBO", "Osmotic-Release Methylphenidate-PBO", "methylphenidate-PBO", "Atomoxetine-PBO"),
                         Comparator = c("PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO", "PBO"), n1 = c(65L, 53L, 83L, 13L, 35L, 151L, 24L, 35L), n2 = c(33L, 53L, 43L, 15L, 34L, 152L, 24L, 35L),
                         n_adhdR_1 = c(16L, 18L, 31L, 10L, 11L, 39L, 18L, 17L), n_adhdR_2 = c(13L, 16L, 5L, 10L, 4L, 35L, 5L, 20L), Time = c("12 weeks", "14 weeks", "13 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks"),
                         Measure = c("CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <2 at end of study", "AISRS - <30% rispetto basline",
                                     "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study"),
                         Age = c(39, 37, 40.4, 32.9, 15.8, 16.5, 37.1, 16.1), Female. = c(42.9, 27, 15.9, 21.4, 15.9, 21.1, 10.4, 21.4), Substance = c("Cocaine", "Cocaine", "Cocaine", "Cannabis", "Mixed", "Alcohol, cannabis, opiate", "Cocaine", "Cannabis, alcohol"),
                         CBT = c("No", "Yes", "Yes", "No", "No", "Yes", "Yes", "Yes")), class = "data.frame", row.names = c(NA, -8L))
#----------------------------Random-effects meta-analysis
pw1 <- metabin(event.e = n_adhdR_1, n.e = n1, event.c = n_adhdR_2, n.c = n2,
               studlab = ID, data = adhd_r, sm = "OR")
pw1
#forestplot
forest(pw1, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Placebo",
       label.left = "Favours placebo", label.right = "Favours active drug", allstudies = F)
#----------------------------------Publication Bias
#Funnel Plot
funnel(pw1, xlab = "Hedges' g")
#Egger's test
#Calculating Effect Size (ES)
adhd_r$n_adhdNR_1 <- (adhd_r$n1 - adhd_r$n_adhdR_1)
adhd_r$n_adhdNR_2 <- (adhd_r$n2 - adhd_r$n_adhdR_2)
adhd_re_data <- escalc(measure = "OR",
                       ai = n_adhdR_1, bi = n_adhdNR_1, n1i = n1,
                       ci = n_adhdR_2, di = n_adhdNR_2, n2i = n2,
                       data = adhd_r)
#Pooling ES
adhd_r_re <- rma(yi = adhd_re_data$yi, vi = adhd_re_data$vi)
adhd_r_re
regtest(adhd_r_re)
#-----Carry out trim-and-fill analysis
adhd_r_taf <- trimfill(adhd_r_re)
adhd_r_taf
#----------------------------------Leave-one-out
adhd_r_loo <- leave1out(adhd_r_re)
adhd_r_loo$or <- exp(adhd_r_loo$estimate)
adhd_r_loo$or_lb <- exp(adhd_r_loo$ci.lb)
adhd_r_loo$or_ub <- exp(adhd_r_loo$ci.ub)
###-------------------------------Meta-regression
#Age
mreg_d_age <- metareg(pw1, ~ Age)
mreg_d_age
#% Female
mreg_d_fem <- metareg(pw1, ~ Female.)
mreg_d_fem
#CBT
mreg_d_cbt <- metareg(pw1, ~ CBT)
mreg_d_cbt

#--------------------------------2. MA ADHD symptoms cont.----------------------------------#
adhd_cont <- structure(list(ID = c("Levin et al. 2015", "Levin et al. 2024", "McRae-Clark et al. 2010", "Riggs et al. 2011", "Thurstone et al. 2010", "Wilens et al. 2008"),
                            t1 = c("AMPH", "AMPH", "ATO", "MPH", "ATO", "ATO"), t2 = c("PBO", "PBO", "PBO", "PBO", "PBO", "PBO"),
                            contrast = c("extended-release mixed amphetamine salts (60mg) - PBO", "extended-release mixed amphetamine salts (80mg) - PBO", "Atomoxetine - PBO", "Osmotic-Release Methylphenidate-PBO", "Atomoxetine-PBO", "Atomoxetine-PBO"),
                            n1 = c(83L, 13L, 24L, 151L, 35L, 72L), n2 = c(43L, 15L, 22L, 152L, 35L, 75L), mean_ADHD_1 = c(18.07, -19.25, 2.63, 20, -18.19, -13.6), sd_ADHD_1 = c(13.79, 10.97, 0.68, 11.91, 14.4, 11.4),
                            mean_ADHD_2 = c(25.78, -15.79, 3.26, 19.4, -19.02, -8.3), sd_ADHD_2 = c(13.94, 9.94, 0.93, 11.95, 15.2, 11.4), Scale = c("AISRS", "AISRS change", "CGI-improvement", "ADHD-RS", "ADHD-RS change", "AISRS change"),
                            Time = c("13 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks"), Age = c(40.4, 32.9, 29.9, 16.5, 16.1, 34.6), Female. = c(15.9, 21.4, 20, 21.1, 21.4, 15),
                            Substance = c("Cocaine", "Cannabis", "Cannabis", "Alcohol, cannabis, opiate", "Cannabis, alcohol", "Alcohol"), CBT = c("Yes", "No", "No", "Yes", "Yes", "No")), row.names = c(NA, -6L), class = "data.frame")
#----------------------------Random-effects meta-analysis
pw2 <- metacont(n.e = n1, mean.e = mean_ADHD_1, sd.e = sd_ADHD_1,
                n.c = n2, mean.c = mean_ADHD_2, sd.c = sd_ADHD_2,
                data = adhd_cont, studlab = ID, sm = "SMD", method.tau = "DL")
pw2
#forestplot
forest(pw2, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Controls",
       label.left = "Favours active drug", label.right = "Favours controls", allstudies = F)
#----------------------------------Publication Bias
#Funnel Plot
funnel(pw2, xlab = "Hedges' g")
#Egger's test
#Calculating Effect Size (ES)
adhdc_re_data <- escalc(measure = "SMD",
                        m1i = mean_ADHD_1, sd1i = sd_ADHD_1, n1i = n1,
                        m2i = mean_ADHD_2, sd2i = sd_ADHD_2, n2i = n2,
                        data = adhd_cont)
#Pooling ES
adhdc_re <- rma(yi = adhdc_re_data$yi, vi = adhdc_re_data$vi)
adhdc_re
regtest(adhdc_re)
#-----Carry out trim-and-fill analysis
adhdc_taf <- trimfill(adhdc_re)
adhdc_taf
#----------------------------------Leave-one-out
adhdc_loo <- leave1out(adhdc_re)
###-------------------------------Meta-regression
#Age
mreg_d_age <- metareg(pw2, ~ Age)
mreg_d_age
#% Female
mreg_d_fem <- metareg(pw2, ~ Female.)
mreg_d_fem
#CBT
mreg_d_cbt <- metareg(pw2, ~ CBT)
mreg_d_cbt

#------------------------3. MA ADHD hyperactivity subscale---------------------------#
hyp <- structure(list(ID = c("Levin et al. 2015", "Wilens et al.2008"), t1 = c("AMPH", "Atomoxetine"), t2 = c("PBO", "PBO"), contrast = c("extended-release mixed amphetamine salts (60mg) - PBO", "Atomoxetine-PBO"),
                      n1 = c(83L, 72L), n2 = c(43L, 75L), mean_ADHD_hyperact_1 = c(14.31, 6.5), sd_ADHD_hyperact_1 = c(13.88, 6), mean_ADHD_hyperact_2 = c(5.42, 3.9), sd_ADHD_hyperact_2 = c(14.92, 5.6),
                      Time = c("13 weeks", "12 weeks"), Measure = c("CAARS-change", "AISRS-change")), row.names = c(NA, -2L), class = "data.frame")
#----------------------------Random-effects meta-analysis
pw3 <- metacont(n.e = n1, mean.e = mean_ADHD_hyperact_1, sd.e = sd_ADHD_hyperact_1,
                n.c = n2, mean.c = mean_ADHD_hyperact_2, sd.c = sd_ADHD_hyperact_2,
                data = hyp, studlab = ID, sm = "SMD", method.tau = "DL")
pw3
#forestplot
forest(pw3, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Controls",
       label.left = "Favours placebo", label.right = "Favours active drug", allstudies = F)

#-------------------------4. MA ADHD inattentive subscale-------------------------#
inatt <- structure(list(ID = c("Levin et al. 2015", "Wilens et al.2008"), t1 = c("extended-release mixed amphetamine salts (60mg)", "Atomoxetine"), t2 = c("PBO", "PBO"),
                        contrast = c("extended-release mixed amphetamine salts (60mg) - PBO", "Atomoxetine-PBO"), n1 = c(83L, 72L), n2 = c(43L, 75L),
                        mean_ADHD_inatt_1 = c(14.86, 7.2), sd_ADHD_inatt_1 = c(15.28, 6.2), mean_ADHD_inatt_2 = c(4.03, 4.4), sd_ADHD_inatt_2 = c(11.66, 6.7),
                        Time = c("13 weeks", "12 weeks"), Measure = c("CAARS-change", "AISRS-change")), row.names = c(NA, -2L), class = "data.frame")
#----------------------------Random-effects meta-analysis
pw4 <- metacont(n.e = n1, mean.e = mean_ADHD_inatt_1, sd.e = sd_ADHD_inatt_1,
                n.c = n2, mean.c = mean_ADHD_inatt_2, sd.c = sd_ADHD_inatt_2,
                data = inatt, studlab = ID, sm = "SMD", method.tau = "DL")
pw4
#forestplot
forest(pw4, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Controls",
       label.left = "Favours placebo", label.right = "Favours active drug", allstudies = F)

#-------------------------5. NMA ADHD parent rating-------------------------#
adhd_parent <- structure(list(ID = c("Riggs et al. 2004", "Riggs et al. 2011", "Thurstone et al. 2010"),
                              t1 = c("Pemoline", "Osmotic-Release Methylphenidate", "Atomoxetine"), t2 = c("PBO", "PBO", "PBO"), contrast = c("Pemoline-PBO", "Osmotic-Release Methylphenidate-PBO", "Atomoxetine-PBO"),
                              n1 = c(33L, 84L, 35L), n2 = c(33L, 68L, 35L), mean_ADHD_1 = c(68.5, 24, -13.82), sd_ADHD_1 = c(21, 11.8, 13.91), mean_ADHD_2 = c(71.2, 30.9, -8.82), sd_ADHD_2 = c(14.18, 13, 16.36),
                              Time = c("12 weeks", "16 weeks", "12 weeks"), Measure = c("CHI", "ADHD-RS", "ADHD-RS change")), row.names = c(NA, -3L), class = "data.frame")
#----------------------------Random-effects meta-analysis
pw5 <- metacont(n.e = n1, mean.e = mean_ADHD_1, sd.e = sd_ADHD_1,
                n.c = n2, mean.c = mean_ADHD_2, sd.c = sd_ADHD_2,
                data = adhd_parent, studlab = ID, sm = "SMD", method.tau = "DL")
pw5
#forestplot
forest(pw5, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Controls",
       label.left = "Favours active drug", label.right = "Favours placebo", allstudies = F)

#---------------------------------6. MA SUD cont.--------------------------------#
sud_cont <- structure(list(ID = c("Levin et al. 2007", "McRae-Clark et al. 2010", "Riggs et al. 2004", "Riggs et al. 2011", "Szobot et al. 2008", "Thurstone et al. 2010"),
                           t1 = c("sustained-release MPH", "Atomoxetine", "Pemoline", "Osmotic-Release Methylphenidate", "methylphenidate-SODAS", "Atomoxetine"), t2 = c("PBO", "PBO", "PBO", "PBO", "PBO", "PBO"),
                           contrast = c("sustained-release MPH- PBO", "Atomoxetine - PBO", "Pemoline-PBO", "Osmotic-Release Methylphenidate-PBO", "methylphenidate-SODAS-PBO", "Atomoxetine-PBO"),
                           n1 = c(53L, 24L, 35L, 151L, 16L, 35L), n2 = c(53L, 22L, 34L, 152L, 16L, 35L), mean_negative_SUD_1 = c(0.73, 60.1, 12.1, 8.2, 5.56, -5.78), 
                           sd_negative_SUD_1 = c(0.29, 31.5, 11.3, 8.78, 2.03, 10.4), mean_negative_SUD_2 = c(0.7, 68.1, 13.7, 9.1, 6, -2.24), sd_negative_SUD_2 = c(0.29, 31.3, 11.5, 9.12, 2.1, 10.3),
                           Time = c("12 weeks", "12 weeks", "12 weeks", "12 weeks", "6 weeks", "12 weeks"), Measure = c("week positive any drugs", "%days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs - change"),
                           Age = c(37, 29.9, 15.8, 16.5, 17.4, 16.1), Female. = c(27, 20, 15.9, 21.1, NA, 21.4), Substance = c("Cocaine", "Cannabis", "Mixed", "Alcohol, cannabis, opiate", "Cannabis", "Cannabis, alcohol"),
                           CBT = c("Yes", "No", "No", "Yes", "No", "Yes")), row.names = c(NA, -6L), class = "data.frame")
#----------------------------Random-effects meta-analysis
pw6 <- metacont(n.e = n1, mean.e = mean_negative_SUD_1, sd.e = sd_negative_SUD_1,
                n.c = n2, mean.c = mean_negative_SUD_2, sd.c = sd_negative_SUD_2,
                data = sud_cont, studlab = ID, sm = "SMD", method.tau = "DL")
pw6
#forestplot
forest(pw6, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Controls",
       label.left = "Favours active drug", label.right = "Favours placebo", allstudies = F)
#----------------------------------Publication Bias
#Funnel Plot
funnel(pw6, xlab = "Hedges' g")
#Egger's test
#Calculating Effect Size (ES)
sudc_re_data <- escalc(measure = "SMD",
                       m1i = mean_negative_SUD_1, sd1i = sd_negative_SUD_1, n1i = n1,
                       m2i = mean_negative_SUD_2, sd2i = sd_negative_SUD_2, n2i = n2,
                       data = sud_cont)
#Pooling ES
sudc_re <- rma(yi = sudc_re_data$yi, vi = sudc_re_data$vi)
sudc_re
regtest(sudc_re)
#-----Carry out trim-and-fill analysis
sudc_taf <- trimfill(sudc_re)
sudc_taf
#----------------------------------Leave-one-out
sudc_loo <- leave1out(sudc_re)
###-------------------------------Meta-regression
#Age
mreg_d_age <- metareg(pw6, ~ Age)
mreg_d_age
#% Female
mreg_d_fem <- metareg(pw6, ~ Female.)
mreg_d_fem
#CBT
mreg_d_cbt <- metareg(pw6, ~ CBT)
mreg_d_cbt

#----------------------------7. MA SUD binary--------------------------------#
sud_bin <- structure(list(ID = c("Levin et al. 2006", "Levin et al. 2007", "Levin et al. 2015", "Levin et al. 2024", "Schubiner et al. 2002"),
                          t1 = c("MPH", "sustained-release MPH", "extended-release mixed amphetamine salts (60mg)", "extended-release mixed amphetamine salts (80mg)", "methylphenidate"), t2 = c("PBO", "PBO", "PBO", "PBO", "PBO"),
                          contrast = c("sustained-release MPH- PBO", "sustained-release MPH- PBO", "extended-release mixed amphetamine salts (60mg) - PBO", "extended-release mixed amphetamine salts (80mg) - PBO", "methylphenidate-PBO"),
                          n1 = c(65L, 53L, 83L, 13L, 24L), n2 = c(33L, 53L, 43L, 15L, 24L), n_SUD_negative_1 = c(5L, 8L, 20L, 2L, 12L), n_SUD_negative_2 = c(5L, 9L, 3L, 0L, 10L),
                          Time = c("12 weeks", "14 weeks", "13 weeks", "12 weeks", "12 weeks"),
                          Measure = c("Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", "Abstinence from using cocaine in the past 3 weeks", "Abstinence from using any drugs in the past 2 weeks", "N cocaine free urine drug test"),
                          Age = c(39, 37, 40.4, 32.9, 37.1), Female. = c(42.9, 27, 15.9, 21.4, 10.4), Substance = c("Cocaine", "Cocaine", "Cocaine", "Cannabis", "Cocaine"), CBT = c("No", "Yes", "Yes", "No", "Yes")), class = "data.frame", row.names = c(NA, -5L))
#----------------------------Random-effects meta-analysis
pw7 <- metabin(event.e = n_SUD_negative_1, n.e = n1, event.c = n_SUD_negative_2, n.c = n2,
               studlab = ID, data = sud_bin, sm = "OR")
pw7
#forestplot
forest(pw7, layout = "RevMan5", digits.sd = 2, random = T, fixed = F,
       label.e = "Active drug", label.c = "Placebo",
       label.left = "Favours placebo", label.right = "Favours active drug", allstudies = F)
#----------------------------------Publication Bias
#Funnel Plot
funnel(pw7, xlab = "Hedges' g")
#Egger's test
#Calculating Effect Size (ES)
sud_bin$n_sudNR_1 <- (sud_bin$n1 - sud_bin$n_SUD_negative_1)
sud_bin$n_sudNR_2 <- (sud_bin$n2 - sud_bin$n_SUD_negative_2)
sud_re_data <- escalc(measure = "OR",
                      ai = n_SUD_negative_1, bi = n_sudNR_1, n1i = n1,
                      ci = n_SUD_negative_2, di = n_sudNR_2, n2i = n2,
                      data = sud_bin)
#Pooling ES
sudd_r_re <- rma(yi = sud_re_data$yi, vi = sud_re_data$vi)
sudd_r_re
regtest(sudd_r_re)
#-----Carry out trim-and-fill analysis
sud_r_taf <- trimfill(sudd_r_re)
sud_r_taf
#----------------------------------Leave-one-out
sud_r_loo <- leave1out(sudd_r_re)
sud_r_loo$or <- exp(sud_r_loo$estimate)
sud_r_loo$or_lb <- exp(sud_r_loo$ci.lb)
sud_r_loo$or_ub <- exp(sud_r_loo$ci.ub)
###-------------------------------Meta-regression
#Age
mreg_d_age <- metareg(pw7, ~ Age)
mreg_d_age
#% Female
mreg_d_fem <- metareg(pw7, ~ Female.)
mreg_d_fem
#CBT
mreg_d_cbt <- metareg(pw7, ~ CBT)
mreg_d_cbt
