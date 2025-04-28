# R code for replicating the dataset and the analyses in the paper #
# Pharmacological treatments for ADHD comorbid with substance use disorder: a systematic review and network meta-analysis #
# R code by Mattia Marchi (mattiamarchimd@gmail.com) 
# April 28, 2025

#----------------------------1. NMA ADHD treatment response--------------------------------#
adhd_r <- structure(list(ID = c("Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2007", "Levin et al. 2007", "Levin et al. 2015", "Levin et al. 2015", "Levin et al. 2024", "Levin et al. 2024", "Riggs et al. 2004", "Riggs et al. 2004", 
                                "Riggs et al. 2011", "Riggs et al. 2011", "Schubiner et al. 2002", "Schubiner et al. 2002", "Thurstone et al. 2010", "Thurstone et al. 2010"),
                         t = c("MPH", "BPR", "PBO", "MPH", "PBO", "AMPH", "PBO", "AMPH", "PBO", "PEM", "PBO", "MPH", "PBO", "MPH", "PBO", "ATO", "PBO"),
                         n = c(32L, 33L, 33L, 53L, 53L, 83L, 43L, 12L, 15L, 35L, 34L, 151L, 152L, 24L, 24L, 35L, 35L), r = c(11L, 16L, 15L, 25L, 29L, 31L, 5L, 10L, 10L, 11L, 4L, 39L, 35L, 18L, 5L, 17L, 20L),
                         Time = c("12 weeks", "12 weeks", "12 weeks", "14 weeks", "14 weeks", "13 weeks", "13 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks"),
                         Measure = c("CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", 
                                     "CGI - CGI ADHD improvement rating <2 at end of study", "CGI - CGI ADHD improvement rating <2 at end of study", "AISRS - <30% than baseline", "AISRS - <30% than baseline", "CGI - CGI ADHD improvement rating <3 at end of study", 
                                     "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study",
                                     "CGI - CGI ADHD improvement rating <3 at end of study", "CGI - CGI ADHD improvement rating <3 at end of study")), class = "data.frame", row.names = c(NA, -17L))
#Organize data/calculate pairwise comparisons
pw1 <- pairwise(treat = t, n = n, event = r,
                studlab = ID, data = adhd_r, sm = "OR")
#Perform standard NMA
net1 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw1,
                common = FALSE, ref = "PBO")
dep_net <- netgraph(net1, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                    cex.points = table(adhd_r$t)*2, col = 1, number.of.studies = TRUE,
                    pos.number.of.studies = 0.5)
print(summary(net1))
forest(net1, label.left = "Favours placebo", label.right = "Favours active drug", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
#To make all the contrasts
plot(net1, ref = c("AMPH", "ATO", "BPR", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net1, small.values = "bad")
set.seed(1909)
ran1 <- rankogram(net1, small.values = "bad")
plot(ran1)#, cumulative.rankprob = T)
set.seed(1909)
sucra1 <- netrank(ran1)
netleague(net1, seq = netrank(net1), ci = T)
#Decompose heterogenity
decomp.design(net1)
netsplit(net1)

#--------------------------------2. NMA ADHD symptoms cont.----------------------------------#
adhd_cont <- structure(list(ID = c("Levin et al. 2015", "Levin et al. 2015", "Levin et al. 2024", "Levin et al. 2024", "McRae-Clark et al. 2010", "McRae-Clark et al. 2010", "Riggs et al. 2011", "Riggs et al. 2011", 
                                   "Thurstone et al. 2010", "Thurstone et al. 2010", "Wilens et al.2008", "Wilens et al.2008"), t1 = c("AMPH", "PBO", "AMPH", "PBO", "ATO", "PBO", "MPH", "PBO", "ATO", "PBO", "ATO", "PBO"),
                            n = c(83L, 43L, 13L, 15L, 24L, 22L, 151L, 152L, 35L, 35L, 72L, 75L), mean = c(18.07, 25.78, -19.25, -15.79, 2.63, 3.26, 20, 19.4, -18.19, -19.02, -13.6, -8.3), sd = c(13.79, 13.94, 10.97, 9.94, 0.68, 0.93, 11.91, 11.95, 14.4, 15.2, 11.4, 11.4),
                            Scale = c("AISRS", "AISRS", "AISRS change (from baseline)", "AISRS change (from baseline)", "CGI-improvement", "CGI-improvement", "ADHD-RS", "ADHD-RS", "ADHD-RS change (pre-post=18.19)", "ADHD-RS change (pre-post=19.2)", "AISRS change (from baseline= -13.6)", "AISRS change (from baseline= -8.3)"),
                            Time = c("13 weeks", "13 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks")), row.names = c(NA, -12L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw2 <- pairwise(treat = t1, n = n, mean = mean, sd = sd,
                studlab = ID, data = adhd_cont, sm = "SMD")
net2 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw2,
                common = FALSE, ref = "PBO")
anx_net <- netgraph(net2, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                    cex.points = table(adhd_cont$t1)*2, col = 1, number.of.studies = TRUE,
                    pos.number.of.studies = 0.5)
print(summary(net2))
forest(net2, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net2)
#To make all the contrasts
plot(net2, ref = c("AMPH", "ATO", "MPH", "PBO"))
#Ranking treatments
netrank(net2, small.values = "good")
set.seed(1909)
ran2 <- rankogram(net2, small.values = "good")
plot(ran2)
set.seed(1909)
sucra2 <- netrank(ran2)
netleague(net2, seq = netrank(net2), ci = F)
#Decompose heterogenity
decomp.design(net2)
netsplit(net2)

#------------------------3. NMA ADHD hyperactivity subscale---------------------------#
hyp <- structure(list(ID = c("Levin et al. 2015", "Levin et al. 2015", "Wilens et al.2008", "Wilens et al.2008"), t1 = c("AMPH", "PBO", "ATO", "PBO"),
                      n = c(83L, 43L, 72L, 75L), mean = c(-14.31, -5.42, -6.5, -3.9), sd = c(13.88, 14.92, 6, 5.6),
                      Time = c("13 weeks", "13 weeks", "12 weeks", "12 weeks"), Measure = c("CAARS-change", "CAARS-change", "AISRS-change", "AISRS-change")), row.names = c(NA, -4L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw6 <- pairwise(treat = t1, n = n, mean = mean, sd = sd,
                studlab = ID, data = hyp, sm = "SMD")
#Perform standard NMA
net6 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw6,
                common = FALSE, ref = "PBO")
da_net <- netgraph(net6, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                   cex.points = table(hyp$t1)*2, col = 1, number.of.studies = TRUE,
                   pos.number.of.studies = 0.5)
print(summary(net6))
forest(net6, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net6)
#To make all the contrasts
plot(net6, ref = c("ATO", "AMPH", "PBO"))
#Ranking treatments
netrank(net6, small.values = "good")
set.seed(1909)
ran6 <- rankogram(net6, small.values = "good")
plot(ran6)
set.seed(1909)
sucra6 <- netrank(ran6)
netleague(net6, seq = netrank(net6), ci = F)
#Decompose heterogenity
decomp.design(net6)
netsplit(net6)

#-------------------------4. NMA ADHD inattentive subscale-------------------------#
inatt <- structure(list(ID = c("Levin et al. 2015", "Levin et al. 2015", "Wilens et al.2008", "Wilens et al.2008"), t1 = c("AMPH", "PBO", "ATO", "PBO"),
                        n = c(83L, 43L, 72L, 75L), mean = c(-14.86, -4.03, -7.2, -4.4), sd = c(15.28, 11.66, 6.2, 6.7),
                        Time = c("13 weeks", "13 weeks", "12 weeks", "12 weeks"), Measure = c("CAARS-change", "CAARS-change", "AISRS-change", "AISRS-change")), row.names = c(NA, -4L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw7 <- pairwise(treat = t1, n = n, mean = mean, sd = sd,
                studlab = ID, data = inatt, sm = "SMD")
#Perform standard NMA
net7 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw7,
                common = FALSE, ref = "PBO")
qol_net <- netgraph(net7, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                    cex.points = table(inatt$t1)*2, col = 1, number.of.studies = TRUE,
                    pos.number.of.studies = 0.5)
print(summary(net7))
forest(net7, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net7)
#To make all the contrasts
plot(net7, ref = c("ATO", "AMPH", "PBO"))
#Ranking treatments
netrank(net7, small.values = "good")
set.seed(1909)
ran7 <- rankogram(net7, small.values = "good")
plot(ran7)
set.seed(1909)
sucra7 <- netrank(ran7)
netleague(net7, seq = netrank(net7), ci = F)
#Decompose heterogenity
decomp.design(net7)
netsplit(net7)

#-------------------------5. NMA ADHD parent rating-------------------------#
adhd_parent <- structure(list(ID = c("Riggs et al. 2004", "Riggs et al. 2004", "Riggs et al. 2011", "Riggs et al. 2011", "Thurstone et al. 2010", "Thurstone et al. 2010"), t1 = c("PEM", "PBO", "MPH", "PBO", "ATO", "PBO"),
                              n = c(33L, 33L, 84L, 68L, 35L, 35L), mean = c(68.5, 71.2, 24, 30.9, -13.82, -8.82), sd = c(21, 14.18, 11.8, 13, 13.91, 16.36), Time = c("12 weeks", "12 weeks", "16 weeks", "16 weeks", "12 weeks", "12 weeks"),
                              Measure = c("CHI", "CHI", "ADHD-RS", "ADHD-RS", "ADHD-RS change", "ADHD-RS change")), row.names = c(NA, -6L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw8 <- pairwise(treat = t1, n = n, mean = mean, sd = sd,
                studlab = ID, data = adhd_parent, sm = "SMD")
#Perform standard NMA
net8 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw8,
                common = FALSE, ref = "PBO")
parent_net <- netgraph(net8, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                       cex.points = table(adhd_parent$t1)*2, col = 1, number.of.studies = TRUE,
                       pos.number.of.studies = 0.5)
print(summary(net8))
forest(net8, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net8)
#To make all the contrasts
plot(net8, ref = c("ATO", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net8, small.values = "good")
set.seed(1909)
ran8 <- rankogram(net8, small.values = "good")
plot(ran8)
set.seed(1909)
sucra8 <- netrank(ran8)
netleague(net8, seq = netrank(net8), ci = F)
#Decompose heterogenity
decomp.design(net8)
netsplit(net8)

#---------------------------------6. NMA SUD cont.--------------------------------#
sud_cont <- structure(list(ID = c("Levin et al. 2007", "Levin et al. 2007", "McRae-Clark et al. 2010", "McRae-Clark et al. 2010", "Riggs et al. 2004", "Riggs et al. 2004", "Riggs et al. 2011", "Riggs et al. 2011", 
                                  "Szobot et al. 2008", "Szobot et al. 2008", "Thurstone et al. 2010", "Thurstone et al. 2010"), t1 = c("MPH", "PBO", "ATO", "PBO", "PEM", "PBO", "MPH", "PBO", "MPH", "PBO", "ATO", "PBO"),
                           n = c(53L, 53L, 24L, 22L, 35L, 34L, 151L, 152L, 16L, 16L, 35L, 35L), mean = c(0.73, 0.7, 60.1, 68.1, 12.1, 13.7, 8.2, 9.1, 5.56, 6, -5.78, -2.24), sd = c(0.29, 0.29, 31.5, 31.3, 11.3, 11.5, 8.78, 9.12, 2.03, 2.1, 10.4, 10.3),
                           Time = c("12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks", "6 weeks", "6 weeks", "12 weeks", "12 weeks"),
                           Measure = c("week positive any drugs", "week positive any drugs", "%days using any drugs", "%days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs", "days using any drugs - change", "days using any drugs - change")),
                      row.names = c(NA, -12L), class = "data.frame")
#Organize data/calculate pairwise comparisons
pw9 <- pairwise(treat = t1, n = n, mean = mean, sd = sd,
                studlab = ID, data = sud_cont, sm = "SMD")
#Perform standard NMA
net9 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw9,
                common = FALSE, ref = "PBO")
sudcont_net <- netgraph(net9, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                        cex.points = table(sud_cont$t1)*2, col = 1, number.of.studies = TRUE,
                        pos.number.of.studies = 0.5)
print(summary(net9))
forest(net9, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Std. Mean Difference
       Random, 95% CI")
plot(net9)
#To make all the contrasts
plot(net9, ref = c("ATO", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net9, small.values = "good")
set.seed(1909)
ran9 <- rankogram(net9, small.values = "good")
plot(ran9)
set.seed(1909)
sucra9 <- netrank(ran9)
netleague(net9, seq = netrank(net9), ci = F)
#Decompose heterogenity
decomp.design(net9)
netsplit(net9)

#----------------------------7. NMA SUD binary--------------------------------#
sud_bin <- structure(list(ID = c("Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2007", "Levin et al. 2007", "Levin et al. 2015", "Levin et al. 2015", "Levin et al. 2024", 
                                 "Levin et al. 2024", "Schubiner et al. 2002", "Schubiner et al. 2002"), t1 = c("MPH", "BPR", "PBO", "MPH", "PBO", "AMPH", "PBO", "AMPH", "PBO", "MPH", "PBO"),
                          n = c(32L, 33L, 33L, 53L, 53L, 83L, 43L, 13L, 15L, 24L, 24L), r = c(3L, 2L, 5L, 8L, 9L, 20L, 3L, 2L, 0L, 12L, 10L),
                          Time = c("12 weeks", "12 weeks", "12 weeks", "14 weeks", "14 weeks", "13 weeks", "13 weeks", "12 weeks", "12 weeks", "12 weeks", "12 weeks"),
                          Measure = c("Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", 
                                      "Abstinence from using cocaine in the past 3 weeks", "Abstinence from using cocaine in the past 3 weeks", "Abstinence from using any drugs in the past 2 weeks", "Abstinence from using any drugs in the past 2 weeks", "N cocaine free urine drug test", "N cocaine free urine drug test")),
                     class = "data.frame", row.names = c(NA, -11L))
#Organize data/calculate pairwise comparisons
pw10 <- pairwise(treat = t1, n = n, event = r,
                 studlab = ID, data = sud_bin, sm = "OR")
#Perform standard NMA
net10 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw10,
                 common = FALSE, ref = "PBO")
sudbin_net <- netgraph(net10, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                       cex.points = table(sud_bin$t)*2, col = 1, number.of.studies = TRUE,
                       pos.number.of.studies = 0.5)
print(summary(net10))
forest(net10, label.left = "Favours placebo", label.right = "Favours active drug", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
#To make all the contrasts
plot(net10, ref = c("AMPH", "BPR", "MPH", "PBO"))
#Ranking treatments
netrank(net10, small.values = "bad")
set.seed(1909)
ran10 <- rankogram(net10, small.values = "bad")
plot(ran10)#, cumulative.rankprob = T)
set.seed(1909)
sucra10 <- netrank(ran10)
netleague(net10, seq = netrank(net10), ci = T)
#Decompose heterogenity
decomp.design(net10)
netsplit(net10)

#--------------------------------------8. NMA Tolerability-------------------------------------#
tol <- structure(list(ID = c("Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2006", "Levin et al. 2007", "Levin et al. 2007", "Levin et al. 2015", "Levin et al. 2015", "Levin et al. 2024", 
                             "Levin et al. 2024", "McRae-Clark et al. 2010", "McRae-Clark et al. 2010", "Riggs et al. 2004", "Riggs et al. 2004", "Riggs et al. 2011", "Riggs et al. 2011", "Schubiner et al. 2002", "Schubiner et al. 2002", 
                             "Szobot et al. 2008", "Szobot et al. 2008", "Thurstone et al. 2010", "Thurstone et al. 2010", "Wilens et al.2008", "Wilens et al.2008"),
                      t1 = c("MPH", "BPR", "PBO", "MPH", "PBO", "AMPH", "PBO", "AMPH", "PBO", "ATO", "PBO", "PEM", "PBO", "MPH", "PBO", "MPH", "PBO", "MPH", "PBO", "ATO", "PBO", "ATO", "PBO"),
                      n = c(32L, 33L, 33L, 53L, 53L, 83L, 43L, 13L, 15L, 24L, 22L, 35L, 34L, 151L, 152L, 24L, 24L, 16L, 16L, 35L, 35L, 72L, 75L),
                      r_doany = c(11L, 10L, 8L, 30L, 29L, 19L, 14L, 3L, 0L, 15L, 15L, 14L, 16L, 33L, 43L, 13L, 10L, 2L, 0L, 3L, 2L, 40L, 27L),
                      r_doae = c(1L, 3L, 3L, 3L, 4L, 0L, 0L, 1L, 0L, 2L, 2L, 2L, 3L, 1L, 4L, 0L, 1L, 2L, 0L, 1L, 1L, 7L, 2L),
                      r_ae = c(1L, 0L, 2L, 1L, 1L, 82L, 16L, 13L, 7L, 19L, 16L, 2L, 3L, 4L, 7L, 17L, 19L, 2L, 0L, 4L, 7L, 31L, 7L)),
                 class = "data.frame", row.names = c(NA, -23L))
#----Dropout any cause
pw3 <- pairwise(treat = t1, n = n, event = r_doany, studlab = ID, data = tol, sm = "OR")
net3 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw3,
                common = FALSE, ref = "PBO")
anycause_net <- netgraph(net3, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                         cex.points = table(tol$t1), col = 1, number.of.studies = TRUE,
                         pos.number.of.studies = 0.5)
print(summary(net3))
forest(net3, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
plot(net3)
#To make all the contrasts
plot(net3, ref = c("AMPH", "ATO", "BPR", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net3, small.values = "good")
set.seed(1909)
ran3 <- rankogram(net3, small.values = "good")
plot(ran3)
set.seed(1909)
sucra3 <- netrank(ran3)
netleague(net3, seq = netrank(net3), ci = F)
#Decompose heterogenity
decomp.design(net3)
netsplit(net3)

#----Dropout severe adverse reactions
pw4 <- pairwise(treat = t1, n = n, event = r_doae, studlab = ID, data = tol, sm = "OR", allstudies = T)
net4 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw4,
                common = FALSE, ref = "PBO")
severeae_net <- netgraph(net4, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                         cex.points = table(tol$t1), col = 1, number.of.studies = TRUE,
                         pos.number.of.studies = 0.5)
print(summary(net4))
forest(net4, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
plot(net4)
#To make all the contrasts
plot(net4, ref = c("AMPH", "ATO", "BPR", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net4, small.values = "good")
set.seed(1909)
ran4 <- rankogram(net4, small.values = "good")
plot(ran4)
set.seed(1909)
sucra4 <- netrank(ran4)
netleague(net4, seq = netrank(net4), ci = F)
#Decompose heterogenity
decomp.design(net4)
netsplit(net4)

#----N adverse effects
pw5 <- pairwise(treat = t1, n = n, event = r_ae, studlab = ID, data = tol, sm = "OR")
net5 <- netmeta(TE, seTE, treat1, treat2, studlab, data = pw5,
                common = FALSE, ref = "PBO")
ae_net <- netgraph(net5, plastic = FALSE, multiarm = FALSE, points = TRUE, seq = "optimal",
                   cex.points = table(tol$t1), col = 1, number.of.studies = TRUE,
                   pos.number.of.studies = 0.5)
print(summary(net5))
forest(net5, label.left = "Favours active drug", label.right = "Favours placebo", sortvar = TE,
       smlab = "Odds Ratio
       Random, 95% CI")
plot(net5)
#To make all the contrasts
plot(net5, ref = c("AMPH", "ATO", "BPR", "MPH", "PEM", "PBO"))
#Ranking treatments
netrank(net5, small.values = "good")
set.seed(1909)
ran5 <- rankogram(net5, small.values = "good")
plot(ran5)
set.seed(1909)
sucra5 <- netrank(ran5)
netleague(net5, seq = netrank(net5), ci = F)
#Decompose heterogenity
decomp.design(net5)
netsplit(net5)