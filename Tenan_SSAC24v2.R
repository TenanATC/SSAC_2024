##SSAC24 script primarily written by Matthew Tenan but reviewed and edited by Andrew Vigotsky

#This is designed to run on a Windows PC with 16 cores, some parallel processing code
#may not work appropriately if you are running Linux or do not have sufficient cores/memory
#These are relatively easily solved, but may require a few 'tweaks' to run without errors

#These are the libraries for the 'core analyses'
#We've already done the hard work of aggregating the various data sources together and
#posted that data to github
library(mgcv)
library(mvtnorm)
library(matrixStats)
library(itsadug)
library(dplyr)
library(WeightIt)
library(cobalt)
library(pbapply)
library(lme4)
library(MASS)
library(r2glmm)
library(ggplot2)
library(ggpubr)



cores <- 16 #if you have a Windows machine and have less than 16 cores, you can just change this

SSAC_dat <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/TENAN_SSAC2024_dat.rds")))
SSAC_dat$elements <- as.factor(SSAC_dat$elements)

#zero-ing out weather variables when elements are '0'
#my guess is that before we submit this for final publication at a peer-reviewed journal, we'll
#change the humidity and temperature numbers to not be '0' because obviously inside games
#aren't played at those actual numbers... but we'll need to pick an empirical number (or number process) that has face validity
SSAC_dat$precip <- ifelse(SSAC_dat$elements == 1, SSAC_dat$precip, 0)
SSAC_dat$humidity_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$humidity_high, 0)
SSAC_dat$temp_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$temp_high, 0)
SSAC_dat$wind_speed <- ifelse(SSAC_dat$elements == 1, SSAC_dat$wind_speed, 0)


#Demonstrating the need for cluster-level propensity scores
#Below code shows the team-level correlation for the treatments.
#While they're relatively low (0.16 and 0.29) they're definitely high enough to bias an analysis,
#so we'll use multilevel model propensity scores
#On the other hand, the actual analysis, has very low ICC (0.015) so we won't do an MLM for final analysis
#Logit based on ICC method from Hox JJ, Moerbeek M, van de Schoot R (2018). 
#Multilevel Analysis: Techniques and Applications. Taylor and Francis. ISBN 9781138121362.  p. 107

lmer_varcor_time <- VarCorr(lmer(time_hrs_round ~ 1 + (1|away_team), data = SSAC_dat))
time_lmer <- as.data.frame(print(lmer_varcor_time, comp=c('Variance')))
time_lmer$vcov[1] / (time_lmer$vcov[1] + time_lmer$vcov[2])

lmer_varcor_tz <- VarCorr(lmer(away_tz_diff ~ 1 + (1|away_team), data = SSAC_dat))
tz_lmer <- as.data.frame(print(lmer_varcor_tz, comp=c('Variance')))
tz_lmer$vcov[1] / (tz_lmer$vcov[1] + tz_lmer$vcov[2])

glmer_logit <- glmer(spread_beat_away ~ 1 + (1|away_team), data = SSAC_dat, family = binomial)
summ_glmer_logit <- summary(glmer_logit)
summ_glmer_logit$varcor$away_team[[1]] / (summ_glmer_logit$varcor$away_team[[1]] + (pi^2 / 3))


#Multilevel model propensity score weights as previously described using both  
#cluster mean (opt1) and marginal (opt2) stabilized weights from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5157938/
#NOTE: These are multivariate PSWs, whereas above reference was for univariate, so our method is a generalization


num.mod.time <-gam(list(time_hrs_round ~ s(away_team, bs= 're'), away_tz_diff ~ s(away_team, bs= 're')),
                   family= mvn(d=2), data = SSAC_dat , method = 'REML')  # numerator model for option 1

#multivariate denominator model
den.mod.time <- gam(list(time_hrs_round ~ s(humidity_high, k=13) + s(away_team, bs= 're')+  grass_type + 
                           s(precip, k=13) +  s(temp_high, k=13) + s(wind_speed, k= 13),
                         away_tz_diff ~ s(humidity_high, k=13) + s(away_team, bs= 're')+  grass_type +  
                           s(precip, k=13) +  s(temp_high, k=13) + s(wind_speed, k=13)),
                    family= mvn(d=2), data= SSAC_dat, method = 'REML')

#create cluster-mean multivariate normal density for treatments
mvt_dat <- as.matrix(cbind(SSAC_dat$time_hrs_round, SSAC_dat$away_tz_diff))
num.p_op1 <- dmvnorm(x= mvt_dat, mean = colMeans(fitted(num.mod.time)),
                     sigma = solve(crossprod(num.mod.time$family$data$R))) 
#create marginal multivariate normal density for treatments, also used to assess covariate balance (line 130 of code)
num.p_op2 <- dmvnorm(x= mvt_dat, mean = colMeans(mvt_dat), sigma = cov(mvt_dat)) # numerator calculation for option 2

#create multivariate denominator
den.p <- dmvnorm(x= mvt_dat, mean = colMeans(fitted(den.mod.time)), sigma = solve(crossprod(den.mod.time$family$data$R))) 

#Now the actual IPW
Tps_MLM_op1 <- num.p_op1/den.p
Tps_MLM_op2 <- num.p_op2/den.p

#Looking at some diagnostics for assessing positivity assumption in the IPW,
#though I'm a little skeptical these mean what they purport to show according to Cole & Hernan
mean(Tps_MLM_op1)
sd(Tps_MLM_op1)
range(Tps_MLM_op1)
mean(Tps_MLM_op2)
sd(Tps_MLM_op2)
range(Tps_MLM_op2)

SSAC_dat$ps1_opt1 <- Tps_MLM_op1
SSAC_dat$ps1_opt2 <- Tps_MLM_op2

#confirming no major issues with the denominator propensity score
#This generally looks very good, remember there is a stochastic element to the basis checking
gam.check(den.mod.time)



##Having fully developed the IPW, we need to validate that is orthogonalizes the 
##covariates to the treatment and "balances" the covariates effectively
#We are unaware of widely accepted methods to check the IPW for the above two things
#in the case of two continuous effect-modification treatments so we have taken two approaches:
#1) demonstrate that the IPW decreases the 'fit' of the effect-modification treatments on each covariate
#by confirming that weighted models have higher AICs than unweighted models and 
#2) Demonstrate that the adjusted relationship between the effect-modification treatments and the
#each individual treatment has an R-squared below 0.10, the common threshold used for univariate
#correlation analysis of covariate balance

#Together, these give us a reasonable level of confidence in our IPW

precipMLM <- lmer(precip ~ away_tz_diff + time_hrs_round:away_tz_diff + (1|away_team), data = SSAC_dat)
precipMLM_wt <- lmer(precip ~ away_tz_diff + time_hrs_round:away_tz_diff + (1|away_team), data = SSAC_dat, weights = ps1_opt1)
aic_precip <- extractAIC(precipMLM)[2]
aic_precipwt <- extractAIC(precipMLM_wt)[2]
(precip_beta <- r2beta(precipMLM, method = 'sgv', partial = F)$Rsq)
(precip_betawt <- r2beta(precipMLM_wt, method = 'sgv', partial = F)$Rsq)

humidMLM <- lmer(humidity_high ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat)
humidMLM_wt <- lmer(humidity_high ~away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat, weights = ps1_opt1)
aic_humid <- extractAIC(humidMLM)[2]
aic_humidwt <- extractAIC(humidMLM_wt)[2]
(humid_beta <- r2beta(humidMLM, method = 'sgv', partial = F)$Rsq)
(humid_betawt <- r2beta(humidMLM_wt, method = 'sgv', partial = F)$Rsq)

tempMLM <- lmer(temp_high ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat)
tempMLM_wt <- lmer(temp_high ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat, weights = ps1_opt1)
aic_temp <- extractAIC(tempMLM)[2]
aic_tempwt <- extractAIC(tempMLM_wt)[2]
(temp_beta <- r2beta(tempMLM, method = 'sgv', partial = F)$Rsq)
(temp_betawt <- r2beta(tempMLM_wt, method = 'sgv', partial = F)$Rsq)

windMLM <- lmer(wind_speed ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat)
windMLM_wt <- lmer(wind_speed ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat, weights = ps1_opt1)
aic_wind <- extractAIC(windMLM)[2]
aic_windwt <- extractAIC(windMLM_wt)[2]
(wind_beta <- r2beta(windMLM, method = 'sgv', partial = F)$Rsq)
(wind_betawt <- r2beta(windMLM_wt, method = 'sgv', partial = F)$Rsq)

#Getting the AIC difference from glmer, but using glmmPQL to get R2 because no sgv method for glmer models
grassMLM <- glmer(as.factor(grass_type) ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat,
                  family = binomial, control = glmerControl(optimizer="Nelder_Mead", boundary.tol = 1e-7))
grassMLM_wt <- glmer(as.factor(grass_type) ~ away_tz_diff + time_hrs_round:away_tz_diff+ (1|away_team), data = SSAC_dat,
                     weights = ps1_opt1, family = binomial, control = glmerControl(optimizer="Nelder_Mead", boundary.tol = 1e-7))
aic_grass <- extractAIC(grassMLM)[2]
aic_grasswt <- extractAIC(grassMLM_wt)[2]

grassPQL <- glmmPQL(as.factor(grass_type) ~ away_tz_diff + time_hrs_round:away_tz_diff, random = ~1|away_team, family = binomial, data = SSAC_dat)
grassPQL_wt <- glmmPQL(as.factor(grass_type) ~ away_tz_diff + time_hrs_round:away_tz_diff, random = ~1|away_team, 
                       family = binomial, weights = ps1_opt1,  data = SSAC_dat)
(grass_beta <- r2beta(grassPQL, method = 'sgv', partial = F)$Rsq)
(grass_betawt <- r2beta(grassPQL_wt, method = 'sgv', partial = F)$Rsq)

#Now Create plots showing effects of IPW on both AIC and R2 for covariates on effect modification
r2_df <- data.frame(var = c('Precipitation', 'Precipitation', 'Humidity', 'Humidity', "High Temperature", "High Temperature",
                             "Wind Speed", "Wind Speed", "Grass Type", "Grass Type"),
                     r_sq = c(precip_beta, precip_betawt, humid_beta, humid_betawt, temp_beta, temp_betawt,
                              wind_beta, wind_betawt, grass_beta, grass_betawt),
                    aic = c(aic_precip, aic_precipwt, aic_humid, aic_humidwt, aic_temp, aic_tempwt, 
                            aic_wind, aic_windwt, aic_grass, aic_grasswt),
                     sample = rep(c("Unadjusted", "Adjusted"), 5))


covariates_r2_plot <- ggplot(r2_df, aes(y= var, x=r_sq, color=sample)) + theme_pubclean() +
                              theme(panel.background = element_rect(fill = "white"),
                                    axis.text.x = element_text(color = "black"),
                                    axis.text.y = element_text(color = "black"),
                                    axis.title.y = element_blank(),
                                    panel.border = element_rect(fill = NA, color = "black"),
                                    plot.background = element_rect(fill = "white"),
                                    legend.background = element_rect(fill = "white"),
                                    legend.key = element_blank()) +
                              geom_jitter(width = 0.001, height = 0.09, size=2.1) + 
                              xlim(0.00, 0.105) + labs(x= bquote(R^2), color= "Sample")+
                              geom_vline(xintercept = .1) + geom_vline(xintercept = 0, alpha= 0.5) +
                              geom_vline(xintercept = .05, linetype=3) 


covariates_r2_plot

#then you can save these figures locally if you want to
#ggsave("covariates_r2_plot.png", covariates_r2_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 1.5)


##Here's what the non-causal analysis would look like it would look like without the weights
##Meaning the p-value is valid but just a non-causal relationship
tzone_l3_nowt <- gam(spread_beat_away ~ s(away_tz_diff, k=5) + ti(away_tz_diff, time_hrs_round, k=7),  
                     data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit")) #
anova(tzone_l3_nowt)
gam.check(tzone_l3_nowt)
#null logit model
tzone_l3_null_nowt <- gam(spread_beat_away ~ 1 ,
                          data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"))
#get F-statistic, showing null effect of those 'treatments'
anova(tzone_l3_null_nowt, tzone_l3_nowt, test='F')


#main causal inference model of interest
tzone_l3_causal <- gam(spread_beat_away ~ s(away_tz_diff, k=5) + ti(away_tz_diff, time_hrs_round, k=7),  
                       data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = ps1_opt1) #
anova(tzone_l3_causal)
gam.check(tzone_l3_causal)
#null logit model
tzone_l3_null <- gam(spread_beat_away ~ 1 ,
                     data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = ps1_opt1) #, weights = ps1_opt1
#get F-stat for incorrect inference, but used to contrast with RI for correct inferences
real_logit_ftest <- anova(tzone_l3_null, tzone_l3_causal, test='F')


##Because we're using 'weights', the Standard Errors are invalid (though the estimates are right)
##Therefore we can't get valid p-values. Bootstrap won't work with splines and
##sampling posterior won't work with weights
##The below code is employing randomization inference is commented out and has already been done and saved to github
away_tz <- SSAC_dat$away_tz
rand_int <- SSAC_dat$away_team


#' #Need to take the away team's usual time zone and use that to determine what their
#' #'possible' time changes could have been for RI (positivity assumption)
#' 
#' get_all_tz <- unique(SSAC_dat$away_tz)
#' eastern <- c(3,2,1,0)
#' central <- c(2,1,0,-1)
#' mtn <- c(1,0,-1,-2)
#' pac <- c(0,-1,-2,-3)
#' 
#' team_tz_randomizations <- vector(mode = 'list', length = 5000)
#' team_time_randomizations <- vector(mode = 'list', length = 5000)
#' 
#' time_teams <- data.frame(table(SSAC_dat$time_hrs_round, SSAC_dat$away_team))
# time_teams2<- time_teams %>% group_by(Var2) %>% mutate(frac = Freq/sum(Freq)) %>% arrange(Var1, .by_group = T)
# 
# 
# #This one creates the proportionally random kickoffs
# #It could be written to run faster but we only need to do it once
# for (g in 1:5000){
#   print(g)
#   tz_vec <- vector(mode = 'integer', length = length(away_tz))
#   time_vec <- vector(mode = 'integer', length = length(rand_int))
#   for (t in 1:length(away_tz)) {
#     tz_vec[t] <- ifelse(away_tz[t] == "America/New_York", sample(x= eastern, size = 1, replace = T),
#                         ifelse(away_tz[t] == "America/Kentucky/Louisville", sample(x= eastern, size = 1, replace = T),
#                                ifelse(away_tz[t] == "America/Indiana/Indianapolis", sample(x= eastern, size = 1, replace = T),
#                                       ifelse(away_tz[t] == "America/Detroit", sample(x= eastern, size = 1, replace = T),
#                                              ifelse(away_tz[t] == "America/Chicago", sample(x= central, size = 1, replace = T),
#                                                     ifelse(away_tz[t] == "America/Denver", sample(x= mtn, size = 1, replace = T),
#                                                            ifelse(away_tz[t] == "America/Phoenix", sample(x= mtn, size = 1, replace = T),
#                                                                   ifelse(away_tz[t] == "America/Boise", sample(x= mtn, size = 1, replace = T),
#                                                                          ifelse(away_tz[t] == "America/Los_Angeles", sample(x= pac, size = 1, replace = T),
#                                                                                 99999)))))))))
#     time_vec[t] <- sample(x = seq(11, 24, by=1), size = 1, replace = T,
#                           prob =time_teams2[ which (time_teams2$Var2== rand_int[t]),]$frac) #this one is for keeping proportions
#   }
#   team_tz_randomizations[[g]] <- tz_vec
#   team_time_randomizations[[g]] <- time_vec
# }
# 
# #This one creates fully random kickoff times
# #Note that it will over-write the above dataset if you run them both... so don't or change the names!
# for (g in 1:5000){
#   print(g)
#   tz_vec <- vector(mode = 'integer', length = length(away_tz))
#   time_vec <- vector(mode = 'integer', length = length(rand_int))
#   for (t in 1:length(away_tz)) {
#     tz_vec[t] <- ifelse(away_tz[t] == "America/New_York", sample(x= eastern, size = 1, replace = T),
#                         ifelse(away_tz[t] == "America/Kentucky/Louisville", sample(x= eastern, size = 1, replace = T),
#                                ifelse(away_tz[t] == "America/Indiana/Indianapolis", sample(x= eastern, size = 1, replace = T),
#                                       ifelse(away_tz[t] == "America/Detroit", sample(x= eastern, size = 1, replace = T),
#                                              ifelse(away_tz[t] == "America/Chicago", sample(x= central, size = 1, replace = T),
#                                                     ifelse(away_tz[t] == "America/Denver", sample(x= mtn, size = 1, replace = T),
#                                                            ifelse(away_tz[t] == "America/Phoenix", sample(x= mtn, size = 1, replace = T),
#                                                                   ifelse(away_tz[t] == "America/Boise", sample(x= mtn, size = 1, replace = T),
#                                                                          ifelse(away_tz[t] == "America/Los_Angeles", sample(x= pac, size = 1, replace = T),
#                                                                                 99999)))))))))
#     time_vec[t] <- sample(x = seq(11, 24, by=1), size = 1, replace = T) #for completely random times
#   }
#   team_tz_randomizations[[g]] <- tz_vec
#   team_time_randomizations[[g]] <- time_vec
# }


team_tz_randomizations <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_tz_randomizations.rds")))
team_time_randomizations <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_time_randomizations.rds")))
team_tz_randomizations_clust <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_tz_randomizations_clust.rds")))
team_time_randomizations_clust <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_time_randomizations_clust.rds")))

##While many Randomization Inference studies use the model coefficient, that does not work for splines
##What makes the most sense in our case is to use the 'drop 1' F-Statistic for Randomization Inference 
##In this case, while we're using this particular type of interaction spline with 'mgcv'
##and we're making the noise data with those interaction variables, the model-wide metric should be valid
##This slightly changes our interpretation, so we're asking if time zone change and game time impact performance,
##Not a question of those individual main effects (though we could run those analyses separately)

outcome_logit <- SSAC_dat$spread_beat_away

#F Statistic for Logit
effect_logit <- real_logit_ftest$F[[2]]

# build null distribution over a bunch of possible randomizations (takes me ~15 min to run)
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl, c('team_tz_randomizations', 'team_time_randomizations', 'outcome_logit')) 

nulls <- pblapply(1:length(team_tz_randomizations), function(i) {
  library(mgcv)
  data_model <- data.frame(temp1 = team_tz_randomizations[[i]], temp2 = team_time_randomizations[[i]], 
                           outcome_logit)
  try(model <- gam(outcome_logit ~ s(temp1, k=5) + ti(temp1, temp2, k=7), data = data_model, 
                   family = quasibinomial(link = "logit"), method = 'REML')) 
  try(null_m <- gam(outcome_logit ~ 1, data = data_model, 
                    family = quasibinomial(link = "logit"), method = 'REML'))
  try(anova(null_m, model, test='F')$F[[2]])
  
},cl = cl)
parallel::stopCluster(cl)

ri_logit <- unlist(nulls) 
# plot the null and observed effect
hist(ri_logit, breaks = 'fd')
abline(v = effect_logit)
# p-value
sum(abs(effect_logit) < abs(ri_logit))/length(ri_logit) # randomization test



# build null distribution over a bunch of possible randomizations with kickof clustered data (takes me ~15 min to run)
cl <- parallel::makePSOCKcluster(cores)
parallel::clusterExport(cl, c('team_tz_randomizations_clust', 'team_time_randomizations_clust', 'outcome_logit')) 

nulls_clust <- pblapply(1:length(team_tz_randomizations_clust), function(i) {
  library(mgcv)
  data_model <- data.frame(temp1 = team_tz_randomizations_clust[[i]], temp2 = team_time_randomizations_clust[[i]], 
                           outcome_logit)
  try(model <- gam(outcome_logit ~ s(temp1, k=5) + ti(temp1, temp2, k=7), data = data_model, 
                   family = quasibinomial(link = "logit"), method = 'REML')) 
  try(null_m <- gam(outcome_logit ~ 1, data = data_model, 
                    family = quasibinomial(link = "logit"), method = 'REML'))
  try(anova(null_m, model, test='F')$F[[2]])
  
},cl = cl)
parallel::stopCluster(cl)

ri_logit_clust <- unlist(nulls_clust) 
# plot the null and observed effect
hist(ri_logit_clust, breaks = 'fd')
abline(v = effect_logit)
# p-value
sum(abs(effect_logit) < abs(ri_logit_clust))/length(ri_logit_clust) # randomization test


#############################################################################################
####Making Figures of above analyses
#############################################################################################
#Making pretty ggplot of Logistic RI histogram


ri_logit_df <- data.frame('F_statistic' = ri_logit)
ri_logit_clust_df <- data.frame('F_statistic' = ri_logit_clust)

ri_logit_plot <- ggplot(ri_logit_df, aes(x=F_statistic)) + geom_histogram(bins = 50, fill='steelblue1', color='black') +
                          geom_vline(xintercept = effect_logit, linewidth=1.2) + 
                          geom_curve(aes(x=effect_logit, y= 310, xend= 3.5, yend= 250), curvature = -0.5, linewidth=1, 
                                     arrow = arrow(length = unit(0.1, 'inches'), ends = 'first', type = 'closed')) +
                          annotate("text", x= 4, y= 230, label= "F-Statistic from Observed Data\np = 0.133") +
                          ggtitle('Randomization Inference for Effect Mediation\nBetween Random Game Time & Hours Gained/Lost') + 
                          ylab('Count') + xlab('F-Statistic Distribution for All Potential Outcomes') + 
                          theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                                                   axis.text=element_text(size=14),
                                                   axis.title=element_text(size=14))

ri_logitclust_plot <- ggplot(ri_logit_clust_df, aes(x=F_statistic)) + geom_histogram(bins = 50, fill='steelblue1', color='black') +
                            geom_vline(xintercept = effect_logit, linewidth=1.2) + 
                            geom_curve(aes(x=effect_logit, y= 310, xend= 3.5, yend= 250), curvature = -0.5, linewidth=1, 
                                       arrow = arrow(length = unit(0.1, 'inches'), ends = 'first', type = 'closed')) +
                            annotate("text", x= 4, y= 230, label= "F-Statistic from Observed Data\np = 0.142") +
                            ggtitle('Randomization Inference for Effect Mediation Between\nProportional Game Time & Hours Gained/Lost') + 
                            ylab('Count') + xlab('F-Statistic Distribution for All Potential Outcomes') + 
                            theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                                                     axis.text=element_text(size=14),
                                                     axis.title=element_text(size=14))

#save to working directory if you want
#ggsave("ri_logit_plot.png", ri_logit_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 2.1)
#ggsave("ri_logitclust_plot.png", ri_logitclust_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 2.1)


#Now the surface plot of logistic regression point estimates saved to working directory
#png("surface_plot.png", width = 4, height = 4, units = 'in', res = 600)
par(mar = c(3, 3, 3, 1), cex.axis= 0.6)
fvisgam(tzone_l3_causal, view = c('away_tz_diff', 'time_hrs_round'), transform = plogis, hide.label = T, 
        add.color.legend = F, xlab='' , ylab='', main= NULL)
mtext(text = "Probability of Away Team Beating\nthe Vegas Spread", side = 3, line = 1, cex = 1)
mtext(text = "Effects Do NOT Exceed Random Noise", side = 3, line = .1, cex = .7)
mtext(text = "Hours Lost/Gained by Away Team in Travel", side = 1, line = 1.8, cex = .7)
mtext(text = "Kickoff Time in Eastern Time Zone", side = 2, line = 1.8, cex = .7)
#dev.off()
