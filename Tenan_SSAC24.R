##SSAC24 script

#This is designed to run on a Windows PC with 16 cores, some parallel processing code
#may not work appropriately if you are running Linux or do not have sufficient cores/memory
#These are relatively easily solved, but may require a few 'tweaks' to run without errors

#These are the libraries for the 'core analyses'
library(mgcv)
library(gratia)
library(itsadug)
library(dplyr)
library(WeightIt)
library(cobalt)
library(pbapply)

cores <- 16 #if you have a Windows machine and just have less than 16 cores, you can just change this

SSAC_dat <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/TENAN_SSAC2024_dat.rds")))


#Demonstrating the need for cluster-level propensity scores
#Below code shows the team-level correlation for the treatments.
#While they're relatively low (0.16 and 0.29) they're definitely high enough to bias an analysis,
#so we'll use multilevel model propensity scores
#On the other hand, the actual analysis, has very low ICC (0.015) so we won't do an MLM for final analysis
#Logit based on ICC method from Hox JJ, Moerbeek M, van de Schoot R (2018). 
#Multilevel Analysis: Techniques and Applications. Taylor and Francis. ISBN 9781138121362.  p. 107
library(lme4)
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
#marginal (opt1) and cluster mean (opt2) stabilized weights
#from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5157938/
num.mod.time <-gam(time_hrs_round ~ s(away_team, bs= 're'), data = SSAC_dat , method = 'REML')  # numerator model for option 1
den.mod.time <- gam(time_hrs_round ~ s(away_team, bs= 're')+  grass_type + 
                      stadium_type + s(precip, by= elements) + s(humidity_high, by=elements) + s(temp_high, by= elements) + 
                      s(wind_speed, by= elements) , data= SSAC_dat, method = 'REML') # denominator model for both options


num.p_op1 <- dnorm((SSAC_dat$time_hrs_round), mean = fitted(num.mod.time), 
                   sd = last(variance_comp(num.mod.time)$std_dev)) # numerator calculation for option 1 
num.p_op2 <- dnorm((SSAC_dat$time_hrs_round), mean = mean((SSAC_dat$time_hrs_round)), 
                   sd = sd((SSAC_dat$time_hrs_round))) # numerator calculation for option 2
den.p <- dnorm((SSAC_dat$time_hrs_round), mean = fitted(den.mod.time), sd = last(variance_comp(den.mod.time)$std_dev)) 
Tps_MLM_op1 <- num.p_op1/den.p
Tps_MLM_op2 <- num.p_op2/den.p

SSAC_dat$ps1_opt1 <- Tps_MLM_op1
SSAC_dat$ps1_opt2 <- Tps_MLM_op2

#confirming no major issues with the denominator propensity score
gam.check(den.mod.time)


#balance diagnostics from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6484705/
covs <- SSAC_dat %>% dplyr::select(grass_type, stadium_type, precip, elements, humidity_high,
                                           temp_high, wind_speed) %>% data.frame()
covs_time <- SSAC_dat %>% dplyr::select(grass_type, stadium_type, precip, elements, humidity_high,
                                                temp_high, wind_speed) %>% data.frame()

ps1_weigh <- as.weightit(SSAC_dat$ps1_opt1, treat= (SSAC_dat$time_hrs_round), estimand='ATT', covs = covs)
ps2_weigh <- as.weightit(SSAC_dat$ps1_opt2, treat= (SSAC_dat$time_hrs_round), estimand='ATT', covs = covs)

bal.tab(ps1_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)
bal.tab(ps2_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)

#Both methods balanced pretty well for time of game, now we'll see how they balanced for
#Hours gained/lost due to time zone changes
ps1_weigh <- as.weightit(SSAC_dat$ps1_opt1, treat= (SSAC_dat$away_tz_diff), estimand='ATT', covs = covs_time)
ps2_weigh <- as.weightit(SSAC_dat$ps1_opt2, treat= (SSAC_dat$away_tz_diff), estimand='ATT', covs = covs_time)
bal.tab(ps1_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)
bal.tab(ps2_weigh, stats = c("c", "m"), un = TRUE, thresholds = c(cor = .1), poly = 3)
#Marginal Stabilized weight PS type balanced the best and is used moving forward

#main model of interest
tzone_l3_causal <- gam(spread_beat_away ~ s(away_tz_diff, k=5) + ti(away_tz_diff, time_hrs_round, k=7),  
                       data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = ps1_opt1) 
anova(tzone_l3_causal)
gam.check(tzone_l3_causal)


#null logit model
tzone_l3_null <- gam(spread_beat_away ~ 1 ,
                     data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = ps1_opt1) 

#get F-stat
real_logit_ftest <- anova(tzone_l3_null, tzone_l3_causal, test='F')


#' ##Because we're using 'weights', the Standard Errors are invalid (though the estimates are right)
#' ##Therefore we can't get valid p-values. Bootstrap won't work with splines and 
#' ##sampling posterior won't work with weights
#' ##The below code is employing randomization inference that has already been done and saved
#' ###Can be un-commented if need to modify for update
#' rand_int <- SSAC_dat$away_team
#' away_tz <- SSAC_dat$away_tz
#' 
#' 
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
#' 
#' 
#' for (g in 1:5000){
#'   print(g)
#'   tz_vec <- vector(mode = 'integer', length = length(away_tz))
#'   time_vec <- vector(mode = 'integer', length = length(rand_int))
#'   for (t in 1:length(away_tz)) {
#'     tz_vec[t] <- ifelse(away_tz[t] == "America/New_York", sample(x= eastern, size = 1, replace = T), 
#'                    ifelse(away_tz[t] == "America/Kentucky/Louisville", sample(x= eastern, size = 1, replace = T),
#'                      ifelse(away_tz[t] == "America/Indiana/Indianapolis", sample(x= eastern, size = 1, replace = T), 
#'                        ifelse(away_tz[t] == "America/Detroit", sample(x= eastern, size = 1, replace = T),
#'                          ifelse(away_tz[t] == "America/Chicago", sample(x= central, size = 1, replace = T), 
#'                            ifelse(away_tz[t] == "America/Denver", sample(x= mtn, size = 1, replace = T),
#'                              ifelse(away_tz[t] == "America/Phoenix", sample(x= mtn, size = 1, replace = T),
#'                                ifelse(away_tz[t] == "America/Boise", sample(x= mtn, size = 1, replace = T),
#'                                  ifelse(away_tz[t] == "America/Los_Angeles", sample(x= pac, size = 1, replace = T),
#'                                  99999)))))))))
#'     
#'     
#'     time_vec[t] <- sample(x = seq(11, 24, by=1), size = 1, replace = T)
#'   }
#'   team_tz_randomizations[[g]] <- tz_vec
#'   team_time_randomizations[[g]] <- time_vec
#'   
#' }
#' 


team_tz_randomizations <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_tz_randomizations.rds")))
team_time_randomizations <- readRDS(gzcon(url("https://github.com/TenanATC/SSAC_2024/raw/main/team_time_randomizations.rds")))




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

#############################################################################################
####Making Figures of above analyses
#############################################################################################

#Making pretty ggplot of Logistic RI histogram
library(ggplot2)
library(ggpubr)

ri_logit_df <- data.frame('F_statistic' = ri_logit)

ggplot(ri_logit_df, aes(x=F_statistic)) + geom_histogram(bins = 50, fill='steelblue1', color='black') +
  geom_vline(xintercept = effect_logit, linewidth=1.2) + 
  geom_curve(aes(x=effect_logit, y= 310, xend= 3.5, yend= 250), curvature = -0.5, linewidth=1, 
             arrow = arrow(length = unit(0.1, 'inches'), ends = 'first', type = 'closed')) +
  annotate("text", x= 3.5, y= 235, label= "F-Statistic from Observed Data\np = 0.39") +
  ggtitle('Randomization Inference for Non-Linear Interaction\nBetween Game Time and Hours Gained/Lost') + 
  ylab('Count') + xlab('F-Statistic Distribution for All Potential Outcomes') + 
  theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                           axis.text=element_text(size=14),
                           axis.title=element_text(size=14))


#Now the surface plot of logistic regression point estimates
fvisgam(tzone_l3_causal, view = c('away_tz_diff', 'time_hrs_round'), transform = plogis, hide.label = T, 
        add.color.legend = F, xlab='' , ylab='', main= NULL)
mtext(text = "Probability of Away Team Beating\nthe Vegas Spread", side = 3, line = 1, cex = 1)
mtext(text = "Effects Do NOT Exceed Random Noise", side = 3, line = .1, cex = .7)
mtext(text = "Hours Lost/Gained by Away Team in Travel", side = 1, line = 1.8, cex = .7)
mtext(text = "Kickoff Time in Eastern Time Zone", side = 2, line = 1.8, cex = .7)

