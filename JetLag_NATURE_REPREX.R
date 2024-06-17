#This is the FINAL ANALYSIS of NCAA Football Jet Lag paper, cleaned for
#Nature Communications submission



library(mgcv)
library(itsadug)
library(dplyr)
library(WeightIt)
library(cobalt)
library(pbapply)
library(lme4)
library(ggplot2)
library(ggpubr)
library(stringr)
library(stringdist)
library(independenceWeights)
library(patchwork)
library(wesanderson)
library(gratia)


#Set your working directory for all of the files, this is my pathway, but not yours
setwd("C:/RStats Files/ProFootballFocus/Nature_REPREX")


cores <- 8 #if you have a Windows machine and have less than 16 cores, you can just change this

SSAC_dat <- readRDS('SSAC_dat.rds')


#zero-ing out weather variables when elements are '0' for precip and wind speed
#public reports state that indoor temp is typically set to 70 degrees, so we'll use that for temperature
# https://www.cbsnews.com/minnesota/news/how-do-they-heat-u-s-bank-stadium/
# for the relative humidity, we'll set it at a random sample 50-60. This seems like a practical choice but not one
# where I have been able to find a good evidence base
SSAC_dat$precip <- ifelse(SSAC_dat$elements == 1, SSAC_dat$precip, 0)
SSAC_dat$humidity_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$humidity_high, sample(x= seq(45,65, by=0.1), 1))
SSAC_dat$temp_high <- ifelse(SSAC_dat$elements == 1, SSAC_dat$temp_high, sample(x= seq(70, 75, by = 0.01), 1))
SSAC_dat$wind_speed <- ifelse(SSAC_dat$elements == 1, SSAC_dat$wind_speed, 0)

SSAC_dat <- data.frame(SSAC_dat)
SSAC_dat$away_team <- as.factor(SSAC_dat$away_team)


#Demonstrating the need for accounting of team-effects in propensity scores
#Below code shows the team-level correlation for the treatments.
#While they're relatively low (0.16 and 0.29) they're definitely high enough to bias an analysis,
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

#non-vegas outcome ICCs, potentially signalling issues with the 
#Roy & Forest analysis we replicate
glmer_wins <- glmer(away_win ~ 1 + (1|away_team), data = SSAC_dat, family = binomial)
summ_glmer_wins <- summary(glmer_wins)
summ_glmer_wins$varcor$away_team[[1]] / (summ_glmer_wins$varcor$away_team[[1]] + (pi^2 / 3))

lmer_varcor_score <- VarCorr(lmer(postgame_spread ~ 1 + (1|away_team), data = SSAC_dat))
score_lmer <- as.data.frame(print(lmer_varcor_score, comp=c('Variance')))
score_lmer$vcov[1] / (score_lmer$vcov[1] + score_lmer$vcov[2])


# ###################Energy Balancing Weights#####################################
covariates <- model.matrix( ~ 0 + precip +
                              humidity_high +
                              temp_high +
                              wind_speed +
                              grass_type + away_team,
                            data = SSAC_dat)
treatment <- model.matrix( ~ 0 + away_tz_diff + time_hrs_round:away_tz_diff,
                           data = SSAC_dat)
# set.seed(112233, kind = "L'Ecuyer-CMRG")
# (dcows <- independenceWeights::independence_weights(treatment, covariates))
# dcow_weights <- dcows$weights
# ggpubr::gghistogram(dcows$weights, bins = 40)
# 
# saveRDS(dcows, "dcows.rds")

#We save and use the saved weights for reproducibility reasons,
# 'independence weights' package is not currently fully reproducible due to slight stochastic in the optimizer,
#the package author states he knows the solution and is working to implement it.
#our note is that the results are not meaningfully different between runs, but fully reproducibility is ideal.
dcows <- readRDS('dcows.rds')
dcow_weights <- dcows$weights

#manually assessing the weights with the data
SSAC_dat$dcow_weights <- dcow_weights
percentiles <- quantile(dcow_weights, probs = c(0.05, 0.95))
SSAC_dat_low <- SSAC_dat %>% filter(dcow_weights <= percentiles[1])
prop.table(table(SSAC_dat_low$away_tz))
prop.table(table(SSAC_dat$away_tz))
prop.table(table(SSAC_dat_low$home_tz))
prop.table(table(SSAC_dat$home_tz))
prop.table(table(SSAC_dat_low$away_team))
x <-  data.frame(table(SSAC_dat_low$home_team.x))
x2 <- data.frame(table(SSAC_dat_low$away_team))

SSAC_dat_high <- SSAC_dat %>% filter(dcow_weights >= percentiles[2])
prop.table(table(SSAC_dat_high$away_tz))
prop.table(table(SSAC_dat$away_tz))
prop.table(table(SSAC_dat_high$home_tz))
prop.table(table(SSAC_dat$home_tz))
x <- data.frame(table(SSAC_dat_high$home_team.x))
x2 <-data.frame(table(SSAC_dat_high$away_team))

x_overall <- data.frame(table(SSAC_dat$home_team.x))
x2_overall <- data.frame(table(SSAC_dat$away_team))
#I don't see any evidence that any teams are being completely over-weighted or under-weighted 

#assessing these weights with cobalt's bal.tab, the 'standard' way of assessing covariate balance
covariates_noteams <- as.data.frame(covariates) %>% select(!starts_with('away_team'))
dcows_weighit_tz <- as.weightit(dcow_weights, covs=as.data.frame(covariates), treat= treatment[,1])
dcows_weighit_tzhrs <- as.weightit(dcow_weights, covs=as.data.frame(covariates), treat= treatment[,2])
dcows_weighit_tz_noteams <- as.weightit(dcow_weights, covs=covariates_noteams, treat= treatment[,1])
dcows_weighit_tzhrs_noteams <- as.weightit(dcow_weights, covs=covariates_noteams, treat= treatment[,2])
summ_weightit_tz <- summary(dcows_weighit_tz)
summary(dcows_weighit_tzhrs)
plot(summ_weightit_tz)
summ_weightit_tzhrs <- summary(dcows_weighit_tzhrs)
plot(summ_weightit_tzhrs)

#evaulating balance manually for reporting
tz_baltab <- as.data.frame(bal.tab(dcows_weighit_tz, stats= c('m', 'cor'), thresholds = c(m=0.1, cor=0.1), un=T)$Balance)
tzhrs_baltab <- as.data.frame(bal.tab(dcows_weighit_tzhrs, stats= c('m', 'cor'), thresholds = c(m=0.1, cor=0.1), un=T)$Balance)

#making love plots for publication
v <- data.frame(old = c("precip", "humidity_high", "temp_high", "wind_speed", 
                        "grass_typeAT", "grass_typeFT", "grass_typeReal"),
new = c("Precipitation", "Humidity", "High Temperature", 
        "Wind Speed", "Artificial Turf Stadium", "Field Turf Stadium", "Real Grass Stadium"))
`-.gg` <- function(plot, layer) {
  if (missing(layer)) {
    stop("Cannot use `-.gg()` with a single argument. Did you accidentally put - on a new line?")
  }
  if (!is.ggplot(plot)) {
    stop('Need a plot on the left side')
  }
  plot$layers = c(layer, plot$layers)
  plot
}
(loveplot_plot1 <- love.plot(dcows_weighit_tz_noteams, stats= c('cor'),
                             var.order = 'unadjusted', 
                             size = 2.5,
                             var.names = v, 
                             colors = wes_palette("FantasticFox1",5,"discrete")[3:4], 
                             title = 'Time Zone',
                             themes = theme_classic()) +
    xlim(-.15,.15) -
    geom_rect(xmin = -0.1,
              xmax = 0.1,
              ymin = -1,
              ymax = 8,
              fill = "lightgrey",
              alpha = 0.05) +
    theme(axis.title.x = element_blank()))

(loveplot_plot2 <- love.plot(dcows_weighit_tzhrs_noteams, 
                             stats= c('cor'),
                             var.order = 'unadjusted',
                             size = 2.5,
                             var.names = v, 
                             colors = wes_palette("FantasticFox1",5,"discrete")[3:4],
                             title = 'Time Zone × Kickoff Time',
                             themes = theme_classic()) + 
    xlim(-.15,.15) -
    geom_rect(xmin = -0.1,
              xmax = 0.1,
              ymin = -1,
              ymax = 8,
              fill = "lightgrey",
              alpha = 0.05))
(loveplot_comb <- loveplot_plot1 / loveplot_plot2 + plot_layout(guides = "collect", heights = c(1,1,0.005)) & 
    theme(legend.position = "bottom", 
          legend.title = element_blank()))

# ggsave("loveplot.png", loveplot_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 1.5)



#main causal inference model of interest 
tzone_l3_causal <- gam(spread_beat_away ~ s(away_tz_diff, k=5) + ti(away_tz_diff, time_hrs_round, k=7),  
                       data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = dcow_weights) #
anova(tzone_l3_causal)
gam.check(tzone_l3_causal)
#null logit model
tzone_l3_null <- gam(spread_beat_away ~ 1 ,
                     data = SSAC_dat, method = 'REML', family = quasibinomial(link = "logit"), weights = dcow_weights) #, weights = dcow_weights
#get F-stat for incorrect inference, but used to contrast with RI for correct inferences
real_logit_ftest <- anova(tzone_l3_null, tzone_l3_causal, test='F')


##Because we're using 'weights', the Standard Errors are invalid (though the estimates are right)
##Therefore we can't get valid p-values. Bootstrap is not pointwise valid with penalized splines and
##sampling posterior won't work with weights
##The below code employing randomization inference is commented out and has already been saved to be reproducible
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


#reading in the RI datasets created with above code
team_tz_randomizations_clust <- readRDS("team_tz_randomizations_clust.rds")
team_time_randomizations_clust <- readRDS("team_time_randomizations_clust.rds")



##While many Randomization Inference studies use the model coefficient, that does not work for splines
##What makes the most sense in our case is to use the 'drop 1' F-Statistic for Randomization Inference
##In this case, while we're using this particular type of interaction spline with 'mgcv'
##and we're making the noise data with those interaction variables, the model-wide metric should be valid
##This slightly changes our interpretation, so we're asking if time zone change and game time impact performance,
##Not a question of those individual main effects (though we could run those analyses separately)

##Vegas Spread of Interest
outcome_logit <- SSAC_dat$spread_beat_away

#F Statistic for Logit
effect_logit <- real_logit_ftest$F[[2]]


# build null distribution over a bunch of possible randomizations with kickof clustered data, takes ~25 min to run on my computer
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
ri_logit_clust_df <- data.frame('F_statistic' = ri_logit_clust)

ri_logitclust_plot <- ggplot(ri_logit_clust_df, aes(x=F_statistic)) + geom_histogram(bins = 50, fill='steelblue1', color='black') +
  geom_vline(xintercept = effect_logit, linewidth=1.2) + 
  geom_curve(aes(x=effect_logit, y= 310, xend= 4.5, yend= 250), curvature = -0.5, linewidth=1, 
             arrow = arrow(length = unit(0.1, 'inches'), ends = 'first', type = 'closed')) +
  annotate("text", x= 5, y= 230, label= "F-Statistic from Observed Data\np = 0.142") +
  ggtitle('Randomization Inference for Effect Mediation on\nBeating the Las Vegas Spread') + 
  ylab('Count') + xlab('F-Statistic Distribution for All Potential Outcomes') + 
  theme_pubclean() + theme(plot.title = element_text(hjust = 0.5, size = 20),
                           axis.text=element_text(size=14),
                           axis.title=element_text(size=14))

(distr <- ggplot() + 
  geom_histogram(data = ri_logit_clust_df, 
                 mapping = aes(x=F_statistic),
                 bins = nclass.FD(ri_logit_clust_df$F_statistic), 
                 fill='darkgray') +
  geom_vline(xintercept = effect_logit,
             color = wes_palette("FantasticFox1",5,"discrete")[4], 
             linewidth = 1.5) +
  annotate("text", x=3.5, y= 230, 
           label= "Observed",
           color = wes_palette("FantasticFox1",5,"discrete")[4]) +
  coord_cartesian(expand = F) +
  xlab("F-statistic") +
  theme_classic() +
  theme(axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank())
  )

#save to working directory if you want
#ggsave("ri_logitclust_plot.png", ri_logitclust_plot, width = 3.5, height = 4, units = 'in', dpi = 500, scale = 2.1)


#Now the surface plot of logistic regression point estimates saved to working directory
# png("surface_plot_spread.png", width = 4, height = 4, units = 'in', res = 600)
# par(mar = c(3, 3, 3, 1), cex.axis= 0.6)
# fvisgam(tzone_l3_causal, view = c('away_tz_diff', 'time_hrs_round'), transform = plogis, hide.label = T, 
#         add.color.legend = F, xlab='' , ylab='', main= NULL)
# 
ds1 <- gratia::data_slice(tzone_l3_causal, away_tz_diff = evenly(away_tz_diff, 50), time_hrs_round = evenly(time_hrs_round, 50))
fv1 <- gratia::fitted_values(tzone_l3_causal, ds1, scale = "response")
pal <- wes_palette("Zissou1", 100, "continuous")
(surface_plot <-
  ggplot(data = fv1,
       aes(x = away_tz_diff,
           y = time_hrs_round,
           z = fitted,
           fill = fitted)) +
  ggrastr::geom_tile_rast() +
  metR::geom_contour2(aes(label = stat(level))) +
  coord_cartesian(xlim=c(-3,3),
                  ylim=c(11,24),
                  expand = F) +
  scale_y_continuous(breaks = seq(12,24,by=4)) +
  scale_fill_gradientn(colors = pal, lim = c(0,1), guide = guide_colorbar(theme = theme(legend.key.width = unit(0.75,"lines")))) +
  labs(fill = "Probability",
       x = "Hours lost/gained by away team travel",
       y = "Kickoff time (Eastern Time)")
)

design <- "
AB
AC
"
ci_results <-
  wrap_elements(full = loveplot_comb) + surface_plot + distr +
  plot_layout(ncol = 2, widths = c(1.2,1), heights = c(1,0.1), 
              design = design) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = 'bold'))

#ggsave("fig2.png", ci_results, width = 8, height = 4.5)




####Trying to re-create the analysis from Smith et al 1999
###"A variable x was defined based on the score of the west coast team minus the score of the east coast team 
###plus the point spread. Thus, the value of x would be positive if the west coast team
###beat the point spread and negative if the west coast team did not beat the point spread. One sample t-tests were performed to
###assess whether the x values were significantly greater than zero."

#overall
smith_west <- SSAC_dat %>% filter(away_tz_diff == -3) %>% mutate(smith_x = (away_score - home_score) - spread_away) %>%
  dplyr::select(home_team.x, away_team, home_score, away_score, smith_x, time_hrs_round)
smith_east <- SSAC_dat %>% filter(away_tz_diff == 3) %>% mutate(smith_x = (home_score - away_score) + spread_away) %>%
  dplyr::select(home_team.x, away_team, home_score, away_score, smith_x, time_hrs_round)
smith_dat <- rbind.data.frame(smith_east, smith_west)
t.test(smith_dat$smith_x)
t.test(smith_dat$smith_x)$stderr

#just evening games
smith_west_evening <- SSAC_dat %>% filter(away_tz_diff == -3) %>% mutate(smith_x = (away_score - home_score) - spread_away) %>%
  filter(time_hrs_round >=20) %>%  dplyr::select(home_team.x, away_team, home_score, away_score, smith_x)
smith_east_evening <- SSAC_dat %>% filter(away_tz_diff == 3) %>% mutate(smith_x = (home_score - away_score) + spread_away) %>%
  filter(time_hrs_round >=20) %>% dplyr::select(home_team.x, away_team, home_score, away_score, smith_x)
smith_evening_dat <- rbind.data.frame(smith_east_evening, smith_west_evening)
t.test(smith_evening_dat$smith_x)
t.test(smith_evening_dat$smith_x)$stderr
#non-significant, so we definitely can't replicate their results at all

gghistogram(smith_dat, x= "smith_x", fill = 'skyblue', bins = 20, add_density = T, 
            xlab = "Number of Points West Coast Team Did or\nDid Not Beat The Las Vegas Point Spread By",
            title = "Inability to Replicate Smith et al.'s Findings",
            ylab= "Count") + 
  geom_vline(xintercept = 0, linewidth=2) + annotate("text", x=19, y=9, label= "By: Matt Tenan PhD ATC")

gghistogram(smith_evening_dat, x= "smith_x", fill = 'skyblue', bins = 10, add_density = T, 
            xlab = "Number of Points West Coast Team Did or\nDid Not Beat The Las Vegas Point Spread By",
            title = "Jet Lag Effects Favor West Coast Teams? Not in the NCAA!") + 
  geom_vline(xintercept = 0, linewidth=2) + annotate("text", x=19, y=9, label= "By: Matt Tenan PhD ATC")


####Now trying to re-create the Roy and Forest NFL analysis
### "The outcomes for each game were: game time (afternoon/17:00 hours Eastern Standard Time or earlier;
### evening/17:30 hours Eastern Standard Time or later); the result of the game 
### (then the winning percentage was calculated for all games); the number of time zones travelled;
### and the direction of travelling (westward, same time zone or eastward)."  ...
### Need to look at their Tables 2 and 3 to replicate results, also report any likely issues.
### Also review Figure 1 for their simple linear regression, reporting any likely issues.

royforest_dat <- SSAC_dat %>% select(season.x, away_team, away_tz, away_win, time_hrs, away_tz_diff)
royforest_dat_regression <- royforest_dat %>% filter(time_hrs >17.5) %>% 
                              mutate(timezone_royforest= away_tz_diff*-1)
royforest_dat_table2 <- royforest_dat %>% mutate(direction = ifelse(away_tz_diff >0, 'Westward',
                                                                    ifelse(away_tz_diff <0, 'Eastward', 'Same')),
                                                 time = ifelse(time_hrs <17, 'Afternoon', 'Evening')) %>%
                              group_by(away_team, direction, time) %>%
                              summarise(winning_pct = 100*(sum(away_win)/n()), .groups='keep') 
royforest_dat_table3 <- royforest_dat %>% mutate(direction = ifelse(away_tz_diff >0, 'Westward',
                                                                    ifelse(away_tz_diff <0, 'Eastward', 'Same'))) %>%
                        filter(time_hrs >17.5) %>%
  group_by(away_team, direction, away_tz_diff) %>%
  summarise(winning_pct = 100*(sum(away_win)/n()), .groups='keep') 
#add data stuff for their t-tests here
table2_tests <- pairwise.t.test(royforest_dat_table2$winning_pct, royforest_dat_table2$direction,
                                p.adjust.method= 'none')
group_by(royforest_dat_table2,time,direction)%>% summarise(tbl2means = mean(winning_pct))

table3_tests <- pairwise.t.test(royforest_dat_table3$winning_pct, as.factor(royforest_dat_table3$away_tz_diff),
                                p.adjust.method= 'none')
group_by(royforest_dat_table3,away_tz_diff,direction)%>% summarise(tbl2means = mean(winning_pct))


#simple regression analysis
royforest_reg <- lm(away_win ~ timezone_royforest, data = royforest_dat_regression)
summary(royforest_reg)
cor.test(royforest_dat_regression$timezone_royforest, royforest_dat_regression$away_win)


