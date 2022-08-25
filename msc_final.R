################################################################################
# Script to analyse environmental data, tick collection data and molecular 
# analysis of Borrelia burgdorferi sensu lato infection data.

# contents
# 1. wrangling environmental data
# 2. wrangling tick data
# 3. merging datasets together
# 4. analysis
#   4.1. calculation of densities
#   4.2. land-use and tick abundance
#   4.3. land-use and nymph abundance
#   4.4. abundances ad environmental parameters

# Author: Aaron Sambrook
# Project: MSc EEC Tick 
################################################################################
# To do:
# sort out the categorical variables i.e. move the light.cat variable
#                                         concatenate the veg variables
# add in the model selection and actual models
# produce figures
################################################################################

# load/download packages
if(!require(tidyverse)){install.packages("tidyverse")} ; require(tidyverse)
if(!require(ggplot2)){install.packages("ggplot2")} ; require(ggplot)
if(!require(lme4)){install.packages("lme4")} ; require(lme4)
if(!require(psych)){install.packages("psych")} ; require(psych)
if(!require(usdm)){install.packages("usdm")} ; require(usdm)
if(!require(MASS)){install.packages("MASS")} ; require(MASS)
if(!require(lubridate)){install.packages("lubridate")} ; require(lubridate)
if(!require(stringi)){install.packages("stringi")} ; require(stringi)
if(!require(verdis)){install.packages("verdis")} ; require(stringi)

# clear the environment
rm(list=ls())
# set the working directory
setwd("/Users/aaron/OneDrive - Imperial College London/University/Project/data")
# read in the environmental dataset
data_env <- read.csv("FieldDataSheet_environment_final_fixed_v2.csv", 
                     header = TRUE)
# read in the tick collection data set
data_ticks <- read.csv("FieldDataSheet_ticks_final.csv",
                     header = TRUE)


################################################################################
# 1. wrangling the environmental dataset
################################################################################

# check the summary and structure of the data
summary(data_env) # summary
str(data_env) # structure

# change the variable names to lower-case
names(data_env) <- tolower(names(data_env))

# create the land-use variable from the transect variable
data_env$landuse <- substring(data_env$transect, 1, 2)

# correct any inconsistencies in names convention of the landuse
data_env <- data_env %>% mutate(landuse = str_replace(landuse, "Wf", "WF")) %>% 
                         mutate(landuse = str_replace(landuse, "Df", "DF")) %>% 
                         mutate(landuse = str_replace(landuse, "Cf", "CF")) %>% 
                         mutate(landuse = str_replace(landuse, "Ug", "UG"))


# i want to see how long the average transect was 

# remove the dung variable because nothing came of this
data_env <- dplyr::select(data_env, -dung)

# calculate saturation deficit using the equation - i should find that 
data_env$saturation.deficit <- (1 - data_env$relative.humidity/100) * 
                                               4.9463 * 
                                 exp(0.062 * data_env$temperature)

# "OL" readings were taken when  light intensity was beyond the capabilities 
# of the equipment so will replace these with "average" daylight lux
data_env$light.intensity[data_env$light.intensity == "OL" |
                           data_env$light.intensity == "O.L"] <- "65000"

# convert all light intensity values to lux
data_env <- data_env %>% 
  mutate(light.intensity = as.numeric(light.intensity)) %>% 
  mutate(light.intensity = ifelse(light.unit == "k", 
                                   light.intensity*1000, 
                                   light.intensity)) 

# rename the sward.length variable to sward.height
data_env <- rename(data_env, sward.height = sward.length)

# create a veg_temp df to deal with categorical variable
veg_temp <- data_env %>% dplyr::select(transect, dominant.vegetation)
# substring to get transect names 
veg_temp$transect <- veg_temp$transect %>% 
                     substring(1, 4)
# create new dom.veg variable by concatenating dominant.vegetation 
veg_temp <- veg_temp %>% 
  group_by(transect) %>% 
  mutate(dom.veg = paste0(dominant.vegetation, collapse = "+")) %>% 
  dplyr::select(-dominant.vegetation)
veg_temp$transect <- gsub("\\..*","",veg_temp$transect)

# get rid of any duplicates
veg_temp <- distinct(veg_temp)

################################################################################
# 2. wrangling the tick dataset
################################################################################

# sum up each developmental stage for each transect
data_ticks <- data_ticks %>% mutate(L = L1 + L2 +L3 + L4 + L5, 
                                    N = N1 + N2 + N3 + N4 + N5, 
                                    A = A1 + A2 + A3 + A4 + A5) %>% 
                                    dplyr::select(c(Transect, L, N, A, Total))


# change the variable names to lower-case 
names(data_ticks) <- tolower(names(data_ticks))

# create a land-use variable by extracting the first 2 characters of the 
# transect variable 
data_ticks$landuse <- substring(data_ticks$transect, 1, 2)

# drop the life stage count by 15 m interval variables 
data_ticks <- data_ticks %>% 
              dplyr::select(c(transect, landuse, l, n , a, total))

# look at the summary and the structure of the dataset 
summary(data_ticks) #summary
str(data_ticks) #structure

# create a vector with the abreviations of the land-use types
landuses <- c("CF", "DF", "UG", "WF")

# create a temp file for loop
temp <- 0

# for loop to take the average of the environmental variables for each transect
for(i in landuses){
  temp <- data_env %>% filter( landuse == i) %>% #filter by land-use
          dplyr::select(c(transect, landuse, relative.humidity, 
              temperature, light.intensity,
              soil.moisture, mat.depth, 
              saturation.deficit, sward.height)) #keep these variables
  temp$transect_no <- temp$transect %>% substring(3) #extract transect no
  temp$transect_no <- gsub("\\..*", "", temp$transect_no)
  temp <- temp %>% dplyr::select(-c(transect, landuse)) %>% 
          group_by(transect_no) %>% # group the variables by transect number 
          summarize(relative.humidity = mean(relative.humidity),  
              temperature = mean(temperature),            # means per transect
              light.intensity = mean(light.intensity),
              soil.moisture = mean(soil.moisture),
              mat.depth = mean(mat.depth),
              saturation.deficit = mean(saturation.deficit),
              sward.height = mean(sward.height)) %>% 
          mutate(landuse = i) %>% 
          mutate(light.cat = case_when( # create a cat. variable for  
               between(light.intensity, 0, 5000) ~ "A", # light intensity
               between(light.intensity, 5001, 10000) ~ "B",
               between(light.intensity, 10001, 15000) ~ "C",
               between(light.intensity, 15001, 20000) ~ "D",
               between(light.intensity, 20001, 25000) ~ "E",
               between(light.intensity, 25001, 30000) ~ "F",
               between(light.intensity, 30001, 35000) ~ "G",
               between(light.intensity, 35001, 40000) ~ "H",
               TRUE ~ "I"))  # create a land-use variable
  assign(paste0("Check_", i), temp) # assign to a data frame
}


# add the datasets back together
data_env_2 <- rbind(Check_CF, Check_DF, Check_UG, Check_WF)

# remove temp files
rm(Check_CF, Check_DF, Check_UG, Check_WF, temp)


###############################################################################
# 3. merge the tick and environmental datasets together
###############################################################################
# recreate transect variable
data_env_2$transect <- paste(data_env_2$landuse, 
                             data_env_2$transect_no, 
                             sep ="")

# merge data sets to get final dataset
data_final <- merge(data_ticks, data_env_2) 
data_final <- merge(data_final, veg_temp)

# remove intermediate datasets and variables
rm(data_env, data_env_2, data_ticks, veg_temp, i, landuses)
data_final <- data_final %>% dplyr::select(-light.intensity) 


# sort out the pseudo replication 

data_final2 <- data_final
data_final2$transect_no <- data_final2$transect %>% substring(3) #extract transect no
temp <- data_final2 %>% filter(landuse == "WF") %>% 
  filter(transect_no == 1 | transect_no == 2 | transect_no == 9 | transect_no == 10)

temp$transect_no2 <- 0
temp <- temp %>% mutate(transect_no2 = ifelse(transect_no == "1" | transect_no == "9",
                                              "1",
                                              "2"))

temp <- temp %>% dplyr::select(-c(transect, landuse)) %>% 
  group_by(transect_no2) %>% # group the variables by transect number 
  summarize(l = mean(l), n = mean(n), total = mean(total),
            temperature = mean(temperature),            # means per transect
            soil.moisture = mean(soil.moisture),
            mat.depth = mean(mat.depth),
            saturation.deficit = mean(saturation.deficit),
            sward.height = mean(sward.height))

temp[,2:4] <- round(temp[,2:4], digits = 0)
temp <- temp %>% rename(transect_no = transect_no2)
temp$landuse <- "WF"

temp$transect <- paste(temp$landuse, 
                       temp$transect_no, 
                       sep ="")

data_final2 <- data_final %>% filter(transect != "WF1" & transect != "WF2" &
                                       transect != "WF9" & transect != "WF10") %>% 
  dplyr::select(1:13)

data_final3 <- merge(data_final2, temp, all.x = TRUE, all.y = TRUE)

###############################################################################
# 4. analysis 
###############################################################################
# contents
# 4.1. calculation of densities 
# 4.2. model selection 



###############################################################################
# 4.1. calculation of densities 
# create variable of tick density per transect
data_final <- data_final %>% mutate(t.density = (total/75) * 100) 

# print the average density of each land-use
t_density <- data_final %>% 
group_by(landuse) %>% 
summarise(t.density.mean = mean(t.density), t.density.sd = sd(t.density))

# summary(t_lm <- lm(t.density ~ landuse, data = data_final))

# create variable of nymph density per transect
data_final <- data_final %>% mutate(n.density = (n/75) * 100) 

# print the average density of each land-use
n_density <- data_final %>% 
group_by(landuse) %>% 
summarise(n.density.mean = mean(n.density), n.density.sd = sd(n.density))

# create variable of nymph density per transect
data_final <- data_final %>% mutate(l.density = (l/75) * 100) 


l_density <- data_final %>% 
  group_by(landuse) %>% 
  summarise(l.density.mean = mean(l.density), l.density.sd = sd(l.density))

# merge the density data frames together
densities <- merge(t_density, n_density) %>% dplyr::select(1, 2, 4)
densities <- merge(densities, l_density) %>% dplyr::select(1:4)

# pivot the density data frame longer
densities_long <- densities %>% 
  pivot_longer(cols = c("n.density.mean", "l.density.mean"), 
               names_to = "dev.stage", 
               values_to = "count") %>% 
              rename(density = count)

# pot the split bar graph
ggplot(data=densities_long, aes(x = landuse, y = density, fill = dev.stage))+
  geom_bar(stat="identity", position=position_dodge())+
  #geom_text(aes(label=dev.stage), vjust=1.6, color="white",
            #position = position_dodge(0.9), size=3.5)+
  scale_fill_brewer(palette="Paired", labels = c("Larvae", "Nymph"))+
  theme_minimal()+
  ylab("Tick Density")+
  xlab("Land Use")+
  guides(fill=guide_legend(title="Developmental 
        Stage")) 
  

###############################################################################
# 4.2. model selection for land-use question (total ticks)
# Poisson %>% quasi-Poisson %>% negative binomial
# shows that DF and WF have sig more ticks than the other two land uses

# poisson
# summary(m1 <- glm(formula = total ~ landuse, 
#                   family = "poisson", 
#                   data = data_final))
# # calculate the pseudo R^2 
# 1 - (m1$deviance / m1$null.deviance)
# # calculate the dispersion parameter 
# m1$deviance / m1$df.residual
# # chi square to discern if model is better than the null
# anova(m1, test = "Chisq")
# # overdispersed

# quasipoisson
# summary(m2 <- glm(formula = total ~ landuse, 
#                   family = "quasipoisson", 
#                   data = data_final))
# # calculate the pseudo R^2 
# 1 - (m2$deviance / m2$null.deviance)
# # calculate the dispersion parameter 
# m2$deviance / m2$df.residual
# # chi square to discern if model is better than the null
# anova(m2, test = "Chisq")
# # overdispersed

# negative binomial
summary(m3 <- glm.nb(formula = total ~ landuse, 
                     link = log, 
                     data = data_final3))

# calculate the pseudo R^2 
1 - (m3$deviance / m3$null.deviance)
# calculate the dispersion parameter 
m3$deviance / m3$df.residual
# chi square to discern if model is better than the null
anova(m3, test = "Chisq")


###############################################################################
# 4.3 model selection land-use and nymph
# poisson
# summary(m4 <- glm(formula = n ~ landuse, 
#                   family = "poisson", 
#                   data = data_final))
# # calculate the pseudo R^2 
# 1 - (m4$deviance / m4$null.deviance)
# # calculate the dispersion parameter 
# m4$deviance / m4$df.residual
# # chi squre to discern if model is better than the null
# anova(m4, test = "Chisq")
# # overdispersed

# quasipoisson
# summary(m5 <- glm(formula = n ~ landuse,
#                   family = "quasipoisson",
#                   data = data_final))
# # calculate the pseudo R^2
# 1 - (m5$deviance / m5$null.deviance)
# # calculate the dispersion parameter
# m5$deviance / m5$df.residual
# # chi square to discern if model is better than the null
# anova(m2, test = "Chisq")
# # overdispersed

# negative binomial
summary(m6 <- glm.nb(formula = n ~ landuse,
                     data = data_final3))


# calculate the pseudo R^2 
1 - (m6$deviance / m6$null.deviance)
# calculate the dispersion parameter 
m6$deviance / m6$df.residual
# chi square to discern if model is better than the null
anova(m3, test = "Chisq")

###############################################################################
# run linear models on the variables and land use

###############################################################################
# check the collinearity of the explanatory variables
pairs.panels(data_final[, c(11:15)])

#################################################################
# overall tick abundance and environmental parameters
m7 <- glm.nb(total ~ soil.moisture + mat.depth + saturation.deficit +
                        sward.height, 
                      data = data_final)
summary(m7)
# pseudo r^2
1 - (m7$deviance / m7$null.deviance)
# calculate the dispersion parameter 
m7$deviance / m7$df.residual

#################################################################
# n and environment
m8 <- glm.nb(n ~ soil.moisture + saturation.deficit, 
                 data = data_final3)
summary(m8)
# pseudo r^2
1 - (m8$deviance / m8$null.deviance)
# calculate the dispersion parameter 
m8$deviance / m8$df.residual

# use predict function to make a graph to show the relationship between
# environmental parameters and nymphal abundance
pred_cf <- data.frame(soil.moisture=seq(from=0.0, to=19.7, length=100), 
                      saturation.deficit = seq(from=8.30, to = 18.1, length = 100),
                      landuse = "CF")

pred_df <- data.frame(soil.moisture=seq(from=0.0, to=19.7, length=100), 
                      saturation.deficit = seq(from=4.24, to = 7.21, length = 100),
                      landuse = "DF")

pred_ug <- data.frame(soil.moisture=seq(from=0.0, to=19.7, length=100), 
                      saturation.deficit = seq(from=5.84, to = 21.1, length = 100),
                      landuse = "UG")

pred_wf <- data.frame(soil.moisture=seq(from=0.0, to=19.7, length=100), 
                      saturation.deficit = seq(from=5.53, to = 16.3, length = 100),
                      landuse = "WF")



sat_def_pred_cf<- predict(m8, newdata = pred_cf, type = "link", se.fit = TRUE) # the type="link" here predicted the fit and se on the log-linear scale. 
sat_def_pred_df<- predict(m8, newdata = pred_df, type = "link", se.fit = TRUE) # the type="link" here predicted the fit and se on the log-linear scale. 
sat_def_pred_ug<- predict(m8, newdata = pred_ug, type = "link", se.fit = TRUE) # the type="link" here predicted the fit and se on the log-linear scale. 
sat_def_pred_wf<- predict(m8, newdata = pred_wf, type = "link", se.fit = TRUE) # the type="link" here predicted the fit and se on the log-linear scale. 

# CF
pred_cf$pred<- sat_def_pred_cf$fit
pred_cf$se<- sat_def_pred_cf$se.fit
pred_cf$upperCI<- pred_cf$pred+(pred_cf$se*1.96)
pred_cf$lowerCI<- pred_cf$pred-(pred_cf$se*1.96)

# DF
pred_df$pred<- sat_def_pred_df$fit
pred_df$se<- sat_def_pred_df$se.fit
pred_df$upperCI<- pred_df$pred+(pred_df$se*1.96)
pred_df$lowerCI<- pred_df$pred-(pred_df$se*1.96)

# UG
pred_ug$pred<- sat_def_pred_ug$fit
pred_ug$se<- sat_def_pred_ug$se.fit
pred_ug$upperCI<- pred_ug$pred+(pred_ug$se*1.96)
pred_ug$lowerCI<- pred_ug$pred-(pred_ug$se*1.96)

# WF
pred_wf$pred<- sat_def_pred_wf$fit
pred_wf$se<- sat_def_pred_wf$se.fit
pred_wf$upperCI<- pred_wf$pred+(pred_wf$se*1.96)
pred_wf$lowerCI<- pred_wf$pred-(pred_wf$se*1.96)

complete<- rbind(pred_cf, pred_df, pred_ug, pred_wf)
# Making the Plot 
ggplot(data = complete, aes(x=soil.moisture, y=exp(pred)))+ 
  geom_line(aes(color=landuse, linetype = landuse), size = 1.1)+
  #geom_ribbon(aes(ymin=exp(lowerCI), ymax=exp(upperCI), fill=factor(landuse)), alpha=0.05, show.legend = FALSE)+ 
  geom_point(data_final, mapping = aes(x=soil.moisture, y=n, color=landuse))+
  labs(y="Nymphal Abundance", x="Soil Moisture (%)")+
  theme_classic()+
  #scale_colour_brewer(palette="Set2", name = "Average 
  #Saturation 
  #Deficit") +
  scale_color_viridis(discrete = TRUE, option = "D")+
  labs(linetype = 
         "Mean 
Saturation 
Deficit", colour = 
         "Mean 
Saturation 
Deficit")

# linear model to test for differences between land uses
summary(lm_soil <- lm(soil.moisture ~ landuse, data = data_final3))
summary(lm_sat <- lm(saturation.deficit ~ landuse, data = data_final3))







