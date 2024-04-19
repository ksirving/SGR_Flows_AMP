## data and model exploration for stem water potential - RELATIVE CHANGE

library(tidyverse)
library(tidylog)
library(lubridate)

library(rfUtilities)
library(randomForest)

library(ks)
library(sm)
library(caret)
library(rpart)

## output file for figures
out.dir <- "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SGR_Flows_AMP/Figures/"


# upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
head(bioMeanQ_long_dist)
names(bioMeanQ_long_dist)

## check damaged trees
## 
damTrees <- bioMeanQ_long_dist %>%
  filter(Damage =="Y")

length(unique(damTrees$PlantID)) ## 18

length(unique(bioMeanQ_long_dist$PlantID)) ## 119

## change blanks to NA
## change blank values in Damage
bioMeanQ_long_dist["Damage"][bioMeanQ_long_dist["Damage"]==''] <- NA

## replace NA with N in damage
bioMeanQ_long_dist$Damage[is.na(bioMeanQ_long_dist$Damage)] <- "N"

unique(bioMeanQ_long_dist$Damage)
sum(bioMeanQ_long_dist$Damage == "Y") ## 504

## make columns wider
bioMeanQ_long_distx1 <- bioMeanQ_long_dist %>%
  group_by(PlantID) %>%
  select(PlantID:Year, SWP, Replacement, Damage) %>%
  distinct() %>%
  pivot_wider(names_from = "Year", values_from = "SWP") %>%
  select(-"NA")

## get info on damage - what year and season damage was noted
damdata <- bioMeanQ_long_dist %>%
  select(PlantID, Season, Year, Damage) %>%
  mutate(DamageYear = ifelse(Damage == "Y", as.character(Year), "None")) %>%
  mutate(DamageSeason = ifelse(Damage == "Y", as.character(Season), "None")) %>%
  distinct() %>%
  drop_na(DamageYear) %>% filter(Damage == "Y")

### join together

# data <- left_join(bioMeanQ_long_distx, damdata, by = c("PlantID", "Season", "Year"))


# Relative and delta change -----------------------------------------------

## remove columns and make distinct
bioMeanQ_long_distx <- bioMeanQ_long_dist %>%
  group_by(PlantID) %>%
  select(PlantID:Year, SWP , Replacement, Damage) %>%
  distinct()
bioMeanQ_long_distx

## if damage then plantid for remaining surveys should be also damaged - for lag vlaues calculation
# ## function to create a sequence, to have re[eated damged trees from first occurence of damage
# create_sequence <- function(x) {
#   first_y <- which.max(x == "Y")
#   ifelse(seq_along(x) <= first_y, "N", "Y")
#   return(sequence)
# }

## add the column with repeated damage
damSeq <- bioMeanQ_long_distx %>%
  drop_na(SWP) %>%
  group_by(PlantID) %>%
  mutate(firstY = which.max(Damage == "Y")) %>% ## index of first damage
  mutate(firstY = ifelse(firstY == 1, 0, firstY)) %>% ## change 1 to 0 
  mutate(Damage2 = ifelse(seq_along(Damage) >= firstY | Damage == "Y", "Y", "N")) %>% ## sequence from first Y to repeat until next group
  mutate(Damage2 = ifelse(firstY == 0, "N", Damage2)) ## change all values with first y as 0 from N to Y (no damage in these groups)
damSeq

### calculate relative change per year
bioMeanQ_long_distRel <- damSeq %>%
  drop_na(SWP) %>%
  group_by(PlantID, Season) %>%
  mutate(lag.SWP = dplyr::lag(SWP, n = 1, default = NA)) %>% ## lag SWP one season
  mutate(SWP.relative = SWP/lag.SWP) %>% ## relative SWP
  mutate(SWP.difference = SWP-lag.SWP) %>%## delta SWP i.e., difference between
  ungroup() %>%
  group_by(PlantID, Damage2, Season) %>%
  mutate(lag.SWPrs = dplyr::lag(SWP, n = 1, default = NA)) %>% ## lag SWP one season, restart from damage observed
  mutate(SWP.relativeRS = SWP/lag.SWPrs) %>% ## relative SWP
  mutate(SWP.differenceRS = SWP-lag.SWPrs) %>%## delta SWP i.e., difference between
  ungroup()

## join back to main df

bioMeanQ_long_distJoin <- full_join(bioMeanQ_long_distRel, bioMeanQ_long_dist, by = c("PlantID","Species", "Season", "Group", "LifeForm", "Year", "SWP", "Replacement", "Damage"))

# SWP model with damaged tree variable ---------------------------------------------------------------

## try relative & difference, with and without restarting
## rearrange for ease into RF - remove observed rainfall - 
names(bioMeanQ_long_distJoin)
rf.data <- bioMeanQ_long_distJoin %>%
  dplyr::select(SWP.relative, SWP.difference, SWP.relativeRS, SWP.differenceRS, SWP, PlantID, Species:LifeForm, Month, Year, Season, POS, Substrate,
                RainFallMod:TempMeanModF, Replacement, Damage2, SJC002_POM001Combined:WN002, ## include damaged trees in all years after damage occurred
                Q, DistToSJC002, Source, lag.SWP, lag.SWPrs) %>%
  mutate(POS = as.factor(POS),
         Species = as.factor(Species),
         Group = as.factor(Group),
         LifeForm = as.factor(LifeForm),
         Month = as.factor(Month),
         Substrate = as.factor(Substrate),
         Damage2 = as.factor(Damage2)) %>%
  # mutate(POS = ifelse(Year == 2018, NA, POS)) %>%
  # mutate(POS = recode_factor(POS, blanks = NA)) %>%
  filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  # drop_na(SWP.relative) %>%
  mutate(SWP.relative = ifelse(is.na(lag.SWP), 1, SWP.relative)) %>%
  mutate(SWP.relativeRS = ifelse(is.na(lag.SWPrs), 1, SWP.relativeRS)) %>%
  mutate(SWP.difference = ifelse(is.na(lag.SWP), 1, SWP.difference)) %>%
  mutate(SWP.differenceRS = ifelse(is.na(lag.SWPrs), 1, SWP.differenceRS)) %>%
  select(-Source, -lag.SWP, - lag.SWPrs)

names(rf.data)
dim(rf.data)

## define plantids 
plantids <- as.data.frame(rf.data$PlantID) %>%
  rename(PlantID = 1)

plantids$SWP <- rf.data$SWP
plantids
## remove plantid from df

rf.data <- rf.data %>% select(-PlantID, -SWP)

# ## change blank values in POS
# rf.data["POS"][rf.data["POS"]==''] <- NA
# 
# rf.data["POS"][rf.data["Year"]=='2018'] <- NA

## make consistant POS values - change LM to L and 

# Random Forest Model -----------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## impute missing values -for each response to test

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")


# Relative  -------------------------------------------------------

## relative per tree
rf.dataR <- rf.data %>%
  select(-c(SWP.relativeRS,SWP.difference, SWP.differenceRS))

rf.data.imputedRel <- rfImpute(SWP.relative ~ ., rf.dataR)
head(rf.data.imputedRel)

## make y data compatible
rf.data.val <- rf.data.imputedRel %>%
  rename(y = SWP.relative) %>% as.data.frame() %>%
  droplevels()

## run the model
rf1 <- randomForest(y~., data=rf.data.val, importance = T)
rf1
mean(rf1$rsq) ##  rel per tree - 76%

## importance quick look
varImpPlot(rf1) ## season
importance(rf1, type = 1) 

## relative per tree and damage to tree
rf.dataRRS <- rf.data %>%
  select(-c(SWP.relative, SWP.difference, SWP.differenceRS))

rf.data.imputedRelRS <- rfImpute(SWP.relativeRS ~ ., rf.dataRRS)
head(rf.data.imputedRelRS)

## make y data compatible
rf.data.valRS <- rf.data.imputedRelRS %>%
  rename(y = SWP.relativeRS) %>% as.data.frame() %>%
  droplevels()

## run the model
rf2 <- randomForest(y~., data=rf.data.valRS, importance = T)
rf2
mean(rf2$rsq) ## re per tree restart - 75%

## importance quick look
varImpPlot(rf2) ## season
importance(rf2, type = 1)

# Delta -------------------------------------------------------------------

## difference per tree
rf.dataD <- rf.data %>%
  select(-c(SWP.relativeRS, SWP.relative, SWP.differenceRS))

rf.data.imputedDiff <- rfImpute(SWP.difference ~ ., rf.dataD)
head(rf.data.imputedRel)

## make y data compatible
rf.data.valD <- rf.data.imputedDiff %>%
  rename(y = SWP.difference) %>% as.data.frame() %>%
  droplevels()

## run the model
rf3 <- randomForest(y~., data=rf.data.valD, importance = T)
rf3
mean(rf3$rsq) ##  diff per tree - 75%

## importance quick look
varImpPlot(rf3) ## season
importance(rf3, type = 1)

## difference per tree and damage to tree
rf.dataDRS <- rf.data %>%
  select(-c(SWP.relativeRS,SWP.difference, SWP.relative))

rf.data.imputedDiffRS <- rfImpute(SWP.differenceRS ~ ., rf.dataDRS)
head(rf.data.imputedRel)

## make y data compatible
rf.data.valDRS <- rf.data.imputedDiffRS %>%
  rename(y = SWP.differenceRS) %>% as.data.frame() %>%
  droplevels()

## run the model
rf4 <- randomForest(y~., data=rf.data.valDRS, importance = T)
rf4
mean(rf4$rsq) ## diff per tree restart - 71%

## importance quick look
varImpPlot(rf4) ## season
importance(rf4, type = 1)


# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/01_var_imp_full_model_SWP_Relative.pdf"), width=12, height=8)

varImpPlot(rf2, type = 1) ## 

dev.off()

# Remove negative importance variables ------------------------------------
## extract importance
imp <- importance(rf1, type = 1)
imp <- as.data.frame(imp)
## get all negative importance vars
vars <- ifelse(imp$`%IncMSE` > 0, rownames(imp), NA)
## remove neg importance vars
vars <- vars[!is.na(vars)]
vars
names(rf.data)
## WN001, WN002, POM001, RainFallMod, TempMeanModF, 

rf.data.red <- rf.data.imputedRel %>%
  select(SWP.relative, all_of(vars), -LifeForm) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

rf.data.val <- rf.data.red %>%
  rename(y = SWP.relative) %>% as.data.frame()

str(rf.data.red)

## run reduced model
rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ##  77%
rf
save(rf, file = "models/04_SWP_final_rf.RData")

## check without damaged trees

## make y data compatible
# rf.data.valDam <- rf.data.val %>%
#   filter(Damage2 == "N") 
# 
# 
# ## run the model
# rf2b <- randomForest(y~., data=rf.data.valDam, importance = T)
# rf2b
# mean(rf2$rsq) ##  75%


# Variable importance -----------------------------------------------------

library("scales")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
## get mean and scale and % to make relative

VarImp <- as.data.frame(rf$importance) 
colnames(VarImp)[1] <- "DecreaseAcc"
VarImp$Variable <- rownames(VarImp)


VarImp

ImpMean <- VarImp %>% 
  # group_by(Variable) %>%
  # summarise(MeanImp = mean(MeanDecreaseAccuracy)) %>%
  mutate(MeanImpScaled = rescale(DecreaseAcc)) %>%
  mutate(MeanImpPerc = (MeanImpScaled/sum(MeanImpScaled))*100)

ImpMean

## humanise variables - check the remote senseing
ImpMean <- ImpMean %>%
  mutate(VariableHuman = case_when(Variable == "POS" ~ "Tree Position",
                                   Variable == "SJC002_POM001Combined" ~ "Discharge: SJC & POM",
                                   Variable == "DischargeSJC002_POM001_WN001Combined" ~ "Discharge: SJC, POM & WN",
                                   Variable == "SJC_002" ~ "Discharge: SJC",
                                   Variable == "Replacement" ~ "Tree Replaced",
                                   Variable == "DistToSJC002" ~ "Distance to SJC",
                                   Variable == "Species" ~ "Species",
                                   Variable == "Group" ~ "Group",
                                   Variable == "Year" ~ "Year",
                                   Variable == "Season" ~ "Season",
                                   Variable == "Substrate" ~ "Substrate",
                                   Variable == "Q" ~ "Stream Flow",
                                   Variable == "WN002" ~ "Discharge: WN2",
                                   Variable == "WN001" ~ "Discharge: WN1",
                                   Variable == "POM001" ~ "Discharge: POM",
                                   Variable == "Damage2" ~ "Damage to Tree",
                                   Variable == "RainfallIntensity" ~ "Rainfall Intensity",
                                   Variable == "TempMeanModF" ~ "Air Temp (°F)"
  )) %>%
  filter(!Variable =="RainFallMod")

## order in increasing values

ImpMean$VariableHuman <- factor(ImpMean$VariableHuman, levels=ImpMean[order(ImpMean$MeanImpPerc,decreasing=F),]$VariableHuman)
ImpMean$VariableHuman

## plot
i1 <- ggplot(ImpMean, aes(x=MeanImpPerc, y=VariableHuman)) +
  geom_point() +
  scale_x_continuous("Relative Importance (%)") +
  scale_y_discrete("") +
  theme_bw()
i1

file.name1 <- "Figures/04_relative_importance_SWP_relative_inc_damage.jpg"
ggsave(i1, filename=file.name1, dpi=300, height=5, width=8)


# Partial Plots -----------------------------------------------------------

str(rf.data.val)

jpeg(paste0( "Figures/04_SWP_relative_Season.jpg"))
partialPlot(rf, rf.data.val, Season)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Distance.jpg"))
partialPlot(rf, rf.data.val, DistToSJC002)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Q.jpg"))
partialPlot(rf, rf.data.val, Q)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Species.jpg"))
partialPlot(rf, rf.data.val, Species)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_SJC.jpg"))
partialPlot(rf, rf.data.val, SJC_002)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Replacement.jpg"))
partialPlot(rf, rf.data.val, Replacement)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Group.jpg"))
partialPlot(rf, rf.data.val, Group)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_POS.jpg"))
partialPlot(rf, rf.data.val, POS)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Substrate.jpg"))
partialPlot(rf, rf.data.val, Substrate)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Year.jpg"))
partialPlot(rf, rf.data.val, Year)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Damage.jpg"))
partialPlot(rf, rf.data.val, Damage2)
dev.off()


# Predictions -------------------------------------------------------------

## predict change in SWP in response to reduction in discharge

## run model with only one of the outfall variables
names(rf.data.red)
rf.data.val <- rf.data.red %>%
  rename(y = SWP.relative) %>% as.data.frame() %>%
  select(-c(DischargeSJC002_POM001_WN001Combined:WN002)) 

range(rf.data.val$SJC002_POM001Combined)
## join SWP to run model on real values

rf.data.valRaw <- cbind(plantids, rf.data.val) %>%
  select(-y, - PlantID) %>%
  drop_na(SWP)
head(rf.data.valRaw)

dim(rf.data.val)
## run reduced model
rfOF <- randomForest(SWP~., data=rf.data.valRaw, importance = T, ntree = 1000)
mean(rfOF$rsq) ## 90% accuracy
rfOF

## plot swp against combined Q
names(rf.data.valRaw)
p0 <- ggplot(data = rf.data.valRaw, aes(x= SJC002_POM001Combined, y = SWP, group=Group, color = Group)) +
  stat_smooth() +
  facet_wrap(~Season)

p0

file.name1 <- "Figures/04_SWP_comb_Q_group.jpg"
ggsave(p0, filename=file.name1, dpi=300, height=5, width=8)

p0 <- ggplot(data = rf.data.valRaw, aes(x= SJC002_POM001Combined, y = SWP, group=Species, color = Species)) +
  stat_smooth() +
  facet_wrap(~Season)

p0

file.name1 <- "Figures/04_SWP_comb_Q_species.jpg"
ggsave(p0, filename=file.name1, dpi=300, height=5, width=8)

c = 2
head(rf.data.val)

## define reductions
percentChange <- c(rev(seq(0,1,0.01)))
percentChange

## make empty df

pred_dfx <- NULL

## loop pver predcitions
for(c in 1:length(percentChange)) {
  
  ## make new data with reduced discharge
  rf.data.valNew <- rf.data.val %>%
    mutate(SJC002_POM001Combined = SJC002_POM001Combined*percentChange[c]) %>%
    select(-y)
  
  ## predict on new data and rename response
  predSWP <- as.data.frame(predict(rfOF, newdata = rf.data.valNew)) %>%
    rename(SWP = 1) 
  
  head(predSWP)

  ## add new data to predicted values
  pred_df <- cbind(predSWP, rf.data.valNew, rf.data[, "Month"])  %>%
    mutate(PercentReduction = 1-percentChange[c]) %>%
    mutate(PlantID = plantids$PlantID)
  
  pred_dfx <- rbind(pred_dfx, pred_df)
  
}

head(pred_df)

## get mean SWP per plant id 

pred_dfx1 <- pred_dfx %>%
  group_by(PlantID, Season, PercentReduction) %>%
  mutate(SWPMean = median(SWP)) %>%
  mutate(PercentReduction = as.factor(PercentReduction))

## save out 

save(pred_dfx1, file = "ignore/04_SWP_reduction_predictions.RData")



# Plot --------------------------------------------------------------------

## reduce df

pred_red <- pred_dfx1 %>%
  select(SWPMean,  PlantID, Season, Species, Group, PercentReduction,  Year) %>%
  filter(Year %in% c(2020, 2021, 2022)) %>%
  distinct()

pred_red

test <- pred_red %>%
  filter(PlantID == "AW-1-1")

p1 <- ggplot(data = pred_red, aes(x = PercentReduction, y = SWPMean, group = Species, color = Species)) +
  stat_smooth() +
  facet_wrap(~Group)

p1

file.name1 <- "Figures/04_predictions_reduction_outfallQ.jpg"
ggsave(p1, filename=file.name1, dpi=300, height=5, width=8)


# Significance tests on reduced discharge ---------------------------------

## compare baseline response to reduced discharge response to test for significant reduction

pred_dfx2 <- pred_dfx %>% 
  mutate(PercentReductionScenario = paste0("Scenario",round(PercentReduction*100, digits = 0)))
pred_dfx2
## define reductions
percentChangeScen <- unique(pred_dfx2$PercentReductionScenario)
percentChangeScen

## get baseline predictions (after reduction)

baseline <- pred_dfx2 %>%
  filter(PercentReductionScenario == percentChangeScen[1], Year %in% c(2018, 2019, 2020, 2021, 2022)) %>%
  select(SWP, PlantID, Year, Season, Species)

## get SWP after 2020 reduction
pred <- pred_dfx2 %>%
  filter(PercentReductionScenario == percentChangeScen[20], Year %in% c(2018, 2019, 2020, 2021, 2022)) %>%
  select(SWP, PlantID, Year, Season) %>%
  rename(PredSWP = SWP)

dat <- cbind(baseline, pred[,"PredSWP"]) %>%
  rename(PredSWP = 6)

dat

datSP <- dat %>%
  filter(Species == "MF")

ttest <- t.test(datSP$SWP, datSP$PredSWP, paired =T)
## no significance
ttest

## empty dataframe
df <- data.frame(matrix(ncol = 9))
colnames(df) <- c("Reduction", "T_Statistic", "PValue", "DegreesOfFreedom", "nBaseline", "nFuture", "meanBaseline", "meanPredicted",  "Scenario")

## empty x df
# dfx <- NULL
# t=100
t=1
  ## loop over reduction scenarios
 for(t in 1:length(percentChangeScen)) {
   
   ## get predictions for specific reduction - only for years after reduction started
   
   pred <- pred_dfx2 %>%
     filter(PercentReductionScenario == percentChangeScen[t], Year %in% c(2018, 2019, 2020, 2021, 2022)) %>%
     select(SWP, PlantID, Year, Season, PercentReduction) %>%
     rename(PredSWP = SWP)
   
   ## join with baseline
   dat <- cbind(baseline, pred[,"PredSWP"]) %>%
     rename(PredSWP = 5)

   ## paired t test
   ttest <- t.test(dat$SWP, dat$PredSWP, paired =T)

   ## get stats
   df[t,1] <- pred$PercentReduction[1]
   df[t,2] <- ttest$statistic
   df[t,3] <- ttest$p.value
   df[t,4] <- ttest$parameter
   df[t,5] <- length(na.omit(dat$SWP))
   df[t,6] <- length(na.omit(dat$PredSWP))
   df[t,7] <- mean(na.omit(dat$SWP))
   df[t,8] <- mean(na.omit(dat$PredSWP))
   df[t,9] <- percentChangeScen[t]
   
   
 }

head(df)
df <- df %>%
  mutate(PValueR = round(PValue, digits = 3))
 
plot(df$Reduction, df$PValueR)

## 31% reduction is where SWP changes significantly
## baseline - 85.23819
## predicted - 85.21370

## discharge at 31%

t=31
names(pred_dfx2)
pred <- pred_dfx2 %>%
  filter(PercentReductionScenario == percentChangeScen[t], Year %in% c(2020, 2021, 2022)) %>%
  select(SWP, PlantID, Year, Season, PercentReduction, SJC002_POM001Combined) %>%
  rename(PredSWP = SWP)

range(pred$SJC002_POM001Combined) ## 0.0000 21.1829


# Cross validation --------------------------------------------------------

## split into training and testing

setsize <- floor(nrow(rf.data.red)*0.8)
index <- sample(1:nrow(rf.data.red), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(rf.data.val)

(rf.SWP <- rf.crossValidation(rf, rf.data.val[, 2:13], 
                             p=0.10, n=99, ntree=501))

# Fit MSE = 0.9600181 
# Fit MSE = 0.9932792 
# Fit percent variance explained = 90.58 
# Median permuted MSE = 1.056655 
# Median permuted percent variance explained = 89.9 
# Median cross-validation RMSE = 1.016558 
# Median cross-validation MBE = -0.003313227 
# Median cross-validation MAE = 0.7115378 
# Range of ks p-values = 3.455415e-09 0.0004783882 
# Range of ks D statistic = 0.06617647 0.1029412 
# RMSE cross-validation error variance = 0.001933153 
# MBE cross-validation error variance = 0.0009840741 
# MAE cross-validation error variance = 0.0006018589 

# Mixed effects model -----------------------------------------------------
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc)
library(sjlabelled)
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library(effects)
library(sjstats) #use for r2 functions
library("scales")
library(performance)
library(lme4)

# install.packages("glmmTMB")
# library(glmmTMB)

## join swp and plantids back, drop nas in SWP
dataall <- cbind(plantids, rf.data.imputedRel) %>%
  drop_na(SWP) %>%
  select(PlantID, SWP, names(rf.data.red))

names(dataall)
corrf <- dataall[,c(10:12, 15:22)]
cor(corrf)

## rescale variables
rf.data.rescale <- dataall %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         SJC002_sc = rescale(SJC_002, to = c(0,1)),
         SJC002_POM001Combined_sc = rescale(SJC002_POM001Combined, to = c(0,1)),
         POM001_sc = rescale(POM001, to = c(0,1)),
         # RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1)))


## mixed effects model

## all combined - relative
mod1 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + 
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Damage2+
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod1)
anova(mod1)
check_singularity(mod1) ## False
icc(mod1, by_group = TRUE)
icc(mod1)
r2_nakagawa(mod1) ## 0.124

## all combined - Raw
mod1a <- lmer(formula = SWP ~  Q_sc +  Replacement + Species + Group + 
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Damage2+ Year +
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               # (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod1a)
anova(mod1a)
check_singularity(mod1) ## False
icc(mod1a, by_group = TRUE)
icc(mod1a)
r2_nakagawa(mod1a) ## 0.596

## pom and sjc combined
mod2 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + 
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Damage2+
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               # (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod2)
anova(mod2)
check_singularity(mod2) ## False
icc(mod2, by_group = TRUE)
icc(mod2)
r2_nakagawa(mod2) ## 0.

## individual outfall discharge

mod3 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + 
               DistToSJC002_sc+ Substrate+ POS+ SJC002_sc + POM001_sc + Damage2+
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               # (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod3)
anova(mod3)
check_singularity(mod3) ## False
icc(mod3, by_group = TRUE)
icc(mod3)
r2_nakagawa(mod3) ## 0.495


# Final Model  ------------------------------------------------------------

## join swp and plantids back, drop nas in SWP
dataall <- cbind(plantids, rf.data.imputedRel) %>%
  drop_na(SWP) %>%
  select(PlantID, SWP, names(rf.data.red))

names(dataall)
corrf <- dataall[,c(10:12, 15:22)]
cor(corrf)

## rescale variables
rf.data.rescale <- dataall %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         SJC002_sc = rescale(SJC_002, to = c(0,1)),
         SJC002_POM001Combined_sc = rescale(SJC002_POM001Combined, to = c(0,1)),
         POM001_sc = rescale(POM001, to = c(0,1)),
         # RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1)))


## all combined - Raw
mod1 <- lmer(formula = SWP ~  Q_sc +  Replacement + Species + Group + 
                DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Damage2+ Year +
                # (1|Year) + ## remove due to low ICC
                # (1|Group) + ## remove due to low ICC
                (1|Species) +
                (1|Season),
              data    = rf.data.rescale) 

summary(mod1)
anova(mod1)
check_singularity(mod1) ## False
icc(mod1, by_group = TRUE)
icc(mod1)
r2_nakagawa(mod1) ## 0.720

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,    #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

ests <- sjPlot::plot_model(mod1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of SWP")

ests
file.name1 <- paste0(out.dir, "effect_sizes_drivers_of_SWP_final.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

### plot random effects
random <- plot_model(mod1, type="re",
                     vline.color="#A9A9A9", dot.size=1.5,
                     show.values=T, value.offset=.2, show.p=TRUE)
random

out.filename <- paste0(out.dir,"00_mixed_mod_SWP_random_effects_inc_damage.jpg")
ggsave(random, file = out.filename, dpi=300, height=10, width=12)

## plot Q per species

## first predict on model to use in ggplot

rf.data.rescalex <- rf.data.rescale %>% 
  mutate(fit.m = predict(mod1, re.form = NA),
         fit.c = predict(mod1, re.form = NULL))

head(rf.data.rescalex)
str(rf.data.rescalex)

m1 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = SJC002_POM001Combined, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Combined Q (MGD)") +
  scale_y_continuous(name = "SWP: Mixed Model Fit")
  # geom_line(aes(col = Species), linewidth = 1)  
 

m1

out.filename <- paste0(out.dir,"04_mixed_mod_combQ_per_species.jpg")
ggsave(m1, file = out.filename, dpi=300, height=10, width=12)

## plot distance per species
m2 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = DistToSJC002, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Distance to SJC 002") +
  scale_y_continuous(name = "SWP: Mixed Model Fit")
# geom_line(aes(col = Species), linewidth = 1)  


m2

out.filename <- paste0(out.dir,"04_mixed_mod_distance_to_SJC_per_species.jpg")
ggsave(m2, file = out.filename, dpi=300, height=10, width=12)

## fixed effects
## gage discharge
Q <- plot_model(mod1, type="pred", terms = c("Q_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_Q.jpg")
ggsave(Q, file = out.filename, dpi=300, height=4, width=6)

## combined outfall discharge
of
of <- plot_model(mod1, type="pred", terms = c("SJC002_POM001Combined_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_combined_outfall.jpg")
ggsave(of, file = out.filename, dpi=300, height=4, width=6)

## replacement
rep
rep <- plot_model(mod1, type="pred", terms = c("Replacement"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_replacement.jpg")
ggsave(rep, file = out.filename, dpi=300, height=4, width=6)

## Year
year
year <- plot_model(mod1, type="pred", terms = c("Year"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_year.jpg")
ggsave(year, file = out.filename, dpi=300, height=4, width=6)

## group
group <- plot_model(mod1, type="pred", terms = c("Group"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_group.jpg")
ggsave(group, file = out.filename, dpi=300, height=4, width=6)

## distance
dist
dist <- plot_model(mod1, type="pred", terms = c("DistToSJC002_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_distance.jpg")
ggsave(dist, file = out.filename, dpi=300, height=4, width=6)

## substrate
sub <- plot_model(mod1, type="pred", terms = c("Substrate"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_substrate.jpg")
ggsave(sub, file = out.filename, dpi=300, height=4, width=10)

# position
pos
pos <- plot_model(mod1, type="pred", terms = c("POS"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_POS.jpg")
ggsave(pos, file = out.filename, dpi=300, height=4, width=6)

## damage
dam <- plot_model(mod1, type="pred", terms = c("Damage2"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_damage.jpg")
ggsave(dam, file = out.filename, dpi=300, height=4, width=6)

# Probability of increase in SWP ------------------------------------------

## join swp and plantids back, drop nas in SWP
dataall <- cbind(plantids, rf.data.imputedRel) %>%
  drop_na(SWP)

## histograms for distribution
hist(dataall$SJC002_POM001Combined)
hist(dataall$SWP)

## get baseline values
names(dataall)

## baseline data years
basedata <- dataall %>%
  filter(Year %in% c(2018, 2019, 2020))

sum(is.na(basedata$SWP))

median(basedata$SWP) ## 9
median(basedata$SJC002_POM001Combined) ## mean - 11.63, median 6.5

basedata %>% group_by(Year, Season) %>% summarise(meds = median(SWP))

# 1 2018  Fall    9.95
# 2 2019  Fall   10   
# 3 2019  Spring  6.5 
# 4 2020  Fall   10   
# 5 2020  Spring  6 

basedata %>% group_by(Year, Season) %>% summarise(meds = median(SJC002_POM001Combined))

# 1 2018  Fall    5.40
# 2 2019  Fall   13.3 
# 3 2019  Spring 13.3 
# 4 2020  Fall    6.50
# 5 2020  Spring  6.50

## means per year

currentdata <- dataall %>%
  filter(!Year %in% c(2018, 2019, 2020))

median(currentdata$SWP) ## 8.5
median(currentdata$SJC002_POM001Combined) ## median - 9.12

currentdata %>% group_by(Year, Season) %>% summarise(meds = median(SWP))

# 1 2021  Fall    9.88
# 2 2021  Spring  8   
# 3 2022  Fall    8.5 
# 4 2022  Spring  7.5 

currentdata %>% group_by(Year, Season) %>% summarise(meds = median(SJC002_POM001Combined))

# 1 2021  Fall   10.0 
# 2 2021  Spring 10.0 
# 3 2022  Fall    8.46
# 4 2022  Spring  8.46

## delta Q

## median monthly Q as baseline

basedatax <- basedata %>%
  ungroup() %>%
  group_by(Month) %>%
  summarise(Baseline = median(SJC002_POM001Combined))

basedatax

## median monthly Q current

currentdatax <- currentdata %>%
  ungroup() %>%
  # select(Species, PlantID, Month, Year, SJC002_POM001Combined) %>% 
  rename(Current = SJC002_POM001Combined) %>%
  distinct() 

dim(currentdatax)

## join together

deltaDat <- full_join(basedatax, currentdatax, by = "Month") %>%
  mutate(Delta = Current-Baseline) %>% select(-SWP) %>% ## calculate delta, remove raw SWP
  pivot_longer(c(Baseline, Current, Delta), names_to = "Type", values_to = "CombQ") %>%
  mutate(Season = case_when(Month %in% c("05","06","07","08","09", "10") ~ "Fall",
                            Month %in% c("11","12","01","02","03","04") ~ "Spring"))

names(deltaDat)
## delta SWP

## median monthly Q as baseline

basedataSWP <- basedata %>%
  ungroup() %>%
  group_by(PlantID, Month, Season) %>%
  ## assign season to months
  summarise(Baseline = median(SWP)) %>%
  mutate(keep = case_when(Month %in% c("05","06","07","08","09", "10") & Season == "Fall" ~ "Y",
                          Month %in% c("11","12","01","02","03","04") & Season == "Spring" ~ "Y")) %>%
  drop_na(keep) %>% select(-keep)

basedataSWP

## median monthly Q current

currentdataSWP <- currentdata %>%
  ungroup() %>%
  # select(Species, PlantID, Month, Season, Year, SWP) %>%
  rename(Current = SWP) %>%
  ## assign season to months
  mutate(keep = case_when(Month %in% c("05","06","07","08","09", "10") & Season == "Fall" ~ "Y",
                          Month %in% c("11","12","01","02","03","04") & Season == "Spring" ~ "Y")) %>%
  drop_na(keep) %>% select(-keep) %>%
  distinct()

currentdataSWP

## join together
## fall = 5,6,7,8,9, 10
## spring = 11,12,1,2,3,4
str(deltaDatSWP)

deltaDatSWP <- full_join(basedataSWP, currentdataSWP, by = c("Month", "Season", "PlantID")) %>%
  mutate(Delta = Current-Baseline) %>% select(-SJC002_POM001Combined) %>% ## calculate delta, remove raw Q
  pivot_longer(c(Baseline, Current, Delta), names_to = "Type", values_to = "SWP")

names(deltaDatSWP)
names(deltaDat)
## join together 

nams <- names(deltaDatSWP)[-c(24)]
nams

# AllDeltaDat <- full_join(deltaDatSWP, deltaDat, by = c("Month", "Type", "Year", "Season", "PlantID", "Species", "SWP.relative")) %>%
#   distinct()

AllDeltaDat <- full_join(deltaDatSWP, deltaDat, by = c(nams)) %>%
  distinct()

names(AllDeltaDat)

## plot

ggplot(data = AllDeltaDat, aes(x = CombQ, y = SWP, group = Type, colour = Type)) +
  geom_point() +
  facet_wrap(~Season)

## plot only delta, by species

p1 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP, group = Species, colour = Species)) +
  geom_smooth()  +
  facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p1

file.name1 <- "Figures/04_delta_SWP_combq_by_species.jpg"
ggsave(p1, filename=file.name1, dpi=300, height=5, width=8)

## plot delta by year

p2 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP, group = Year, colour = Year)) +
  geom_smooth()  +
  facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p2

file.name1 <- "Figures/04_delta_SWP_combq_by_year.jpg"
ggsave(p2, filename=file.name1, dpi=300, height=5, width=8)

## plot delta all data

p3 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP)) +
  geom_smooth()  +
  facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p3

file.name1 <- "Figures/04_delta_SWP_combq_all.jpg"
ggsave(p3, filename=file.name1, dpi=300, height=5, width=8)


# Delta models ------------------------------------------------------------

## we want probability of any increase in SWP 

## transform to 0 = decrease (less stress), 1 = increase (more stress)
AllDeltaDat
binData <- AllDeltaDat %>%
  filter(Type == "Delta", !Group == "G5") %>%
  mutate(Stress = ifelse(SWP > 0,1,0)) %>%
  select(PlantID, Species, Group, Month, Year, SWP, CombQ, Stress, Season) %>%
  # mutate(Species = as.character(Species)) %>%
  drop_na(SWP)

## find NAs
# ind <- which(is.na(binData$SWP))
# binData[ind,] ## where tree has been replaced between baseline and current, can remove

## empty DFs

df <- data.frame(matrix(ncol = 4)) 
names(df) <- c("Species", "AIC", "Pval", "McFadsR2")
spDataPx <- NULL
## loop over species
s=1
## define species
species <- unique(binData$Species)

species
## 95% confidence intervals
critval <- 1.96 ## approx 95% CI

## start loop
s = 2
for(s in 1:length(species)) {

  spData <- binData %>%
    filter(Species == species[s])

  mod <- glm(Stress~CombQ, family=binomial(link="logit"), data = spData)
  ## summary
  modsum <- summary(mod)
  ## get vals into df
  df[s,1] <- paste(species[s])
  df[s,2] <- mod$aic
  df[s,3] <- modsum$coefficients[8]
  df[s,4] <- 1-modsum$deviance/modsum$null.deviance ## 0.06
  
  ## predict glm
  preds <- as.data.frame(predict.glm(mod, type = "response", se.fit = T)) %>%
    rename(ProbabilityOfStress = 1) %>%
    mutate(upr = ProbabilityOfStress + (critval * se.fit),
           lwr = ProbabilityOfStress - (critval * se.fit)) 

  ## join with data
  
  spDataP <- cbind(spData, preds)
  
  ## accumulate DF
  spDataPx <- rbind(spDataPx, spDataP)
  
}

head(df)
head(spDataPx)

## plot

s1 <- ggplot(spDataPx, aes(x = CombQ, y = ProbabilityOfStress, group = Species, colour = Species)) +
  geom_line()  +
  # facet_wrap(~Season, scale = "free_x") +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress")

s1

file.name1 <- "Figures/04_probability_of_stress_per_species.jpg"
ggsave(s1, filename=file.name1, dpi=300, height=5, width=8)


## model altogether
mod <- glm(Stress~CombQ, family=binomial(link="logit"), data = binData)
## summary
modsum <- summary(mod)
## get vals into df
mod$aic
modsum$coefficients[8]
1-modsum$deviance/modsum$null.deviance ## 0.09

## predict glm
preds <- as.data.frame(predict(mod, type = "response",  se.fit = TRUE)) %>%
  rename(ProbabilityOfStress = 1) %>% 
  ## confidence intervals
  mutate(upr = ProbabilityOfStress + (critval * se.fit),
         lwr = ProbabilityOfStress - (critval * se.fit)) 

head(preds)

## join with data

binDataP <- cbind(binData, preds)

## plot

s2 <- ggplot(binDataP, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  # facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress")

s2

file.name1 <- "Figures/04_probability_of_stress_overall.jpg"
ggsave(s2, filename=file.name1, dpi=300, height=5, width=8)

## plot with secondary axis to show observed Q values

names(binDataP)
names(binDataC)
## get current values of combined Q
binDataC <- AllDeltaDat %>%
  filter(!Group == "G5") %>%
  # mutate(TypeQ = Type) %>%
  select(PlantID, Species, Group, Month, Year, Season, Type, CombQ) %>%
  distinct() %>%
  pivot_wider(names_from = "Type", values_from = "CombQ") %>%
  rename(BaselineQ = Baseline, CurrentQ = Current, DeltaQ = Delta) #%>%
  # mutate(Species = as.character(Species)) %>%
  # drop_na(SWPCurrent)

binDataC
## join with predictions

binDataA <- full_join(binDataC, binDataP, by = c("PlantID", "Species", "Group", "Month", "Year",   "Season", "DeltaQ" = "CombQ")) %>%
  drop_na(SWP)

## get multiplier for secondary axis

min(binDataA$CurrentQ) -min(binDataA$DeltaQ) ## 22.9

## plot with secondary axis

s3 <- ggplot(binDataA, aes(x = CurrentQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  # geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  # facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name="Combined Q (MGD)") +
  # scale_x_continuous(sec.axis = ~. + 22.9, name="Combined Q (MGD)" )+
  scale_y_continuous(name = "Probability of Stress")

s3

file.name1 <- "Figures/04_probability_of_stress_overall_currentQ.jpg"
ggsave(s2, filename=file.name1, dpi=300, height=5, width=8)

## delta doesn't go hand in hand with current Q
binDataA %>% group_by(Group) %>% summarise(MeanBLQ = mean(BaselineQ))
##12-13

binDataA %>% group_by(Group) %>% summarise(MeanBLQ = mean(CurrentQ))
## 10-11

# mean delta of ~ 2-3

newd <- data.frame(CombQ=seq(-5,-2, 0.5))

pnewd <- as.data.frame(predict(mod, newdata = newd, type = "response"))

dnew <- cbind(newd, pnewd)

dnew

write.csv(dnew, "output_data/04_predicted_probs_current_delta.csv")

## plot with current delta ranges -2:-3

s3 <- ggplot(binDataP, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  annotate("rect", xmin = -3, xmax = -2, ymin = -Inf, ymax = Inf,
           alpha = .2, fill='red') +
  geom_text(aes(x=0, label="Baseline Discharge", y=0.70), colour="forestgreen", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=-2, label="Current Discharge Range", y=0.73), colour="red", angle=90, vjust = 1.3, size=4)+
  geom_vline(xintercept = 0, col = "forestgreen", lty = "dashed") +
  # facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s3

file.name1 <- "Figures/04_probability_of_stress_overall_current_delta_range.jpg"
ggsave(s3, filename=file.name1, dpi=300, height=5, width=8)

# Mixed model on prob of stress -------------------------------------------



## mixed model
## rescale variables
rf.data.rescale <- binData %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         SJC002_sc = rescale(SJC_002, to = c(0,1)),
         # SJC002_POM001Combined_sc = rescale(SJC002_POM001Combined, to = c(0,1)),
         POM001_sc = rescale(POM001, to = c(0,1)),
         # RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1)))

## all combined - Raw
mod1 <- glmer(formula = Stress ~  Q_sc + Species + Group + DistToSJC002_sc+ Substrate+
                Substrate+ POS+ CombQ + Damage2+ Year + 
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               # (1|Year) +
               (1|Season),
             data = rf.data.rescale,
             family = binomial)
summary(mod1)
# +  Replacement + Species + Group + 
#   DistToSJC002+ Substrate+ POS+ CombQ + Damage2+ Year +
summary(mod1)
anova(mod1)
check_singularity(mod1) ## False
icc(mod1, by_group = TRUE)
icc(mod1)
r2_nakagawa(mod1) ## 0.323

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,    #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size

ests <- sjPlot::plot_model(mod1, 
                           show.values=TRUE, show.p=TRUE,
                           title="Drivers of SWP")

ests
file.name1 <- paste0(out.dir, "effect_sizes_drivers_of_SWP_relative_inc_damage.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

