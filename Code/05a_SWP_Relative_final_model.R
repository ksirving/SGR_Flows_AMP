## SWP relative final models

library(tidyverse)
library(tidylog)
library(lubridate)

## random forest
library(rfUtilities)
library(randomForest)
library(ks)
library(sm)
library(caret)
library(rpart)
library(scales)

## mixed model
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc)
library(sjlabelled)
library(lmerTest) ## from Luke et al 2016 - evaluating significance in linear mixed-effects models in R
library(effects)
library(sjstats) #use for r2 functions
library(performance)
library(lme4)

## output file for figures
out.dir <- "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SGR_Flows_AMP/Figures/"

# upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups_distance.RData")

## check damaged trees

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
sum(bioMeanQ_long_dist$Damage == "Y") ## 505a_relative

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

# Relative  -----------------------------------------------

## remove columns and make distinct
bioMeanQ_long_distx <- bioMeanQ_long_dist %>%
  group_by(PlantID) %>%
  select(PlantID:Year, SWP , Replacement, Damage) %>%
  distinct()
bioMeanQ_long_distx

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
  ungroup() 

## join back to main df

bioMeanQ_long_distJoin <- full_join(bioMeanQ_long_distRel, bioMeanQ_long_dist, by = c("PlantID","Species", "Season", "Group", "LifeForm", "Year", "SWP", "Replacement", "Damage"))

# SWP model on relative change---------------------------------------------------------------

rf.data <- bioMeanQ_long_distJoin %>%
  dplyr::select(SWP.relative, SWP, PlantID, Species:LifeForm, Month, Year, Season, POS, Substrate,
                RainFallMod:TempMeanModF, Replacement, Damage2, SJC002_POM001Combined:WN002, ## include damaged trees in all years after damage occurred
                Q, DistToSJC002, Source, lag.SWP) %>%
  mutate(POS = as.factor(POS),
         Species = as.factor(Species),
         Group = as.factor(Group),
         LifeForm = as.factor(LifeForm),
         Month = as.factor(Month),
         Substrate = as.factor(Substrate),
         Damage2 = as.factor(Damage2)) %>%
  filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  mutate(SWP.relative = ifelse(is.na(lag.SWP), 1, SWP.relative)) %>%
  select(-Source, -lag.SWP)

names(rf.data)
dim(rf.data)

## define plantids 
plantids <- as.data.frame(rf.data$PlantID) %>%
  rename(PlantID = 1)

plantids$SWP <- rf.data$SWP
plantids

## remove plantid from df

rf.data <- rf.data %>% select(-PlantID, -SWP)


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
rf.data.imputedRel <- rfImpute(SWP.relative ~ ., rf.data)
head(rf.data.imputedRel)

## make y data compatible
rf.data.val <- rf.data.imputedRel %>%
  rename(y = SWP.relative) %>% as.data.frame() %>%
  droplevels()

## run the model
rf1 <- randomForest(y~., data=rf.data.val, importance = T)
rf1
mean(rf1$rsq) ##  - 76%

## importance quick look
varImpPlot(rf1) ## season
importance(rf1, type = 1) 


# Remove negative importance variables ------------------------------------

## extract importance
imp <- as.data.frame(importance(rf1, type = 1))

## get all negative importance vars
vars <- ifelse(imp$`%IncMSE` > 0, rownames(imp), NA)
## remove neg importance vars
vars <- vars[!is.na(vars)]

## remove from dataset
rf.data.red <- rf.data.imputedRel %>%
  select(SWP.relative, all_of(vars), -LifeForm) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

## format data
rf.data.val <- rf.data.red %>%
  rename(y = SWP.relative) %>% as.data.frame()

## run reduced model
rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ##  77%
rf
save(rf, file = "models/05a_SWP_relative_final_rf.RData")

load(file = "models/05a_SWP_relative_final_rf.RData")
rf
# Relative importance -----------------------------------------------------

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

file.name1 <- "Figures/05a_relative_importance_SWP_relative_final_model.jpg"
ggsave(i1, filename=file.name1, dpi=300, height=5, width=8)


# Partial Plots -----------------------------------------------------------

jpeg(paste0( "Figures/05a_relative_SWP_relative_Season.jpg"))
partialPlot(rf, rf.data.val, Season)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Distance.jpg"))
partialPlot(rf, rf.data.val, DistToSJC002)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Q.jpg"))
partialPlot(rf, rf.data.val, Q)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Species.jpg"))
partialPlot(rf, rf.data.val, Species)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_SJC.jpg"))
partialPlot(rf, rf.data.val, SJC_002)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Replacement.jpg"))
partialPlot(rf, rf.data.val, Replacement)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Group.jpg"))
partialPlot(rf, rf.data.val, Group)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_POS.jpg"))
partialPlot(rf, rf.data.val, POS)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Substrate.jpg"))
partialPlot(rf, rf.data.val, Substrate)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Year.jpg"))
partialPlot(rf, rf.data.val, Year)
dev.off()

jpeg(paste0( "Figures/05a_relative_SWP_relative_Damage.jpg"))
partialPlot(rf, rf.data.val, Damage2)
dev.off()


# Mixed models ------------------------------------------------------------

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


## pom and sjc combined
mod2 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + 
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Damage2+
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod2)
anova(mod2)
check_singularity(mod2) ## False
icc(mod2, by_group = TRUE)
icc(mod2)
r2_nakagawa(mod2) ## 0.12

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
r2_nakagawa(mod3) ## 0.061


# Cross Validation --------------------------------------------------------

## split into training and testing

setsize <- floor(nrow(rf.data.red)*0.8)
index <- sample(1:nrow(rf.data.red), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(rf.data.val)

(rf.SWP <- rf.crossValidation(rf, rf.data.val[, 2:13], ## check indexing!!!!
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




