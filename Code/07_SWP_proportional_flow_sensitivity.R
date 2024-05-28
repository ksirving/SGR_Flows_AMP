## final models for SWP

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

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## range function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# Upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
head(bioMeanQ_long_dist)
names(bioMeanQ_long_dist)

sum(is.na(bioMeanQ_long_dist$Proportion)) ## 6

## remove columns and make distinct
bioMeanQ_long_distx <- bioMeanQ_long_dist %>%
  group_by(PlantID) %>%
  select(PlantID:Year, SWP , Replacement, Damage) %>%
  distinct()
bioMeanQ_long_distx

## if damage then plantid for remaining surveys should be also damaged

## add the column with repeated damage
damSeq <- bioMeanQ_long_distx %>%
  drop_na(SWP) %>%
  group_by(PlantID) %>%
  mutate(firstY = which.max(Damage == "Y")) %>% ## index of first damage
  mutate(firstY = ifelse(firstY == 1, 0, firstY)) %>% ## change 1 to 0 
  mutate(Damage2 = ifelse(seq_along(Damage) >= firstY | Damage == "Y", "Y", "N")) %>% ## sequence from first Y to repeat until next group
  mutate(Damage2 = ifelse(firstY == 0, "N", Damage2)) ## change all values with first y as 0 from N to Y (no damage in these groups)
damSeq

## join back to main df

bioMeanQ_long_distJoin <- full_join(damSeq, bioMeanQ_long_dist, by = c("PlantID","Species", "Season", "Group", "LifeForm", "Year", "SWP", "Replacement", "Damage"))


# Random Forest model on Raw values--------------------------------------------------------------

## take all columns needed, change categroical data to factors

## need to run through proportions, filter here

pros <- unique(bioMeanQ_long_distJoin$Proportion)
pros
p=1

p

## define empty df
dfx <- NULL

for(p in 1:length(pros)) {
  
  rf.data <- bioMeanQ_long_distJoin %>%
    filter(Proportion == pros[p]) %>%
    dplyr::select(SWP, PlantID, Species:LifeForm, Month, Year, Season, POS, Substrate,
                  RainFallMod:TempMeanModF, Replacement, Damage2, SJC002_POM001Combined:WN002, ## include damaged trees in all years after damage occurred
                  Q, DistToSJC002, Source) %>%
    mutate(POS = as.factor(POS),
           Species = as.factor(Species),
           Group = as.factor(Group),
           LifeForm = as.factor(LifeForm),
           Month = as.factor(Month),
           Substrate = as.factor(Substrate),
           Damage2 = as.factor(Damage2)) %>%
    filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
    select(-Source)
  
  
  ## remove plantid from df
  
  rf.data <- rf.data %>% ungroup() %>% select(-PlantID) %>% drop_na(SWP)
  dim(rf.data)
  
  # Impute Values -----------------------------------------------------------
  
  ## impute missing values
  rf.data.imputed <- rfImpute(SWP ~., rf.data)
  head(rf.data.imputed)
  dim(rf.data.imputed)
  
  save(rf.data.imputed, file = paste0("ignore/07_rf_data_imputed_SWP_raw_", pros[p], ".RData"))
  

# Random Forest Model -----------------------------------------------------

## make y data compatible
rf.data.val <- rf.data.imputed %>%
  # select(- PlantIDNum) %>%
  rename(y = SWP) %>% as.data.frame() %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

## run model
rf <- randomForest(y~., data=rf.data.val, importance = T)

# Remove negative importance variables ------------------------------------

## get importance from model
imp <- as.data.frame(importance(rf, type = 1))

## get all variables over 0 (NA negative variables)
vars <- ifelse(imp$`%IncMSE` > 0, rownames(imp), NA) 
vars <- vars[!is.na(vars)] ## remove NAs - neg vars

## only keep positive importance vars from main df
rf.data.red <- rf.data.imputed %>%
  select(SWP, all_of(vars), -LifeForm) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

# save(rf.data.red, file = "ignore/07_SWP_red_imputed_rf_data.RData")

## format data for reduced model
rf.data.val <- rf.data.red %>%
  rename(y = SWP) %>% as.data.frame()

## run model
rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)

# Variable importance -----------------------------------------------------
## get mean and scale and % to make relative

VarImp <- as.data.frame(rf$importance) 
colnames(VarImp)[1] <- "DecreaseAcc"
VarImp$Variable <- rownames(VarImp)

## get perc means
ImpMean <- VarImp %>% 
  # group_by(Variable) %>%
  # summarise(MeanImp = mean(MeanDecreaseAccuracy)) %>%
  mutate(MeanImpScaled = rescale(DecreaseAcc)) %>%
  mutate(MeanImpPerc = (MeanImpScaled/sum(MeanImpScaled))*100)

## humanise variables
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

df <- ImpMean %>%
  mutate(RFVarExpl = mean(rf$rsq),
         Proportion = pros[p])

dfx <- bind_rows(dfx, df)

}

head(dfx)

write.csv(dfx, "ignore/07_proportional_varImp_varexpl.csv")

dfx <- read.csv("ignore/07_proportional_varImp_varexpl.csv")

## upload original variable importance
origVarImp <- read.csv("output_data/05_SWP_var_imp.csv") %>%
  mutate(Proportion = "Constant (0.001)")

## add to sensitivity analysis

head(dfx)
head(origVarImp)

dfx <- bind_rows(dfx, origVarImp)

## plot proportions Vs variable importance
unique(dfx$Proportion)

# factorize proportiona flow

dfx <- dfx %>%
  mutate(ProportionF = factor(Proportion, levels = c("Constant (0.001)", "QPro1",  "QPro2",  "QPro5",  "QPro10", "QPro20"),
                              labels = c("0.01", "1%", "2%", "5%", "10%", "20%")))

dfx$VariableHuman <- factor(dfx$VariableHuman, levels=ImpMean[order(dfx$MeanImpPerc,decreasing=F),]$VariableHuman)

dfx
## plot all variables according to proportion of flow

i1 <- ggplot(dfx, aes(y=MeanImpPerc, x=ProportionF)) +
  geom_point() +
  scale_x_discrete("Proportional Flow to Group 5") +
  scale_y_continuous("Relative Importance (%)") +
  theme_bw() +
  facet_wrap(~VariableHuman, scales = "free")
i1

file.name1 <- "Figures/07_relative_importance_SWP_sens_analysis.jpg"
ggsave(i1, filename=file.name1, dpi=300, height=5, width=8)

## plot var figure - colour points as proprtion 
i2 <- ggplot(dfx, aes(x=MeanImpPerc, y=VariableHuman)) +
  geom_point(aes(col = ProportionF, size = RFVarExpl)) +
  scale_x_continuous("Relative Importance (%)") +
  scale_y_discrete("") +
  theme_bw() ## add legend title
i2

file.name1 <- "Figures/07_relative_importance_SWP_comparison_plot.jpg"
ggsave(i2, filename=file.name1, dpi=300, height=5, width=8)


