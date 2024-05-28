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

# Upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
head(bioMeanQ_long_dist)
names(bioMeanQ_long_dist)

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

rf.data <- bioMeanQ_long_distJoin %>%
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

## make plantID a numeric

rf.data$PlantIDNum <- as.numeric(factor(rf.data$PlantID))

## define plantids 
plantids <- rf.data %>%
  select(PlantID, Species, Group, PlantIDNum) 
  

plantidsD <- rf.data %>%
  select(PlantID, Species, Group, PlantIDNum) %>%
  distinct()

save(plantidsD, file = "ignore/05_plantids.RData")

load(file = "ignore/05_plantids.RData")

str(rf.data)
## remove plantid from df

rf.data <- rf.data %>% ungroup() %>% select(-PlantID) %>% drop_na(SWP)
dim(rf.data)

# Impute Values -----------------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## impute missing values
rf.data.imputed <- rfImpute(SWP ~., rf.data)
head(rf.data.imputed)
dim(rf.data.imputed)

save(rf.data.imputed, file = "ignore/05_rf_data_imputed_SWP_raw.RData")

load(file = "ignore/05_rf_data_imputed_SWP_raw.RData") #### change WN to 0 for grps 1,2,3

head(rf.data.imputed)

rf.data.imputed <- rf.data.imputed %>%
  mutate(WN001 = case_when(Group %in% c("G1", "G2", "G3") ~ 0, .default = WN001),
         WN002 = case_when(Group %in% c("G1", "G2", "G3") ~ 0, .default = WN002))

# Random Forest Model -----------------------------------------------------

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

## make y data compatible
rf.data.val <- rf.data.imputed %>%
  select(- PlantIDNum) %>%
  rename(y = SWP) %>% as.data.frame() %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

## run model
rf <- randomForest(y~., data=rf.data.val, importance = T)
rf ## 90
mean(rf$rsq) ## 0.90

## importance
varImpPlot(rf)
importance(rf, type = 1) ## season



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

save(rf.data.red, file = "ignore/05_SWP_red_imputed_rf_data.RData")

## format data for reduced model
rf.data.val <- rf.data.red %>%
  rename(y = SWP) %>% as.data.frame()

## run model
rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ## 0.92
rf

## save final reduced model

# save(rf, file = "models/05_SWP_final_rf_Raw_updated.RData")

save(rf, file = "models/05_SWP_final_rf_Raw.RData")

# Random Forest without damaged trees -------------------------------------

head(rf.data.imputed)

rf.data.valD <- rf.data.red %>%
  # select(- PlantIDNum) %>%
  rename(y = SWP) %>% as.data.frame() %>%
  filter(Damage2 == "N") %>% ## 7% of trees removed, 18 altogether 
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

## run model
rfD <- randomForest(y~., data=rf.data.valD, importance = T)
rfD ## 92.05
mean(rfD$rsq) ## 0.92

## importance
varImpPlot(rfD)
importance(rfD, type = 1) ## season

# Variable importance -----------------------------------------------------

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
## add mod peformance
ImpMean <- ImpMean %>%
  mutate(RFVarExpl = mean(rf$rsq))

write.csv(ImpMean, "output_data/05_SWP_var_imp.csv")

## plot
i1 <- ggplot(ImpMean, aes(x=MeanImpPerc, y=VariableHuman)) +
  geom_point() +
  scale_x_continuous("Relative Importance (%)") +
  scale_y_discrete("") +
  theme_bw()
i1

file.name1 <- "Figures/05_relative_importance_SWP_raw_final_model.jpg"
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


# Probability of increase in SWP ------------------------------------------
load(file = "ignore/05_rf_data_imputed_SWP_raw.RData")
load(file = "ignore/05_plantids.RData")

dim(rf.data.imputed)
## join swp and plantids back, drop nas in SWP
dataall <- inner_join(plantidsD, rf.data.imputed, by = c("PlantIDNum", "Group", "Species")) %>%
  drop_na(SWP) %>%
  distinct()

## histograms for distribution
hist(dataall$SJC002_POM001Combined)
hist(dataall$SWP)

## get baseline values
names(dataall)

## baseline data years - dry years only
basedata <- dataall %>%
  filter(Year %in% c(2020))
basedata

median(basedata$SWP) ## 8.75
median(basedata$SJC002_POM001Combined) ## mean - 11.63, median 6.5

## means per year

currentdata <- dataall %>%
  filter(Year %in% c(2021, 2022))

median(currentdata$SWP) ## 8.5
median(currentdata$SJC002_POM001Combined) ## median - 9.12

## delta Q

## median monthly Q as baseline

basedatax <- basedata %>%
  ungroup() %>%
  group_by(Month) %>%
  summarise(Baseline = median(SJC002_POM001Combined))

basedatax

## median monthly Q current
currentdata
currentdatax <- currentdata %>%
  ungroup() %>%
  dplyr::select(Species, PlantID, Month, Year, SJC002_POM001Combined, Group) %>%
  rename(Current = SJC002_POM001Combined) %>%
  distinct() 

dim(currentdatax)
names(currentdatax)
## join together

deltaDat <- full_join(basedatax, currentdatax, by = "Month") %>%
  # select(-c(POS:Damage2, DischargeSJC002_POM001_WN001Combined:DistToSJC002)) %>%
  mutate(Delta = Current-Baseline) %>%  ## calculate delta, remove raw SWP
  pivot_longer(c(Baseline, Current, Delta), names_to = "Type", values_to = "CombQ") %>%
  ## assign season to months
  mutate(Season = case_when(Month %in% c("05","06","07","08","09", "10") ~ "Fall",
                            Month %in% c("11","12","01","02","03","04") ~ "Spring")) %>%
  distinct()


## median flows
## Baseline
deltaDat %>% group_by(Season) %>% filter(Type == "Baseline") %>% summarise(med = median(CombQ)) 

# 1 Fall    7.3
# 2 Spring 24.1

# Current
deltaDat %>% group_by(Season) %>% filter(Type == "Current") %>% summarise(med = median(CombQ)) 

# 1 Fall    8.46
# 2 Spring 10.2 

## delta SWP

## median monthly Q as baseline

basedataSWP <- basedata %>%
  ungroup() %>%
  dplyr::select(Species, PlantID, Month, Season, Year, SWP, Group) %>%
  group_by(PlantID, Month, Season) %>%
  ## assign season to months
  summarise(Baseline = median(SWP)) %>%
  mutate(keep = case_when(Month %in% c("05","06","07","08","09", "10") & Season == "Fall" ~ "Y",
                          Month %in% c("11","12","01","02","03","04") & Season == "Spring" ~ "Y")) %>%
  drop_na(keep) %>% dplyr::select(-keep) %>% ungroup() %>%  distinct()


basedataSWP

## median monthly Q current

currentdataSWP <- currentdata %>%
  ungroup() %>%
  dplyr::select(Species, PlantID, Month, Season, Year, SWP, Group) %>%
  rename(Current = SWP) %>%
  # ## assign season to months
  mutate(keep = case_when(Month %in% c("05","06","07","08","09", "10") & Season == "Fall" ~ "Y",
                          Month %in% c("11","12","01","02","03","04") & Season == "Spring" ~ "Y")) %>%
  drop_na(keep) %>% dplyr::select(-keep)  %>% distinct()


currentdataSWP
unique(currentdataSWP$Year)

## join together
## fall = 5,6,7,8,9, 10
## spring = 11,12,1,2,3,4
str(deltaDatSWP)

deltaDatSWP <- full_join(basedataSWP, currentdataSWP, by = c("Season", "PlantID", "Month")) %>%
  mutate(Delta = Current-Baseline) %>%  ## calculate delta, remove raw Q
  pivot_longer(c(Baseline, Current, Delta), names_to = "Type", values_to = "SWP")

names(deltaDatSWP)
names(deltaDat)

deltaDatSWP
deltaDat
## join together 
# AllDeltaDat <- full_join(deltaDatSWP, deltaDat, by = c("Month", "Type", "Year", "Season", "PlantID", "Species", "SWP.relative")) %>%
#   distinct()

AllDeltaDat <- inner_join(deltaDatSWP, deltaDat, by = c("PlantID", "Species", "Season","Year","Month", "Type", "Group")) %>%
  distinct() ## 6930

names(AllDeltaDat)



## plot only delta, by species

p1 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP, group = Species, colour = Species)) +
  geom_smooth()  +
  facet_grid(cols = vars(Season), rows = vars(Year), scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p1

file.name1 <- "Figures/05_delta_SWP_combq_by_species_dry_years.jpg"
ggsave(p1, filename=file.name1, dpi=300, height=5, width=8)

## plot delta by year

p2 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP, group = Year, colour = Year)) +
  geom_smooth()  +
  facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p2

file.name1 <- "Figures/05_delta_SWP_combq_by_year_dry_years.jpg"
ggsave(p2, filename=file.name1, dpi=300, height=5, width=8)

## plot delta all data

p3 <- ggplot(data = subset(AllDeltaDat, Type == "Delta"), aes(x = CombQ, y = SWP)) +
  geom_smooth()  +
  facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Delta SWP")

p3

file.name1 <- "Figures/05_delta_SWP_combq_all.jpg"
ggsave(p3, filename=file.name1, dpi=300, height=5, width=8)

## save delta data

save(AllDeltaDat, file = "ignore/06_delta_SWP.Rdata")

# Delta models ------------------------------------------------------------

## we want probability of any increase in SWP 

load(file = "ignore/06_delta_SWP.Rdata")

length(unique(binData$PlantID))

## transform to 0 = decrease (less stress), 1 = increase (more stress)
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

df <- data.frame(matrix(ncol = 6)) 
names(df) <- c("Species","Year",  "AIC", "Pval", "McFadsR2", "n")
dfx <- NULL

spDataPx <- NULL
## loop over species
s=1
## define species
species <- unique(binData$Species)

species

## define years

years <- unique(binData$Year)

years
## 95% confidence intervals
critval <- 1.96 ## approx 95% CI

## start loop
s = 2
y=1

for(y in 1:length(years)) {
  
  spData <- binData %>%
    filter(Year == years[y])
  
for(s in 1:length(species)) {
  
  spDatax <- spData %>%
    filter(Species == species[s])
  
    mod <- glm(Stress~CombQ, family=binomial(link="logit"), data = spDatax)
    ## summary
    modsum <- summary(mod)
    ## get vals into df
    
    df[s,1] <- paste(species[s])
    df[s,2] <- paste(years[y])
    df[s,3] <- mod$aic
    df[s,4] <- modsum$coefficients[8]
    df[s,5] <- 1-modsum$deviance/modsum$null.deviance ## 
    df[s,6] <- length(spDatax$SWP)
    
    ## predict glm
    preds <- as.data.frame(predict.glm(mod, type = "response", se.fit = T)) %>%
      rename(ProbabilityOfStress = 1) %>%
      mutate(upr = ProbabilityOfStress + (critval * se.fit),
             lwr = ProbabilityOfStress - (critval * se.fit)) 
    
    ## join with data
    
    spDataP <- cbind(spDatax, preds)
    
    ## accumulate DF
    spDataPx <- rbind(spDataPx, spDataP)
    
      }
  
  dfx <- rbind(dfx, df)
  
  }
  


head(dfx)
dfx
head(spDataPx)

write.csv(dfx, "output_data/05_SWP_pred_stress_per_species_dry_years.csv")

## plot

s1 <- ggplot(spDataPx, aes(x = CombQ, y = ProbabilityOfStress, group = Species, colour = Species)) +
  geom_line()  +
  facet_wrap(~Year, scale = "free_x") +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress")

s1

file.name1 <- "Figures/05_probability_of_stress_per_species_dry_years.jpg"
ggsave(s1, filename=file.name1, dpi=300, height=5, width=8)


## model altogether - dry years
## 2021
binData2021 <- binData %>%
  filter(Year == 2021)

mod2021 <- glm(Stress~CombQ, family=binomial(link="logit"), data = binData2021)
## summary
modsum <- summary(mod2021)
## get vals into df
mod2021$aic
modsum$coefficients[8]
1-modsum$deviance/modsum$null.deviance ## 0.07

## predict glm
preds2021 <- as.data.frame(predict(mod2021, type = "response",  se.fit = TRUE)) %>%
  rename(ProbabilityOfStress = 1) %>% 
  ## confidence intervals
  mutate(upr = ProbabilityOfStress + (critval * se.fit),
         lwr = ProbabilityOfStress - (critval * se.fit)) 

head(preds2021)

## join with data

binData2021 <- cbind(binData2021, preds2021)

## 2022
binData2022 <- binData %>%
  filter(Year == 2022)

mod2022 <- glm(Stress~CombQ, family=binomial(link="logit"), data = binData2022)
## summary
modsum <- summary(mod2022)
## get vals into df
mod2022$aic
modsum$coefficients[8]
1-modsum$deviance/modsum$null.deviance ## 0.05

## predict glm
preds2022 <- as.data.frame(predict(mod2022, type = "response",  se.fit = TRUE)) %>%
  rename(ProbabilityOfStress = 1) %>% 
  ## confidence intervals
  mutate(upr = ProbabilityOfStress + (critval * se.fit),
         lwr = ProbabilityOfStress - (critval * se.fit)) 

head(preds2022)

## join with data

binData2022 <- cbind(binData2022, preds2022)

## join 2021/22 together
binDataP <- rbind(binData2022, binData2021)

binDataP

## plot

s2 <- ggplot(binDataP, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  facet_wrap(~Year, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress")

s2

file.name1 <- "Figures/04_probability_of_stress_overall_dry_years.jpg"
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

min(binDataA$CurrentQ) -min(binDataA$DeltaQ) ## 26.37

## plot with secondary axis

s3 <- ggplot(binDataA, aes(x = CurrentQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  # geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  facet_wrap(~Year, scale = "free_x") +
  scale_x_continuous(name="Combined Q (MGD)") +
  # scale_x_continuous(sec.axis = ~. + 22.9, name="Combined Q (MGD)" )+
  scale_y_continuous(name = "Probability of Stress")

s3

file.name1 <- "Figures/05_probability_of_stress_overall_currentQ_dry_years.jpg"
ggsave(s3, filename=file.name1, dpi=300, height=5, width=8)

## delta doesn't go hand in hand with current Q
binDataA %>% group_by(Group) %>% summarise(MeanBLQ = mean(BaselineQ))
## 12-13
binDataA %>% group_by( Season) %>% summarise(MeanBLQ = mean(BaselineQ))

binDataA %>%  summarise(MeanBLQ = mean(BaselineQ))

binDataA %>% group_by(Group) %>% summarise(MeanBLQ = mean(CurrentQ))
## 10-11
binDataA %>% group_by(Season) %>% summarise(MeanBLQ = mean(CurrentQ))

binDataA %>%  summarise(MeanBLQ = mean(CurrentQ))

## median
binDataA %>% group_by(Group) %>% summarise(MedianBLQ = median(BaselineQ))
## 8-8.5
binDataA %>% group_by(Group,Season) %>% summarise(medianBLQ = median(BaselineQ))
## 4 - 24
binDataA %>%  summarise(medianBLQ = median(BaselineQ))
## 8

binDataA %>% group_by(Group) %>% summarise(medianBLQ = median(CurrentQ))
## 10
binDataA %>% group_by(Group, Season) %>% summarise(medianBLQ = median(CurrentQ))
## 9-11.6
binDataA %>%  summarise(medianBLQ = median(CurrentQ))
## 10.2

# median delta of ~ fall = +5, spring = - 12

newd <- data.frame(CombQ=seq(-5,-2, 0.5))

pnewd <- as.data.frame(predict(mod, newdata = newd, type = "response"))

dnew <- cbind(newd, pnewd)

dnew

write.csv(dnew, "output_data/05_predicted_probs_current_delta.csv")

## plot with current delta ranges -2:-3

s3 <- ggplot(binDataP, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
           alpha = .2, fill='red') +
  geom_text(aes(x=0, label="Baseline Discharge", y=0.70), colour="gray30", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=2.5, label="Current Discharge Range", y=0.70), colour="red", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=5, label="Fall Current Discharge", y=0.70), colour="forestgreen", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=-12, label="Spring Current Discharge", y=0.70), colour="forestgreen", angle=90, vjust = 1.3, size=4)+
  geom_vline(xintercept = 0, col = "gray30", lty = "dashed") +
  geom_vline(xintercept = 5, col = "forestgreen", lty = "dashed") + ## fall line
  geom_vline(xintercept = -12, col = "forestgreen", lty = "dashed") + ## spring line
  facet_wrap(~Year, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s3

file.name1 <- "Figures/05_probability_of_stress_overall_current_delta_range_dry_years.jpg"
ggsave(s3, filename=file.name1, dpi=300, height=5, width=8)


# Mixed Model -------------------------------------------------------------

## join swp and plantids back, drop nas in SWP
dataall <- inner_join(plantidsD, rf.data.imputed, by = c("PlantIDNum", "Group", "Species")) %>%
  drop_na(SWP) %>%
  distinct() %>%
  select(PlantID, SWP, names(rf.data.red)) %>%
  distinct()

names(dataall)
corrf <- dataall[,c(9:10, 13:20)]
cor(corrf)

dim(dataall)

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
mod1 <- lmer(formula = SWP ~  Q_sc +  Replacement +  Group + RainfallIntensity +
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
r2_nakagawa(mod1) ## 0.653

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
file.name1 <- paste0(out.dir, "05_effect_sizes_drivers_of_SWP_final.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

### plot random effects
random <- plot_model(mod1, type="re",
                     vline.color="#A9A9A9", dot.size=1.5,
                     show.values=T, value.offset=.2, show.p=TRUE)
random

## species plot
speciesR <- random[[1]]
speciesR
## save
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_random_effects_inc_damage_Species.jpg")
ggsave(speciesR, file = out.filename, dpi=300, height=10, width=12)


## Season plot
seasonR <- random[[2]]
seasonR
## save
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_random_effects_inc_damage_Season.jpg")
ggsave(seasonR, file = out.filename, dpi=300, height=10, width=12)


## plot Q per species

## first predict on model to use in ggplot
rf.data.rescalex <- rf.data.rescale %>% 
  ungroup() %>%
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

out.filename <- paste0(out.dir,"05_mixed_mod_combQ_per_species.jpg")
ggsave(m1, file = out.filename, dpi=300, height=10, width=12)

## plot distance per species
m2 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = DistToSJC002, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Distance to SJC 002") +
  scale_y_continuous(name = "SWP: Mixed Model Fit")
# geom_line(aes(col = Species), linewidth = 1)  


m2

out.filename <- paste0(out.dir,"05_mixed_mod_distance_to_SJC_per_species.jpg")
ggsave(m2, file = out.filename, dpi=300, height=4, width=6)

## plot gage Q per species
m3 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = Q, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Gage Q (MGD)") +
  scale_y_continuous(name = "SWP: Mixed Model Fit")
# geom_line(aes(col = Species), linewidth = 1)  


m3

out.filename <- paste0(out.dir,"05_mixed_mod_gage_Q_per_species.jpg")
ggsave(m3, file = out.filename, dpi=300, height=10, width=12)

## fixed effects
## gage discharge
Q <- plot_model(mod1, type="pred", terms = c("Q_sc"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_Q.jpg")
ggsave(Q, file = out.filename, dpi=305, height=4, width=6)

## combined outfall discharge
of <- plot_model(mod1, type="pred", terms = c("SJC002_POM001Combined_sc"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_combined_outfall.jpg")
ggsave(of, file = out.filename, dpi=305, height=4, width=6)

## replacement

rep <- plot_model(mod1, type="pred", terms = c("Replacement"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_replacement.jpg")
ggsave(rep, file = out.filename, dpi=305, height=4, width=6)
rep
## Year

year <- plot_model(mod1, type="pred", terms = c("Year"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_year.jpg")
ggsave(year, file = out.filename, dpi=305, height=4, width=6)
year
## group
group <- plot_model(mod1, type="pred", terms = c("Group"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_group.jpg")
ggsave(group, file = out.filename, dpi=305, height=4, width=6)
group
## distance

dist <- plot_model(mod1, type="pred", terms = c("DistToSJC002_sc"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_distance.jpg")
ggsave(dist, file = out.filename, dpi=305, height=4, width=6)
dist
## substrate
sub <- plot_model(mod1, type="pred", terms = c("Substrate"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_substrate.jpg")
ggsave(sub, file = out.filename, dpi=305, height=4, width=10)
sub
# position

pos <- plot_model(mod1, type="pred", terms = c("POS"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_POS.jpg")
ggsave(pos, file = out.filename, dpi=305, height=4, width=6)
pos
## damage
dam <- plot_model(mod1, type="pred", terms = c("Damage2"))
out.filename <- paste0(out.dir,"05_mixed_mod_SWP_relationships_damage.jpg")
ggsave(dam, file = out.filename, dpi=305, height=4, width=6)
dam


# Mixed model without damaged trees ---------------------------------------

## rescale variables
rf.data.rescale <- dataall %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         SJC002_sc = rescale(SJC_002, to = c(0,1)),
         SJC002_POM001Combined_sc = rescale(SJC002_POM001Combined, to = c(0,1)),
         POM001_sc = rescale(POM001, to = c(0,1)),
         # RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1))) %>%
  filter(Damage2 == "N")


## all combined - Raw
mod1 <- lmer(formula = SWP ~  Q_sc +  Replacement +  Group + RainfallIntensity +
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + Year +
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
r2_nakagawa(mod1) ## 0.653

set_theme(base = theme_classic(), #To remove the background color and the grids
          theme.font = 'serif',   #To change the font type
          axis.title.size = 1.2,  #To change axis title size
          axis.textsize.x = 1,    #To change x axis text size
          # axis.angle.x = 60,    #To change x axis text angle
          # axis.hjust.x = 1,
          # axis.ticksmar = 50,
          axis.textsize.y = 1)  #To change y axis text size


## first predict on model to use in ggplot
rf.data.rescalex <- rf.data.rescale %>% 
  ungroup() %>%
  mutate(fit.m = predict(mod1, re.form = NA),
         fit.c = predict(mod1, re.form = NULL))

## plot gage Q per species
m3 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = Q, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Gage Q (MGD)") +
  scale_y_continuous(name = "SWP: Mixed Model Fit")
# geom_line(aes(col = Species), linewidth = 1)  


m3

out.filename <- paste0(out.dir,"05_mixed_mod_gage_Q_per_species_no_damaged_trees.jpg")
ggsave(m3, file = out.filename, dpi=300, height=10, width=12)

# Cross Validation --------------------------------------------------------

## split into training and testing

setsize <- floor(nrow(rf.data.val)*0.8)
index <- sample(1:nrow(rf.data.val), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(rf.data.val)

(rf.SWP <- rf.crossValidation(rf, rf.data.val[, 2:19], ## check indexing
                              p=0.10, n=99, ntree=501))

# Fit MSE = 0.8166054 
# Fit percent variance explained = 92.18 
# Median permuted MSE = 0.8840166 
# Median permuted percent variance explained = 91.56 
# Median cross-validation RMSE = 0.931385 
# Median cross-validation MBE = -0.008685373 
# Median cross-validation MAE = 0.6468283 
# Range of ks p-values = 3.573079e-07 0.0008054739 
# Range of ks D statistic = 0.06407563 0.09033613 
# RMSE cross-validation error variance = 0.00128322 
# MBE cross-validation error variance = 0.0008042036 
# MAE cross-validation error variance = 0.0003677655 
