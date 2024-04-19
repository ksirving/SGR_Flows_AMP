## data and model exploration for stem water potential - RELATIVE CHANGE

library(tidyverse)
library(tidylog)
library(lubridate)
# install.packages("MixRF")
library(rfUtilities)
library(randomForest)

# library(kernlab)
# library(rgl)
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

data <- left_join(bioMeanQ_long_distx1, damdata, by = c("PlantID", "Season", "Year"))


# Relative and delta change -----------------------------------------------

## remove columns and make distinct
bioMeanQ_long_distx <- bioMeanQ_long_dist %>%
  group_by(PlantID) %>%
  select(PlantID:Year, SWP , Replacement, Damage) %>%
  distinct()
bioMeanQ_long_distx

## remove damaged grees entitely from data

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
bioMeanQ_long_distRel <-damSeq %>%
  drop_na(SWP) %>%
  filter(!Damage2 == "Y") %>%
  group_by(PlantID) %>%
  mutate(lag.SWP = dplyr::lag(SWP, n = 1, default = NA)) %>% ## lag SWP one season
  #mutate(SWPLag = c(SWP[1], SWP[-793])) %>%
  # mutate(Year = as.numeric(Year)) %>%
  mutate(SWP.relative = SWP/lag.SWP) %>% ## relative SWP 
  mutate(SWP.difference = SWP-lag.SWP) %>%## delta SWP i.e., difference between
  ungroup()
names(bioMeanQ_long_distRel)
## join back to main df

bioMeanQ_long_distJoin <- full_join(bioMeanQ_long_distRel, bioMeanQ_long_dist, by = c("PlantID","Species", "Season", "Group", "LifeForm", "Year", "SWP", "Replacement", "Damage"))

# SWP model with damaged tree variable ---------------------------------------------------------------

## substrate has NAs, try model without substrate, then with substrate but removing the NAs
## rearrange for ease into RF - remove observed rainfall
names(bioMeanQ_long_distJoin)
rf.data <- bioMeanQ_long_distJoin %>%
  dplyr::select(SWP.relative, Species:LifeForm, Month, Year, Season, POS, Substrate,
                RainFallMod:TempMeanModF, Replacement, SJC002_POM001Combined:WN002, ## include damaged trees in all years after damage occurred
                Q, DistToSJC002, Source) %>%
  mutate(POS = as.factor(POS),
         Species = as.factor(Species),
         Group = as.factor(Group),
         LifeForm = as.factor(LifeForm),
         Month = as.factor(Month),
         Substrate = as.factor(Substrate)) %>%
  # mutate(POS = ifelse(Year == 2018, NA, POS)) %>%
  # mutate(POS = recode_factor(POS, blanks = NA)) %>%
  filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  drop_na(SWP.relative) %>%
  select(-Source)

str(rf.data)

str(rf.data)
names(rf.data)

## change blank values in POS
rf.data["POS"][rf.data["POS"]==''] <- NA

rf.data["POS"][rf.data["Year"]=='2018'] <- NA


# Random Forest Model -----------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## impute missing values
rf.data.imputed <- rfImpute(SWP.relative ~ ., rf.data)
head(rf.data.imputed)

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

## make y data compatible
rf.data.val <- rf.data.imputed %>%
  rename(y = SWP.relative) %>% as.data.frame() %>%
  # mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

## run the model
rf <- randomForest(y~., data=rf.data.val, importance = T)
rf
mean(rf$rsq) ## seasonal - 0.70, survey - 0.79 ,no damage = 81.36

## importance quick look
varImpPlot(rf)
importance(rf, type = 1)

# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/01_var_imp_full_model_SWP_Relative.pdf"), width=12, height=8)

varImpPlot(rf, type = 1) ## damage number 

dev.off()

# Remove negative importance variables ------------------------------------
## extract importance
imp <- importance(rf, type = 1)
imp <- as.data.frame(imp)
## get all negative importance vars
vars <- ifelse(imp$`%IncMSE` > 0, rownames(imp), NA)
## remove neg importance vars
vars <- vars[!is.na(vars)]
vars
names(rf.data)
## WN001, WN002, POM001, RainFallMod, TempMeanModF, 

rf.data.red <- rf.data.imputed %>%
  select(SWP.relative, all_of(vars), -LifeForm) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

rf.data.val <- rf.data.red %>%
  rename(y = SWP.relative) %>% as.data.frame()
str(rf.data.red)


## run reduced model
rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ## season - 0.75, survey - 0.85, damage = 85.2
rf


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
                                   Variable == "Damage2" ~ "Damage to Tree"
  ))

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

file.name1 <- "Figures/04_relative_importance_SWP_relative_not_inc_damage.jpg"
ggsave(i1, filename=file.name1, dpi=300, height=5, width=8)


# Partial Plots -----------------------------------------------------------

str(rf.data.val)

jpeg(paste0( "Figures/04_SWP_relative_Season_no_dam.jpg"))
partialPlot(rf, rf.data.val, Season)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Distance_no_dam.jpg"))
partialPlot(rf, rf.data.val, DistToSJC002)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Q_no_dam.jpg"))
partialPlot(rf, rf.data.val, Q)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Species_no_dam.jpg"))
partialPlot(rf, rf.data.val, Species)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_SJC_no_dam.jpg"))
partialPlot(rf, rf.data.val, SJC_002)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Replacement_no_dam.jpg"))
partialPlot(rf, rf.data.val, Replacement)
dev.off()
jpeg(paste0( "Figures/04_SWP_relative_Group_no_dam.jpg"))
partialPlot(rf, rf.data.val, Group)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_POS_no_dam.jpg"))
partialPlot(rf, rf.data.val, POS)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Substrate_no_dam.jpg"))
partialPlot(rf, rf.data.val, Substrate)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Year_no_dam.jpg"))
partialPlot(rf, rf.data.val, Year)
dev.off()

jpeg(paste0( "Figures/04_SWP_relative_Damage_no_dam.jpg"))
partialPlot(rf, rf.data.val, Damage2)
dev.off()

# Cross validation --------------------------------------------------------

## split into training and testing

setsize <- floor(nrow(rf.data.red)*0.8)
index <- sample(1:nrow(rf.data.red), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(rf.data.val)

(rf.cv <- rf.crossValidation(rf, rf.data.val[, 2:13], 
                             p=0.10, n=99, ntree=501) )

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
library(glmmTMB)

names(rf.data.val)

corrf <- rf.data.val[,c(9:14)]
cor(corrf)

## rescale variables
rf.data.rescale <- rf.data.red %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         SJC002_sc = rescale(SJC_002, to = c(0,1)),
         SJC002_POM001Combined_sc = rescale(SJC002_POM001Combined, to = c(0,1)),
         POM001_sc = rescale(POM001, to = c(0,1)),
         # RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1)))


## mixed effects model

## all combined
mod1 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + Year +
               DistToSJC002_sc+ Substrate+ POS+ DischargeSJC002_POM001_WN001Combined_sc +
               # (1|Year) + ## remove due to low ICC
               # (1|Group) + ## remove due to low ICC
               # (1|Year) +
               (1|Season),
             data    = rf.data.rescale) 

summary(mod1)
anova(mod1)
check_singularity(mod1) ## False
icc(mod1, by_group = TRUE)
icc(mod1)
r2_nakagawa(mod1) ## 0.517

## pom and sjc combined
mod2 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + Year +
               DistToSJC002_sc+ Substrate+ POS+ SJC002_POM001Combined_sc + 
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
r2_nakagawa(mod2) ## 0.517

## individual outfall discharge

mod3 <- lmer(formula = SWP.relative ~  Q_sc +  Replacement + Species + Group + Year +
               DistToSJC002_sc+ Substrate+ POS+ SJC002_sc + POM001_sc + 
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
r2_nakagawa(mod3) ## 0.517

### final model plots

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
file.name1 <- paste0(out.dir, "effect_sizes_drivers_of_SWP_relative_not_inc_damage.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

### plot random effects
random <- plot_model(mod1, type="re",
                     vline.color="#A9A9A9", dot.size=1.5,
                     show.values=T, value.offset=.2, show.p=TRUE)
random

out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relative_random_effects_not_inc_damage.jpg")
ggsave(random, file = out.filename, dpi=300, height=10, width=12)

## fixed effects
## gage discharge
Q <- plot_model(mod1, type="pred", terms = c("Q_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_Qnot_inc_damage.jpg")
ggsave(Q, file = out.filename, dpi=300, height=4, width=6)

## combined outfall discharge
of
of <- plot_model(mod1, type="pred", terms = c("SJC002_POM001Combined_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_combined_outfallnot_inc_damage.jpg")
ggsave(of, file = out.filename, dpi=300, height=4, width=6)

## replacement

rep <- plot_model(mod1, type="pred", terms = c("Replacement"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_replacementnot_inc_damage.jpg")
ggsave(rep, file = out.filename, dpi=300, height=4, width=6)

## Year
year
year <- plot_model(mod1, type="pred", terms = c("Year"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_yearnot_inc_damage.jpg")
ggsave(year, file = out.filename, dpi=300, height=4, width=6)

## group
group <- plot_model(mod1, type="pred", terms = c("Group"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_groupnot_inc_damage.jpg")
ggsave(group, file = out.filename, dpi=300, height=4, width=6)

## distance
dist
dist <- plot_model(mod1, type="pred", terms = c("DistToSJC002_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_distancenot_inc_damage.jpg")
ggsave(dist, file = out.filename, dpi=300, height=4, width=6)

## substrate
sub <- plot_model(mod1, type="pred", terms = c("Substrate"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_substratenot_inc_damage.jpg")
ggsave(sub, file = out.filename, dpi=300, height=4, width=10)

# position
pos
pos <- plot_model(mod1, type="pred", terms = c("POS"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_POSnot_inc_damage.jpg")
ggsave(pos, file = out.filename, dpi=300, height=4, width=6)

## damage
dam <- plot_model(mod1, type="pred", terms = c("Damage2"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_damagenot_inc_damage.jpg")
ggsave(dam, file = out.filename, dpi=300, height=4, width=6)

