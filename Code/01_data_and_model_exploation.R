## data and model exploration for stem water potential

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
# 
# data(sleepstudy)
# head(sleepstudy)
# 
# tmp = MixRF(Y = sleepstudy$Reaction, X = as.data.frame(sleepstudy$Days), 
#             random = "(Days|Subject)", data = sleepstudy, initialRandomEffects = 0, 
#             ErrorTolerance = 0.01, MaxIterations = 100)
# tmp

## output file for figures
out.dir <- "/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/SGR_Flows_AMP/Figures/"


# upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
head(bioMeanQ_long_dist)
names(bioMeanQ_long_dist)


# SWP model without substrate---------------------------------------------------------------

## substrate has NAs, try model without substrate, then with substrate but removing the NAs
## rearrange for ease into RF - remove observed rainfall
names(rf.data)
rf.data <- bioMeanQ_long_dist %>%
  dplyr::select(SWP, Species:LifeForm, Month, Year, Season, POS, Substrate,
                RainFallMod:Replacement, SJC002_POM001Combined:WN002,
                Q, DistToSJC002, Source) %>%
  # rename(SJC002_POM001Combined = "DischargeSJC002&POM001Combined(MGD)",
  #        DischargeSJC002_POM001_WN001Combined = "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #        SJC_002 = "SJC-002(MGD)",
  #        POM001 = "POM-001",
  #        WN001 = "WN-001", 
  #        WN002 = "WN-002(Zone1Ditch)") %>%
  mutate(POS = as.factor(POS),
         Species = as.factor(Species),
         Group = as.factor(Group),
         LifeForm = as.factor(LifeForm),
         Month = as.factor(Month),
         Substrate = as.factor(Substrate)) %>%
  # mutate(POS = ifelse(Year == 2018, NA, POS)) %>%
  # mutate(POS = recode_factor(POS, blanks = NA)) %>%
  filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  drop_na(SWP) %>%
  select(-Source)

str(rf.data)

## change blank values in POS
rf.data["POS"][rf.data["POS"]==''] <- NA

rf.data["POS"][rf.data["Year"]=='2018'] <- NA

sum(is.na(rf.data$POS))
dim(rf.data)



# Random Forest Model -----------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## impute missing values
rf.data.imputed <- rfImpute(SWP ~ ., rf.data)
head(rf.data.imputed)

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

## make y data compatible
rf.data.val <- rf.data.imputed %>%
  rename(y = SWP) %>% as.data.frame() %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

# str(rf.data.val)
# str(rf.data)
# levels(rf.data.val$Year)
# unique(rf.data.val$Year)

# sum(is.na(rf.data)) ## 0
# 
# ind <- which(is.na(rf.data$Substrate))
# rf.data[ind,]
# 
# trcontrol = trainControl(method='cv', number=10, savePredictions = T,
#                          classProbs = F,summaryFunction = twoClassSummary,returnResamp="all")

rf <- randomForest(y~., data=rf.data.val, importance = T)
mean(rf$rsq) ## 0.85
varImpPlot(rf)
importance(rf, type = 1)

# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/01_var_imp_full_model_SWP.pdf"), width=12, height=8)

varImpPlot(rf, type = 1)

dev.off()

# % Var explained: 61.08
# 
# conMat <- confusionMatrix(predict(rf,testing),testing$y)
# 
# save(model, file=paste0(gridFile, "validation.RData"))
# save(conMat, file=paste0(gridFile, "confusion_matrix.RData"))


# pdf(paste0(gridFile, "Trees_full_model.pdf"), width=25, height=15)
# 
# full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data.val)
# plot(full_tree)
# text(full_tree)
# 
# dev.off()


# Reducing model variables ------------------------------------------------

## non important and remove combined doscharge
names(rf.data)
## put importance in a df
imp <- as.data.frame(importance(rf, type = 1))

imp
## NA all vars with negative importance 
vars <- ifelse(imp$`%IncMSE` > 0, rownames(imp), NA)
## remove Nas
vars <- vars[!is.na(vars)]
vars ## vars with posituve importance

## month, lifeform 
rf.data.red <- rf.data.imputed %>%
  select(SWP, all_of(vars), -LifeForm) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group), Year = as.factor(Year)) %>% 
  droplevels()

rf.data.val <- rf.data.red %>%
  rename(y = SWP) %>% as.data.frame()
str(rf.data.red)

## split into training and testing

setsize <- floor(nrow(rf.data.red)*0.8)
index <- sample(1:nrow(rf.data.red), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(training)
# 
# trcontrol = trainControl(method='cv', number=10, savePredictions = T,
#                          classProbs = F,summaryFunction = twoClassSummary,returnResamp="all")

rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ## 0.90
rf
varImpPlot(rf)
importance(rf, type = 1)

# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/var_imp_reduced_model.pdf"), width=25, height=15)

varImpPlot(rf)

dev.off()

## partial plots
# 
# data(iris)
# head(iris)
# set.seed(543)
# iris.rf <- randomForest(Species~., iris)
# partialPlot(iris.rf, iris, Petal.Width, "versicolor")
str(rf.data.val)

partialPlot(rf, rf.data.val, Season)
partialPlot(rf, rf.data.val, Year)
partialPlot(rf, rf.data.val, Q)
# partialPlot(rf, rf.data.val, LACDPWG44B)
partialPlot(rf, rf.data.val, POM001)
partialPlot(rf, rf.data.val, RainFallMod)
partialPlot(rf, rf.data.val, SJC_002)
partialPlot(rf, rf.data.val, TempMeanModF)
# partialPlot(rf, rf.data.val, USGSGauge11087020)
partialPlot(rf, rf.data.val, WN002)
partialPlot(rf, rf.data.val, Replacement)
partialPlot(rf, rf.data.val, WN001)

partialPlot(rf, rf.data.val, Species)
partialPlot(rf, rf.data.val, Group)

imp <- importance(rf, type = 1)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op <- par(mfrow=c(2, 3))

for (i in seq_along(impvar)) {
  partialPlot(rf, training, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]))
}
par(op)

# Variable importance -----------------------------------------------------

library("scales")

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
importance(rf, type = 1)
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
                                   Variable == "WN002" ~ "Discharge: WN",
                                   Variable == "RainfallIntensity" ~ "Observed Rainfall",
                                   Variable == "POM001" ~ "Discharge: POM",
                                   Variable == "WN002" ~ "Discharge: WN (zone ditch)",
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

file.name1 <- "Figures/01_relative_importance_SWP.jpg"
ggsave(i1, filename=file.name1, dpi=300, height=5, width=8)


# PLOT VARIABLE IMPORTANCE
# pdf(paste0( "Figures/var_imp_reduced_model_SWP_2.pdf"), width=12, height=8)
# 
# varImpPlot(rf, type = 1, main = "")
# 
# dev.off()


# Partial Plots -----------------------------------------------------------



str(rf.data.val)

jpeg(paste0( "Figures/02_SWP_Season.jpg"))
partialPlot(rf, rf.data.val, Season)
dev.off()

jpeg(paste0( "Figures/02_SWP_Distance.jpg"))
partialPlot(rf, rf.data.val, DistToSJC002)
dev.off()

jpeg(paste0( "Figures/02_SWP_Q.jpg"))
partialPlot(rf, rf.data.val, Q)
dev.off()
jpeg(paste0( "Figures/02_SWP_Species.jpg"))
partialPlot(rf, rf.data.val, Species)
dev.off()
jpeg(paste0( "Figures/02_SWP_SJC.jpg"))
partialPlot(rf, rf.data.val, SJC_002)
dev.off()
jpeg(paste0( "Figures/02_SWP_Replacement.jpg"))
partialPlot(rf, rf.data.val, Replacement)
dev.off()
jpeg(paste0( "Figures/02_SWP_Group.jpg"))
partialPlot(rf, rf.data.val, Group)
dev.off()

jpeg(paste0( "Figures/02_SWP_POS.jpg"))
partialPlot(rf, rf.data.val, POS)
dev.off()

jpeg(paste0( "Figures/02_SWP_Substrate.jpg"))
partialPlot(rf, rf.data.val, Substrate)
dev.off()



# partialPlot(rf, rf.data.val, LACDPWG44B)
# partialPlot(rf, rf.data.val, POM001)
# partialPlot(rf, rf.data.val, RainFallMod)

# partialPlot(rf, rf.data.val, TempMeanModF)
# partialPlot(rf, rf.data.val, USGSGauge11087020)
# partialPlot(rf, rf.data.val, WN002)
# partialPlot(rf, rf.data.val, WN001)

# imp <- importance(rf, type = 1)
# impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
# op <- par(mfrow=c(2, 3))
# 
# for (i in seq_along(impvar)) {
#   partialPlot(rf, training, impvar[i], xlab=impvar[i],
#               main=paste("Partial Dependence on", impvar[i]))
# }
# par(op)


# Remove species ----------------------------------------------------------

impvar 
imp
names(rf.data)
## WN001, WN002, POM001, RainFallMod, TempMeanModF, 

rf.data.red <- rf.data %>%
  select(-LifeForm, -Month, -SJC002_POM001Combined, -WN002, -WN001, -POM001, -RainFallMod, -TempMeanModF, -Species, -Season) %>%
  mutate( Group = as.factor(Group))

rf.data.val <- rf.data.red %>%
  rename(y = SWP) %>% as.data.frame()
str(rf.data.red)
# sum(is.na(rf.data)) ## 0
# 
# ind <- which(is.na(rf.data$Substrate))
# rf.data[ind,]

## split into training and testing

setsize <- floor(nrow(rf.data.red)*0.8)
index <- sample(1:nrow(rf.data.red), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(training)
# 
# trcontrol = trainControl(method='cv', number=10, savePredictions = T,
#                          classProbs = F,summaryFunction = twoClassSummary,returnResamp="all")

rf <- randomForest(y~., data=rf.data.val, importance = T, ntree = 1000)
mean(rf$rsq) ##
rf
varImpPlot(rf)
importance(rf, type = 1)

# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/var_imp_reduced_model_important_vars_no_species_no_season.pdf"), width=25, height=15)

varImpPlot(rf)

dev.off()

## partial plots
# 
# data(iris)
# head(iris)
# set.seed(543)
# iris.rf <- randomForest(Species~., iris)
# partialPlot(iris.rf, iris, Petal.Width, "versicolor")
str(rf.data.val)

partialPlot(rf, rf.data.val, Season)
partialPlot(rf, rf.data.val, Year)
partialPlot(rf, rf.data.val, Q)
# partialPlot(rf, rf.data.val, LACDPWG44B)
# partialPlot(rf, rf.data.val, POM001)
# partialPlot(rf, rf.data.val, RainFallMod)
partialPlot(rf, rf.data.val, SJC_002)
# partialPlot(rf, rf.data.val, TempMeanModF)
# partialPlot(rf, rf.data.val, USGSGauge11087020)
# partialPlot(rf, rf.data.val, WN002)
partialPlot(rf, rf.data.val, Replacement)
# partialPlot(rf, rf.data.val, WN001)

partialPlot(rf, rf.data.val, Species)
partialPlot(rf, rf.data.val, Group)

imp <- importance(rf, type = 1)
impvar <- rownames(imp)[order(imp[, 1], decreasing=TRUE)]
op <- par(mfrow=c(2, 3))

for (i in seq_along(impvar)) {
  partialPlot(rf, training, impvar[i], xlab=impvar[i],
              main=paste("Partial Dependence on", impvar[i]))
}
par(op)


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

## Feb 2024
# Fit MSE = 1.043827 
# Fit percent variance explained = 90.15 
# Median permuted MSE = 2.622232 
# Median permuted percent variance explained = 74.9 
# Median cross-validation RMSE = 1.614194 
# Median cross-validation MBE = 0.002431555 
# Median cross-validation MAE = 1.13129 
# Range of ks p-values = 1.332268e-15 2.747508e-05 
# Range of ks D statistic = 0.07668067 0.1355042 
# RMSE cross-validation error variance = 0.003508198 
# MBE cross-validation error variance = 0.002995246 
# MAE cross-validation error variance = 0.001178508

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

## data
# load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
# head(bioMeanQ_long_dist)
# names(bioMeanQ_longx)
# 
# rf.data <- bioMeanQ_long_dist %>%
#   dplyr::select(SWP, Species:LifeForm, Month, Year, Season, 
#                 "DischargeSJC002&POM001Combined(MGD)":"RainFallMod",
#                 "POM-001", "TempMeanModF",
#                 "SJC-002(MGD)":"Q", DistToSJC002) %>%
#   rename(SJC002_POM001Combined = "DischargeSJC002&POM001Combined(MGD)",
#          DischargeSJC002_POM001_WN001Combined = "DischargeSJC002,POM001,&WN001Combined(MGD)",
#          SJC_002 = "SJC-002(MGD)",
#          POM001 = "POM-001",
#          WN001 = "WN-001", 
#          WN002 = "WN-002(Zone1Ditch)") %>%
#   filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
#   drop_na(SWP, DistToSJC002) %>%
#   select(-Source) %>%
#   mutate(Year = as.factor(Year),
#           Group = as.factor(Group)) %>% 
#   droplevels() 
# 
# head(rf.data)
names(rf.data.red)
str(rf.data.red)

corrf <- rf.data.red[,c(8,10:16)]
str(corrf)
cor(corrf)

## rescale variables
rf.data.rescale <- rf.data.red %>%
  mutate(Q_sc = rescale(Q, to = c(0,1)),
         DischargeSJC002_POM001_WN001Combined_sc = rescale(DischargeSJC002_POM001_WN001Combined, to = c(0,1)),
         RainfallIntensity_sc = rescale(RainfallIntensity, to = c(0,1)),
         DistToSJC002_sc = rescale(DistToSJC002, to = c(0,1)))

## variables from RF model to include
?rescale
## removed  DischargeSJC002_POM001_WN001Combined
#  

mod1 <- lmer(formula = SWP ~  RainfallIntensity_sc + Q_sc + DischargeSJC002_POM001_WN001Combined_sc + Replacement + Year + Group + DistToSJC002_sc+ Substrate + POS+ 
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
r2_nakagawa(mod1) ## 0.63

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
file.name1 <- paste0(out.dir, "effect_sizes_drivers_of_SWP.jpg")
ggsave(ests, filename=file.name1, dpi=300, height=8, width=10)

### plot random effects
random <- plot_model(mod1, type="re",
           vline.color="#A9A9A9", dot.size=1.5,
           show.values=T, value.offset=.2, show.p=TRUE)

r1 <- random[[1]]

r1

# out.filename <- paste0(out.dir,"00_mixed_mod_SWP_random_effects_species.jpg")
# ggsave(r1, file = out.filename, dpi=300, height=4, width=6)

r2 <- random[[2]]

r2

r3 <- plot_grid(list(r2, r1))

out.filename <- paste0(out.dir,"00_mixed_mod_SWP_random_effects.jpg")
ggsave(r3, file = out.filename, dpi=300, height=10, width=12)

### plot realtionships
# results <- plot_model(mod1, type="pred",
#            vline.color="#A9A9A9", dot.size=1.5,
#            show.values=T, value.offset=.2, show.p=TRUE, ci.style = c("whisker"))
# 
mod1
## gage discharge
Q <- plot_model(mod1, type="pred", terms = c("Q_sc"))
out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_Q.jpg")
ggsave(Q, file = out.filename, dpi=300, height=4, width=6)

## combined outfall discharge
of
of <- plot_model(mod1, type="pred", terms = c("DischargeSJC002_POM001_WN001Combined_sc"))
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

# allFigs <- plot_grid(list(results[[3]],results[[4]], results[[5]],
#                           results[[6]], results[[7]], results[[8]]))
# 
# out.filename <- paste0(out.dir,"00_mixed_mod_SWP_relationships_ALL.jpg")
# ggsave(allFigs, file = out.filename, dpi=300, height=15, width=21)
