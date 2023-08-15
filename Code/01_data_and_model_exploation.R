## data and model exploration

library(tidyverse)
library(tidylog)
library(lubridate)
# install.packages("MixRF")
# library(MixRF)
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


# upload data -------------------------------------------------------------

load(file = "output_data/00_bio_Q_matched_groups.RData")
head(bioMeanQ_longx)
names(bioMeanQ_longx)


# SWP model without substrate---------------------------------------------------------------

## substrate has NAs, try model without substrate, then with substrate but removing the NAs
## rearrange for ease into RF - remove observed rainfall

rf.data <- bioMeanQ_longx %>%
  dplyr::select(SWP, Species:LifeForm, Month, Year, Season, 
                "DischargeSJC002&POM001Combined(MGD)":"RainFallMod",
                  "POM-001", "TempMeanModF",
                "SJC-002(MGD)":"Q") %>%
  rename(SJC002_POM001Combined = "DischargeSJC002&POM001Combined(MGD)",
         DischargeSJC002_POM001_WN001Combined = "DischargeSJC002,POM001,&WN001Combined(MGD)",
         SJC_002 = "SJC-002(MGD)",
         POM001 = "POM-001",
         WN001 = "WN-001", 
         WN002 = "WN-002(Zone1Ditch)") %>%
    filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  drop_na(SWP) %>%
  select(-Source)

  unique(rf.data$Source)
  rf.data
# Random Forest Model -----------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

## make y data compatible
rf.data.val <- rf.data %>%
  rename(y = SWP) %>% as.data.frame() %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group))

str(rf.data)
# sum(is.na(rf.data)) ## 0
# 
# ind <- which(is.na(rf.data$Substrate))
# rf.data[ind,]

## split into training and testing

setsize <- floor(nrow(rf.data)*0.8)
index <- sample(1:nrow(rf.data), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

names(training)
# 
# trcontrol = trainControl(method='cv', number=10, savePredictions = T,
#                          classProbs = F,summaryFunction = twoClassSummary,returnResamp="all")

rf <- randomForest(y~., data=training, importance = T)
mean(rf$rsq) ## 0.6
varImpPlot(rf)
importance(rf, type = 1)

# PLOT VARIABLE IMPORTANCE
pdf(paste0( "Figures/var_imp_full_model.pdf"), width=25, height=15)

varImpPlot(rf)

dev.off()

# % Var explained: 61.08
# 
# conMat <- confusionMatrix(predict(rf,testing),testing$y)
# 
# save(model, file=paste0(gridFile, "validation.RData"))
# save(conMat, file=paste0(gridFile, "confusion_matrix.RData"))


pdf(paste0(gridFile, "Trees_full_model.pdf"), width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data.val)
plot(full_tree)
text(full_tree)

dev.off()


# Reducing model variables ------------------------------------------------

## non important and remove combined doscharge
names(rf.data)
## month, lifeform 

rf.data.red <- rf.data %>%
  select(-LifeForm, -Month, -DischargeSJC002_POM001_WN001Combined, -SJC002_POM001Combined) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group))

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
mean(rf$rsq) ## 0.62
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


# Remove negative importance variables ------------------------------------
impvar 
imp
names(rf.data)
## WN001, WN002, POM001, RainFallMod, TempMeanModF, 

rf.data.red <- rf.data %>%
  select(-LifeForm, -Month, -SJC002_POM001Combined, -WN002, -WN001, -POM001, -RainFallMod, -TempMeanModF) %>%
  mutate(Species = as.factor(Species), Group = as.factor(Group))

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
mean(rf$rsq) ## 0.65
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
mean(rf$rsq) ## 0.22
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


