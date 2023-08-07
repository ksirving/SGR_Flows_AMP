## data and model exploration

library(tidyverse)
library(tidylog)
library(lubridate)
install.packages("MixRF")
library(MixRF)
library(randomForest)

# library(kernlab)
library(rgl)
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

load(file = "output_data/00_bio_meanQ_data.RData")
head(bioMeanQ)
names(bioMeanQ)

## rearrange for ease into RF

rf.data <- bioMeanQ %>%
  dplyr::select(SWP,CV, Species:LifeForm, Month, Year, Season, Substrate, "DischargeSJC002&POM001Combined(MGD)":"LACDPWF313B(cfs)") %>%
          drop_na(SWP) 


# Random Forest Model -----------------------------------------------------

set.seed(234) ## reproducibility

# NUMBER OF BOOTSTRAP REPLICATES
b=10001

sp=0.6 ## for dependence plots

## get path for functions
source("/Users/katieirving/Library/CloudStorage/OneDrive-SCCWRP/Documents - Katie’s MacBook Pro/git/RB9_Vulnerability_Arroyo_Toad/original_model/Current/randomForests/PARTITIONING/DATA3/Functions.R")

## make y data compatible
rf.data.val <- rf.data %>%
  mutate(y = SWP)

## split into training and testing

setsize <- floor(nrow(rf.data)*0.8)
index <- sample(1:nrow(rf.data), size = setsize)
training <- rf.data.val[index,]
testing <- rf.data.val[-index,]

trcontrol = trainControl(method='cv', number=10, savePredictions = T,
                         classProbs = TRUE,summaryFunction = twoClassSummary,returnResamp="all")

model = train(y ~ . , data=training, method = "rf", trControl = trcontrol,metric="ROC", na.action = na.roughfix) 

conMat <- confusionMatrix(predict(model,testing),testing$y)

save(model, file=paste0(gridFile, "validation.RData"))
save(conMat, file=paste0(gridFile, "confusion_matrix.RData"))


pdf(paste0(gridFile, "Trees_no_clim_gridded.pdf"), width=25, height=15)

full_tree <- rpart(y~., method = "class", control = rpart.control(cp = 0, minsplit = 2), data = rf.data)
plot(full_tree)
text(full_tree, use.n = T)

dev.off()


# Coeficients and importance ----------------------------------------------

# PLOT VARIABLE IMPORTANCE
pdf(paste0(gridFile, "var_imp_no_clim_gridded.pdf"), width=25, height=15)

varImpPlot(rf.final)

dev.off()
