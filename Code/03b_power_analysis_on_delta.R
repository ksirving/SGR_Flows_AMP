## power analysis on delta

### power analysis

# install.packages("WebPower")
# install.packages("pwr")
library(pwr)
library(WebPower)
library(tidyverse)
library(tidylog)
library(rfUtilities)
library(randomForest)


# Upload data -------------------------------------------------------------

load(file = "ignore/06_delta_SWP.Rdata") ## swp
head(AllDeltaDat)

## check the NAs - they happen when a tree, usually a replacement is not present in baseline and current timeframes
# ind <- which(is.na(AllDeltaDat$SWP))
# test <- AllDeltaDat[ind,]

AllDeltaDatswp <- AllDeltaDat %>%
  filter(Type == "Delta", !Group =="G5") %>%
  drop_na(SWP)


load(file = "ignore/06_delta_CV.Rdata") ## swp
head(AllDeltaDat)

AllDeltaDatcv <- AllDeltaDat %>%
  filter(Type == "Delta", !Group =="G5") %>%
  drop_na(CV)


numGroup <- AllDeltaDatswp %>% group_by(Group) %>% summarise(NperGroup = length(unique(PlantID)))
numGroup
sum(numGroup$NperGroup)

numGroup %>% filter(!Group =="G5") %>% summarise(tots = sum(NperGroup))

## 29 trees removed = group 5

# Find info from SGR data -------------------------------------------------


## correlation of SWP and sjc002

cor(AllDeltaDatswp$SWP, AllDeltaDatswp$CombQ) ## -0.33

cor(AllDeltaDatcv$CV, AllDeltaDatcv$CombQ) ##  -0.067


## correlation per year
length(unique(AllDeltaDatswp$PlantID)) ## 67

datacor <- AllDeltaDatswp %>%
  dplyr::select(PlantID, SWP, CombQ, Year, Season) %>%
  distinct() %>%
  group_by(Year, Season) %>%
  summarise(corCombQ = round(cor(SWP, CombQ), digits = 4),
            n = length(PlantID)/6) ## number of trees sampled per season

head(datacor)

sum(datacor$n)
mean(datacor$n)

datacorcv <- AllDeltaDatcv %>%
  dplyr::select(PlantID, CV, CombQ, Year, Season) %>%
  distinct() %>%
  group_by(Year, Season) %>%
  summarise(corCombQ = cor(CV, CombQ),
            n = length(CV)/6) ## sample size?

datacorcv
sum(datacorcv$n)
mean(datacorcv$n)

mean(datacor$corCombQ)
mean(datacor$corQ) ## -0.0336666

mean(datacorcv$corCombQ)
mean(datacorcv$corQ) ## 0.03059078

# Power analysis ----------------------------------------------------------
example
## Q

## define empty df
dfx <- NULL

## pvals to loop through
pvals <- c(0.05,0.1,0.15, 0.2, 0.25, 0.3)

p=1

for (p in 1:length(pvals)) {
  
  example=wp.correlation(n=seq(132,132*15,132), r=-0.33, alpha = pvals[p], alternative = "two.sided") 
  df<-as.data.frame(matrix(ncol = 4, nrow=15))
  
  colnames(df) <- c("nobvs", "r", "alpha", "Power")
  df[,1] <- example$n
  df[,2] <- example$r
  df[,3] <- example$alpha
  df[,4] <- example$power
  
  # x <- rep(c("Spring", "Fall"),times=10)
  x1 <- seq(2023, 2037, 1)
  x1
  # df$Season <- x
  df$Year <- as.factor(x1)
  
  df
  
  df$Stressor <- "SWP"
  
  dfx <- bind_rows(dfx, df)
  
}

head(dfx)

write.csv(dfx, "output_data/03b_pvals_sens_swp_power_analysis")



## CV

## define empty df
dfxCV <- NULL

for (p in 1:length(pvals)) {

examplecv=wp.correlation(n=seq(132,132*15,132), r=-0.067,alpha = pvals[p], alternative = "two.sided") 
dfcv<-as.data.frame(matrix(ncol = 4, nrow=15))


colnames(dfcv) <- c("nobvs", "r", "alpha", "Power")
dfcv[,1] <- examplecv$n
dfcv[,2] <- examplecv$r
dfcv[,3] <- examplecv$alpha
dfcv[,4] <- examplecv$power

# x <- rep(c("Spring", "Fall"),times=10)
x1cv <- seq(2023, 2037, 1)
x1cv
# df$Season <- x
dfcv$Year <- as.factor(x1cv)

dfcv$Stressor <- "CV"

dfxCV <- bind_rows(dfxCV, dfcv)

}

write.csv(dfxCV, "output_data/03b_pvals_sens_CV_power_analysis")


alldf <- bind_rows(dfx,dfxCV)
alldf

alldf <- alldf %>%
  filter(alpha == 0.05)

p1 <- ggplot(alldf, aes(x=Year, y = Power, group = Stressor, col = Stressor)) +
  geom_point(aes(group = Stressor, col = Stressor)) +
  facet_wrap(~Stressor)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=0.8, col = "red", linetype="dashed") +
  # geom_vline(xintercept="2028", col = "red", linetype="dashed") +
  theme(legend.position = "none")

p1

file.name1 <- "Figures/03b_SWP_CV_power_analysis_CombQ_sensitivity.jpg"
ggsave(p1, filename=file.name1, dpi=600, height=5, width=8)

## double the frequencey
examplecv=wp.correlation(n=seq(264,264*15,264), r=-0.067, alternative = "two.sided") 
dfcv<-as.data.frame(matrix(ncol = 4, nrow=15))

dfcv
colnames(dfcv) <- c("nobvs", "r", "alpha", "Power")
dfcv[,1] <- examplecv$n
dfcv[,2] <- examplecv$r
dfcv[,3] <- examplecv$alpha
dfcv[,4] <- examplecv$power

# x <- rep(c("Spring", "Fall"),times=10)
x1cv <- seq(2023, 2037, 1)
x1cv
# df$Season <- x
dfcv$Year <- as.factor(x1cv)

dfcv$Stressor <- "CV"


## SJC - too low correlation to calculate - e.g., 3.70 x 10^-21

# example=wp.correlation(n=seq(1050,81900,2100), r=-4.496284e-22, alternative = "less") 
# df<-as.data.frame(matrix(ncol = 4, nrow=78))
# example
# colnames(df) <- c("nobvs", "r", "alpha", "power")
# df[,1] <- example$n
# df[,2] <- example$r
# df[,3] <- example$alpha
# df[,4] <- example$power
# 
# # x <- rep(c("Spring", "Fall"),times=10)
# x1 <- seq(2023,2100,1)
# length(x1)
# # df$Season <- x
# df$Year <- as.factor(x1)


# Power analysis with Cohen's d -------------------------------------------

## swp data
load(file = "ignore/05_rf_data_imputed_SWP_raw.RData")
rf.datax <- rf.data.imputed

## cv data 
load(file = "ignore/06_rf_data_imputed_CV_raw.RData")

rf.data.imputed <- rf.data.imputed %>%
  dplyr::select(PlantIDNum, CV:Season)

## join together

dataall <- full_join(rf.datax, rf.data.imputed, by = c("PlantIDNum", "Species", "Group", "LifeForm", "Month", "Year"))

## upload plantids identifier
load(file = "ignore/05_plantids.RData")
plantidsD

## join
rf.data <- inner_join(plantidsD, dataall, by = c("PlantIDNum", "Group", "Species")) %>%
  drop_na(SWP, CV) %>%
  distinct()


## began reducing discharge in December 2020
df2020 <- rf.data %>%
  filter(Year == 2020) %>%
  dplyr::select(PlantID, SWP, CV) %>% distinct()

## most recent full year of data
df2022 <- rf.data %>%
  filter(Year == 2022) %>%
  dplyr::select(PlantID, SWP, CV) %>% distinct()

length(unique(df2022$PlantID)) ## 98 plants per season per year
# Calculate means
mean1 <- mean(df2020$SWP)
mean2 <- mean(df2022$SWP)
mean2
# Calculate standard deviations
sd1 <- sd(df2020$SWP)
sd2 <- sd(df2022$SWP)

# Calculate pooled standard deviation
pooled_sd <- sqrt(((sd1^2 + sd2^2) / 2))

# Calculate Cohen's d
cohen_d <- (mean1 - mean2) / pooled_sd

cat("Cohen's d:", cohen_d, "\n") ## 0.1375473

# Set parameters for the power analysis
effect_size <- cohen_d   # Desired effect size (Cohen's d)
alpha <- 0.05        # Significance level
power <- 0.80        # Desired power level
n <- NULL            # Sample size (we'll calculate this)

# Perform the power analysis
pwr.t.test(d = effect_size, sig.level = alpha, power = power) # n = 994.571 

## years
830/98
## 8.469388

## CV
# Calculate means
mean1 <- mean(df2020$CV)
mean2 <- mean(df2022$CV)

# Calculate standard deviations
sd1 <- sd(df2020$CV)
sd2 <- sd(df2022$CV)

# Calculate pooled standard deviation
pooled_sd <- sqrt(((sd1^2 + sd2^2) / 2))

# Calculate Cohen's d
cohen_d <- (mean1 - mean2) / pooled_sd

cat("Cohen's d:", cohen_d, "\n") ## 0.09267874

# Set parameters for the power analysis
effect_size <- cohen_d   # Desired effect size (Cohen's d)
alpha <- 0.05        # Significance level
power <- 0.80        # Desired power level
n <- NULL            # Sample size (we'll calculate this)

# Perform the power analysis
pwr.t.test(d = effect_size, sig.level = alpha, power = power) # n = 3899.839 

## years
1828/98
## 18

