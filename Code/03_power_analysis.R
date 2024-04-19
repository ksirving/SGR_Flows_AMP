### power analysis

install.packages("WebPower")
install.packages("pwr")
library(pwr)
library(WebPower)
library(tidyverse)
library(tidylog)
library(rfUtilities)
library(randomForest)

# Find info from SGR data -------------------------------------------------

## upload data
load(file = "output_data/00_bio_Q_matched_groups_distance.RData")
head(bioMeanQ_long_dist)
names(bioMeanQ_long_dist)


rf.data <- bioMeanQ_long_dist %>%
  dplyr::select(SWP, CV, PlantID:LifeForm, Month, Year, Season, POS, Substrate,
                SJC002_POM001Combined:WN002,
                Q, "TempMeanModF",
                RainFallMod:TempMeanModF, DistToSJC002, Source) %>%
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
         Substrate = as.factor(Substrate),
         PlantID = as.factor(PlantID)) %>%
  # mutate(POS = ifelse(Year == 2018, NA, POS)) %>%
  # mutate(POS = recode_factor(POS, blanks = NA)) %>%
  filter(Source %in% c("USGSGauge11087020(MGD)", "LACDPWG44B(MGD)", "LACDPWF313B(MGD)")) %>%
  drop_na(SWP, CV) %>%
  select(-Source)

sum(is.na(rf.data$POS))
dim(rf.data)

str(rf.data)

## define and remove plantids 
plantids <- rf.data %>%
  select(PlantID)

rf.data <- rf.data %>%
  select(-PlantID)

rf.datax <- cbind(plantids, rf.data)
rf.datax

## impute missing values
rf.data.imputed <- rfImpute(SWP ~ ., rf.data)
rf.data.imputed

rf.data.imputedcv <- rfImpute(CV ~ ., rf.data)
rf.data.imputedcv

## correlation of SWP and sjc002

cor(rf.data.imputed$SWP, rf.data.imputed$SJC_002) ## -0.03062607
cor(rf.data.imputedcv$CV, rf.data.imputed$SJC_002) ## 0.0206011

## correlation per year
length(unique(rf.data$PlantID)) ## 113
datacor <- rf.datax %>%
  select(PlantID, SWP, SJC_002, Q, Year, Season) %>%
  group_by(Year, Season) %>%
  summarise(corQ = round(cor(SWP, Q), digits = 4),
            corSJC = round(cor(SWP, SJC_002), digits = 4),
            n = length(SWP))

datacorcv <- rf.data %>%
  select(PlantID, CV, SJC_002, Q, Year, Season) %>%
  group_by(Year, Season) %>%
  summarise(corQ = round(cor(CV, Q),digits=6),
            corSJC = cor(CV, SJC_002),
            n = length(CV))

sum(datacorcv$n)
mean(datacorcv$n)

mean(datacor$corSJC)
mean(datacor$corQ) ## -0.0336666

mean(datacorcv$corSJC)
mean(datacorcv$corQ) ## 0.03059078

# Power analysis ----------------------------------------------------------
example
## Comb Q
?wp.correlation
example=wp.correlation(n=seq(1152,1152*15,1152), r=-0.034, alternative = "two.sided") 
df<-as.data.frame(matrix(ncol = 4, nrow=15))

df
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

## CV
examplecv=wp.correlation(n=seq(1152,1152*15,1152), r=-0.031, alternative = "two.sided") 
dfcv<-as.data.frame(matrix(ncol = 4, nrow=15))

dfcv
colnames(dfcv) <- c("nobvs", "r", "alpha", "Power")
dfcv[,1] <- example$n
dfcv[,2] <- example$r
dfcv[,3] <- example$alpha
dfcv[,4] <- example$power

# x <- rep(c("Spring", "Fall"),times=10)
x1cv <- seq(2023, 2037, 1)
x1cv
# df$Season <- x
dfcv$Year <- as.factor(x1cv)

dfcv$Stressor <- "CV"

alldf <- bind_rows(df,dfcv)

p1 <- ggplot(alldf, aes(x=Year, y = Power, group = Stressor, col = Stressor)) +
  geom_point(aes(group = Stressor, col = Stressor)) +
  facet_wrap(~Stressor)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=0.8, col = "red", linetype="dashed") +
  geom_vline(xintercept="2028", col = "red", linetype="dashed") +
  theme(legend.position = "none")

p1

file.name1 <- "Figures/03_SWP_CV_power_analysis_Q.jpg"
ggsave(p1, filename=file.name1, dpi=600, height=5, width=8)

## SJC - too low correlation to calculate - e.g., 3.70 x 10^-21

example=wp.correlation(n=seq(1050,81900,2100), r=-4.496284e-22, alternative = "less") 
df<-as.data.frame(matrix(ncol = 4, nrow=78))
example
colnames(df) <- c("nobvs", "r", "alpha", "power")
df[,1] <- example$n
df[,2] <- example$r
df[,3] <- example$alpha
df[,4] <- example$power

# x <- rep(c("Spring", "Fall"),times=10)
x1 <- seq(2023,2100,1)
length(x1)
# df$Season <- x
df$Year <- as.factor(x1)


# Power analysis with Cohen's d -------------------------------------------
rf.data.imputed

## began reducing discharge in December 2020
df2020 <- rf.data %>%
  filter(Year == 2020) %>%
  select(PlantID, SWP, CV) %>% distinct()

## most recent full year of data
df2022 <- rf.data %>%
  filter(Year == 2022) %>%
  select(PlantID, SWP, CV) %>% distinct()

length(unique(df2022$PlantID)) ## 186 plants per season per year
# Calculate means
mean1 <- mean(df2020$SWP)
mean2 <- mean(df2022$SWP)

# Calculate standard deviations
sd1 <- sd(df2020$SWP)
sd2 <- sd(df2022$SWP)

# Calculate pooled standard deviation
pooled_sd <- sqrt(((sd1^2 + sd2^2) / 2))

# Calculate Cohen's d
cohen_d <- (mean1 - mean2) / pooled_sd

cat("Cohen's d:", cohen_d, "\n")

# Set parameters for the power analysis
effect_size <- cohen_d   # Desired effect size (Cohen's d)
alpha <- 0.05        # Significance level
power <- 0.80        # Desired power level
n <- NULL            # Sample size (we'll calculate this)

# Perform the power analysis
pwr.t.test(d = effect_size, sig.level = alpha, power = power) # n = 994.571 

## years
1058/186
## 5.6

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

cat("Cohen's d:", cohen_d, "\n")

# Set parameters for the power analysis
effect_size <- cohen_d   # Desired effect size (Cohen's d)
alpha <- 0.05        # Significance level
power <- 0.80        # Desired power level
n <- NULL            # Sample size (we'll calculate this)

# Perform the power analysis
pwr.t.test(d = effect_size, sig.level = alpha, power = power) # n = 3899.839 

## years
3899/186
## 21

