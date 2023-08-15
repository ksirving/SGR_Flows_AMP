## format raw data

library(tidyverse)
library(tidylog)
library(lubridate)



# Upload data -------------------------------------------------------------

data <- read.csv("input_data/RawData_SGR.csv") %>%
  rename(Substrate = Susbstrate)
head(data)


### make substrate values consistent and create separate df as no season year

unique(data$Substrate)

data <- data %>%
  mutate(Substrate = case_when(Substrate == "C,F" ~ "Cobble/Fines",
                   Substrate == "F" ~ "Fines",
                   Substrate == "C,F,G" ~ "Cobble/Fines/Gravel",
                   Substrate == "C" ~ "Cobble",
                   Substrate == "F,G" ~ "Fines/Gravel",
                   Substrate == "B,C,F" ~ "Bedrock/Cobble/Fines",
                   Substrate == "F,C" ~ "Fines/Cobble",
                   Substrate == "B,F" ~ "Bedrock/Fines",
                   Substrate == "Fines/Grave;" ~ "Fines/Gravel",
                   Substrate == "Cobble/Fines" ~ "Cobble/Fines",
                   Substrate == "Cobble" ~ "Cobble",
                   Substrate == "Fines/cobble" ~ "Fines/Cobble",
                   Substrate == "Fines" ~ "Fines"))

subs <- data %>%
  select(PlantID:LifeForm, Substrate)
  


### make longer to arrange variables, create variable, season and year columns, make variables columns
names(data)
str(data)

data_long <- data %>%
  select(-Substrate) %>%  ## remove substrate
  mutate(across(where(is.double), as.character), across(where(is.integer), as.character)) %>%
  pivot_longer(SWP_F18:Recruits_S21, names_to = "Variable", values_to = "Values") %>%
  separate(Variable, into = c("VariableName", "SeasonYear")) %>% ## warning message ok: refers to substrate as has no season/year
  separate(SeasonYear, 
           into = c("Season", "Year"), 
           sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  mutate(Season = case_when(Season == "F" ~ "Fall",
                              Season == "S" ~ "Spring")) %>%
  pivot_wider(names_from = VariableName, values_from = Values)

## add substrate back

alldata <- data_long %>%
  full_join(subs, by = c("PlantID", 'Species', "Group", "LifeForm")) %>%
  mutate(Year = as.numeric(Year)) %>%
  mutate(Year = Year+2000)

head(alldata)
str(alldata)

# Plot some data ----------------------------------------------------------

## prep data 

## season and year as factor, swp, cv as numeric

alldata <- alldata %>%
  mutate(SWP = as.numeric(SWP),
         CV = as.numeric(CV),
         Season = as.factor(Season),
         Year = as.factor(Year)) %>%
  pivot_longer(SWP:CV, names_to = "Variable", values_to = "Values")

write.csv(alldata, "output_data/00_bio_data_long.csv")


p1 <- ggplot(alldata, aes(x=Year, y=Values)) + 
  geom_boxplot() +
  facet_wrap(~Variable, scales= "free_y")
  

p1

p2 <- ggplot(alldata, aes(x=Season, y=Values)) + 
  geom_boxplot() +
  facet_wrap(~Variable, scales= "free_y")

p2

p3 <- ggplot(alldata, aes(x=Year, y=Values)) + 
  geom_boxplot() +
  facet_grid(vars(Season), vars(Variable), scales= "free")

p3


# Discharge ---------------------------------------------------------------

hydro <- read.csv("input_data/Daily_Data_Dec2022.csv", header = T) %>% as_tibble()
## get row 2 as headings
colnames(hydro) <- unlist(hydro[row.names(hydro)==2,])
dim(hydro)

hydro <- hydro[, 1:14] ## remove weird/not needed columns

## remove row 2
hydro <- hydro[!row.names(hydro) ==2,]

## also remove row 1 as this is an average
hydro <- hydro[!row.names(hydro) ==1,]
hydro
### make time series of discharge

## format data
## format date, make long, make unit and gauge columns, add impacting group
names(hydro)


hydrox <- hydro %>%
  mutate(Date = lubridate::mdy(Date)) %>%
  pivot_longer(cols = "USGS Avg Flow (cfs)":"Discharge SJC002, POM001, & WN001 Combined (MGD)", 
               names_to = "Source", values_to = "Q") %>%
  mutate("Source" = gsub( " ", "", Source)) #%>%

cfs <- unique(hydrox$Source)[c(1,3,5)] ## get cfs q
cfs
unique(hydrox$Source)

hydrox <- hydrox %>%
  mutate(Units = ifelse(Source %in% cfs, "CFS", "MGD")) %>%
  # mutate(Units = ifelse(Source == "RainfallIntensity(in)", "INCHES", Units)) %>%
  mutate(Q = as.numeric(Q)) %>%
  mutate(Group = ifelse(Source %in% c("USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Groups 1,2,3", "Group 4")) %>%
  mutate(Group = ifelse(Source %in% c("LACDPWF313B(MGD)", "LACDPWF313B(cfs)"), "Groups 5", Group)) %>%
  mutate(Group = ifelse(Source %in% c( "POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                      "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                      "SJC-002(MGD)"), "All Groups ", Group)) %>%
  rename(RainfallIntensity = "Rainfall Intensity (in)") %>%
  mutate(RainfallIntensity = as.numeric(RainfallIntensity))

write.csv(hydrox, "output_data/00_daily_Q_long.csv")
unique(hydrox$Source)


# Plot discharge ----------------------------------------------------------

str(hydrox)

q1 <- ggplot(subset(hydrox, Group == "Group 4"), aes(x=Date, y=Q), group = Source, col = Source) + 
  geom_line(aes(group = Source, col = Source)) +
  # scale_x_date(date_breaks = "1 month", date_labels = "%M")+
  facet_wrap(~Units, scales= "free_y")


q1
str(hydrox)
coeffcfs <- 5000/3
coeff

## plot rainfall with CFS gages
cfsSources <- -c("USGSAvgFlow(cfs)","LACDPWF313B(cfs)", "LACDPWG44B(cfs)" )


q2 <- ggplot(subset(hydrox, Source %in% cfsSources), aes(x=Date, y=Q)) + 
  geom_line(aes(y=RainfallIntensity*coeffcfs), col = "lightblue") +
  geom_line(aes(y=Q), col = "blue") +
  scale_y_continuous(
    # Features of the first axis
    name = "Discharge (cfs)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeff, name="Rainfall (in)")) +
  facet_grid(~Source, scales= "free_y")




q2


## plot rainfall with MGD 
mgdSources <- unique(hydrox$Source)[-c(1,3,5)]

coeffmgd <- 3000/3
coeffmgd

q3 <- ggplot(subset(hydrox, Source %in% mgdSources), aes(x=Date, y=Q)) + 
  geom_line(aes(y=RainfallIntensity*coeffmgd), col = "lightblue") +
  geom_line(aes(y=Q), col = "blue") +
  scale_y_continuous(
    # Features of the first axis
    name = "Discharge (MGD)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeffmgd, name="Rainfall (in)")) +
  facet_grid(~Source, scales= "free_y")

q3

## plot only effluent

effsources <- unique(hydrox$Source)[c(7:12)]
effsources

coeffeff <- 60/3
coeffeff

q4 <- ggplot(subset(hydrox, Source %in% effsources), aes(x=Date, y=Q)) + 
  geom_line(aes(y=RainfallIntensity*coeffeff), col = "lightblue") +
  geom_line(aes(y=Q), col = "blue") +
  scale_y_continuous(
    # Features of the first axis
    name = "Discharge (MGD)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./coeffeff, name="Rainfall (in)")) +
  facet_grid(~Source, scales= "free_y")

q4



# Climate data from PRISM -------------------------------------------------
head(hydrox)
## https://prism.oregonstate.edu/explorer/

clims <- read.csv("ignore/PRISM_ppt_tmean_stable_4km_20140101_20230101_34.0283_-118.0331.csv") %>%
  rename(RainFallMod = 2, TempMeanModF = 3)

head(clims)

## format date same as hydro

climsx <- clims %>%
  mutate(Date = lubridate::mdy(Date))

head(climsx)

## join with hydro

hydroxclim <- inner_join(climsx, hydrox, by="Date")
head(hydroxclim)

### plot
## both rain

r1 <- ggplot(hydroxclim, aes(x=Date, y=RainFallMod)) + 
  geom_line(aes(y=RainfallIntensity), col = "lightblue") +
  geom_line(aes(y = RainFallMod), col = "blue") +
  scale_y_continuous(
    # Features of the first axis
    name = "Modelled Rainfall (Inches)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~., name="Observed Rainfall (Inches)")) 
r1


# Creating Monthly data ----------------------------------------------------

names(hydroxclim)

## split date and make rainfall in with source
hydroxDate <- hydroxclim %>%
  select(-Group,-Units) %>%
  pivot_wider(names_from = Source, values_from = Q) %>%
  pivot_longer(cols = c("RainFallMod":"DischargeSJC002,POM001,&WN001Combined(MGD)"), names_to = "Source", values_to = "Q") %>%
  mutate(Q = as.numeric(Q)) %>%
  separate(Date, into = c("Year", "Month", "Day"))

## monthly summaries and add groups and units
### need repeated values for each group

unique(hydroxDate$Source)

# test <- hydroxDate %>%
#   filter(Source == "RainfallIntensity",
#          Year == 2018)
# 
# test2 <- test %>%
#   group_by(Year, Month) %>%
#   summarise(MinQ = mean(na.omit(Q)))

# unique(test$Year)

# monthlyQ <- hydroxDate %>%
#   group_by(Year, Month, Source) %>%
#   summarise(MinQ = min(na.omit(Q)), MaxQ = max(na.omit(Q)), MedianQ = median(na.omit(Q)), MeanQ = mean(na.omit(Q))) %>%
#   mutate(Units = ifelse(Source %in% cfs, "CFS", "MGD")) %>%
#   mutate(Units = ifelse(Source == "RainfallIntensity", "INCHES", Units)) %>%
#   mutate(Group1 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
#                                         "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
#                                         "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
#   mutate(Group2 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
#                                         "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
#                                         "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
#   mutate(Group3 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
#                                         "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
#                                         "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
#   mutate(Group4 = ifelse(Source %in% c("RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
#                                        "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
#                                        "SJC-002(MGD)","LACDPWG44B(MGD)", "LACDPWG44B(cfs)"), "Yes", NA)) %>%
#   mutate(Group5 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
#                                        "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
#                                        "SJC-002(MGD)", "LACDPWF313B(MGD)", "LACDPWF313B(cfs)"), "Yes", NA)) %>%
#   pivot_longer(Group1:Group5, names_to = "Group", values_to = "Value") %>%
#   mutate(GroupCheck = ifelse(Value =="Yes", Group, NA)) %>%
#   drop_na(GroupCheck) %>% select(-GroupCheck, -Value)   %>% 
#   mutate(GroupQ = gsub("roup", "", Group)) %>%
#   mutate(Year = as.factor(Year))

cfs <- c( "USGSAvgFlow(cfs)", "LACDPWF313B(cfs)", "LACDPWG44B(cfs)")
inches <- c("RainFallMod", "RainfallIntensity")

monthlyQ <- hydroxDate %>%
  group_by(Year, Month, Source) %>%
  summarise(MinQ = min(na.omit(Q)), MaxQ = max(na.omit(Q)), MedianQ = median(na.omit(Q)), MeanQ = mean(na.omit(Q))) %>%
  mutate(Units = ifelse(Source %in% cfs, "CFS", "MGD")) %>%
  mutate(Units = ifelse(Source %in% inches, "INCHES", Units)) %>%
  mutate(Units = ifelse(Source == "TempMeanModF", "Fahrenheit", Units)) %>%
  # mutate(Group1 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
  #                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #                                       "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  # mutate(Group2 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
  #                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #                                       "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  # mutate(Group3 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
  #                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #                                       "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  # mutate(Group4 = ifelse(Source %in% c("RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
  #                                      "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #                                      "SJC-002(MGD)","LACDPWG44B(MGD)", "LACDPWG44B(cfs)"), "Yes", NA)) %>%
  # mutate(Group5 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
  #                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
  #                                       "SJC-002(MGD)", "LACDPWF313B(MGD)", "LACDPWF313B(cfs)"), "Yes", NA)) %>%
  # pivot_longer(Group1:Group5, names_to = "Group", values_to = "Value") %>%
  # mutate(GroupCheck = ifelse(Value =="Yes", Group, NA)) %>%
  # drop_na(GroupCheck) %>% select(-GroupCheck, -Value)   %>% 
  # mutate(GroupQ = gsub("roup", "", Group)) %>%
  mutate(Year = as.factor(Year))

head(monthlyQ)


# Join with Bio data ------------------------------------------------------

head(monthlyQ)
head(alldata)

sum(is.na(alldata$Variable))

## use inner join for now, can extend to full join after model test, and include lags etc
bioMeanQ <- inner_join(alldata, monthlyQ, by = c("Year"), relationship = "many-to-many") %>%
  drop_na(Values) %>% ### NA values in response, remove - due to replacements of trees
  pivot_wider(names_from = Variable, values_from = Values) %>%
  select(-c(MinQ, MaxQ, MedianQ, Units)) %>%
  pivot_wider(names_from = Source, values_from = MeanQ) %>%
  mutate(Replacement = grepl("R", PlantID)) %>%
  mutate(Replacement = as.factor(Replacement), PlantID = as.factor(PlantID))
names(bioMeanQ)
str(bioMeanQ)

## how many tree replacements
bioMeanQ %>% select(PlantID, Replacement) %>% distinct() %>% group_by(Replacement) %>% tally()
# 1 FALSE          96
# 2 TRUE           17


## rainfall has lots of NAs!!! - happens when no values (ie all NAs) for one month, maybe find an alternative? 
## using modelled data from PRISM

save(bioMeanQ, file = "output_data/00_bio_meanQ_data.RData")

cor(bioMeanQ$RainFallMod, bioMeanQ$RainfallIntensity, use = "complete.obs")
## 0.96


# Related Q to Group ------------------------------------------------------

## one column of Q from gages - related to each group
## groups 1,2,3 = USGS Gauge 11087020 
##. group 4 = LACDPW G44B 
## group 5 = LACDPW F313B 
head(bioMeanQ)

bioMeanQ_long <- bioMeanQ %>%
  pivot_longer(c("LACDPWF313B(MGD)":"LACDPWG44B(cfs)", "USGSAvgFlow(cfs)":"USGSGauge11087020(MGD)"), 
               names_to = "Source", values_to = "Q") %>%
  mutate(Group1 = ifelse(Source %in% c("USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
                mutate(Group2 = ifelse(Source %in% c("USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
                mutate(Group3 = ifelse(Source %in% c("USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
                mutate(Group4 = ifelse(Source %in% c("LACDPWG44B(MGD)", "LACDPWG44B(cfs)"), "Yes", NA)) %>%
                mutate(Group5 = ifelse(Source %in% c("LACDPWF313B(MGD)", "LACDPWF313B(cfs)"), "Yes", NA)) %>%
              pivot_longer(Group1:Group5, names_to = "GroupGage", values_to = "Value") %>%
              mutate(GroupCheck = ifelse(Value =="Yes", GroupGage, NA)) %>%
              drop_na(GroupCheck)  %>%
              mutate(GroupQ = gsub("roup", "", GroupGage))

### remove row if plant group does not match Q group
bioMeanQ_longx <- bioMeanQ_long %>%
  mutate(GroupKeep = ifelse(Group == GroupQ, "Yes", "No")) %>%
   select(-GroupCheck, - Value)  %>% filter(GroupKeep == "Yes") 

save(bioMeanQ_longx, file = "output_data/00_bio_Q_matched_groups.RData")
