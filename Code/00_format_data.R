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
  full_join(subs, by = c("PlantID", 'Species', "Group", "LifeForm"))



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

cfs <- unique(hydrox$Source)[c(2,4,6)]
cfs
unique(hydrox$Source)
hydrox <- hydro %>%
  mutate(Date = lubridate::mdy(Date)) %>%
  pivot_longer(cols = "USGS Avg Flow (cfs)":"Discharge SJC002, POM001, & WN001 Combined (MGD)", 
               names_to = "Source", values_to = "Q") %>%
  mutate("Source" = gsub( " ", "", Source)) %>%
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


# Creating Monthly data ----------------------------------------------------

names(hydroxDate)

## split date and make rainfall in with source
hydroxDate <- hydrox %>%
  select(-Group,-Units) %>%
  pivot_wider(names_from = Source, values_from = Q) %>%
  pivot_longer(cols = c("RainfallIntensity":"DischargeSJC002,POM001,&WN001Combined(MGD)"), names_to = "Source", values_to = "Q") %>%
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

unique(test$Year)

monthlyQ <- hydroxDate %>%
  group_by(Year, Month, Source) %>%
  summarise(MinQ = min(na.omit(Q)), MaxQ = max(na.omit(Q)), MedianQ = median(na.omit(Q)), MeanQ = mean(na.omit(Q))) %>%
  mutate(Units = ifelse(Source %in% cfs, "CFS", "MGD")) %>%
  mutate(Units = ifelse(Source == "RainfallIntensity", "INCHES", Units)) %>%
  mutate(Group1 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                        "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                        "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  mutate(Group2 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                        "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                        "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  mutate(Group3 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                        "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                        "SJC-002(MGD)","USGSGauge11087020(MGD)", "USGSAvgFlow(cfs)"), "Yes", NA)) %>%
  mutate(Group4 = ifelse(Source %in% c("RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                       "SJC-002(MGD)","LACDPWG44B(MGD)", "LACDPWG44B(cfs)"), "Yes", NA)) %>%
  mutate(Group5 = ifelse(Source %in% c( "RainfallIntensity","POM-001", "WN-002(Zone1Ditch)", "WN-001",
                                       "DischargeSJC002&POM001Combined(MGD)", "DischargeSJC002,POM001,&WN001Combined(MGD)",
                                       "SJC-002(MGD)", "LACDPWF313B(MGD)", "LACDPWF313B(cfs)"), "Yes", NA)) %>%
  pivot_longer(Group1:Group5, names_to = "Group", values_to = "Value") %>%
  mutate(GroupCheck = ifelse(Value =="Yes", Group, NA)) %>%
  drop_na(GroupCheck) %>% select(-GroupCheck, -Value)   %>% 
  mutate(Group = gsub("roup", "", Group)) %>%
  mutate(Year = gsub("20", "", Year))


head(monthlyQ)


# Join with Bio data ------------------------------------------------------

head(monthlyQ)
head(alldata)

## use inner join for now, can extend to full join after model test, and include lags etc
bioMeanQ <- inner_join(alldata, monthlyQ, by = c("Year", "Group"), relationship = "many-to-many") %>%
  pivot_wider(names_from = Variable, values_from = Values) %>%
  select(-c(MinQ, MaxQ, MedianQ, Units)) %>%
  pivot_wider(names_from = Source, values_from = MeanQ)
names(bioMeanQ)

## rainfall has lots of NAs!!! - happens when no values (ie all NAs) for one month, maybe find an alternative?

save(bioMeanQ, file = "output_data/00_bio_meanQ_data.RData")
