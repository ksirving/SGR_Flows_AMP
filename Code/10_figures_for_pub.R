## Figures for publication.

library(tidyverse)
library(tidylog)

#Figure 5 to 14

## increase dpi
## font size as standard


# SWP ---------------------------------------------------------------------

# Figure 5 - variable importance

## data 
ImpMean <- read.csv("output_data/05_SWP_var_imp.csv")
## order in increasing values

ImpMean$VariableHuman <- factor(ImpMean$VariableHuman, levels=ImpMean[order(ImpMean$MeanImpPerc,decreasing=F),]$VariableHuman)
## add mod peformance
ImpMean <- ImpMean %>%
  mutate(RFVarExpl = mean(rf$rsq))
ImpMean

## plot

i1 <- ggplot(ImpMean, aes(x = MeanImpPerc, y = VariableHuman)) +
  geom_point() +
  scale_x_continuous("Relative Importance (%)") +
  scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
    # plot.title = element_text(size = 13),
    # legend.text = element_text(size = 9),
    # legend.title = element_text(size = 10)
  )
i1

file.name1 <- "FinalFigures/Figure5_SWP_varImp.jpg"
ggsave(i1, filename=file.name1, dpi=600, height=5, width=8)

## Figure 6 - mixed model 
rf.data.rescalex <- read.csv("FigureData/05_swp_mixedModel_for_figures.csv")

m1 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = SJC002_POM001Combined, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Combined Q (MGD)") +
  scale_y_continuous(name = "SWP: Mixed Model Fit (atm)") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  )
# geom_line(aes(col = Species), linewidth = 1)  

m1

file.name1 <- "FinalFigures/Figure6_mixedSWP.jpg"
ggsave(m1, file = file.name1, dpi=600, height=10, width=12)

# Figure 9 - probably of stress all 

## data
binData <- read.csv( "FigureData/09_SWP_prob_stress_figure.csv")
head(binData)

## col blind pallette
cols <- palette.colors()
cols

# binData <- binData %>%
#   # filter(Type == "Overall") %>%
#   select(Type, CombQ, ProbabilityOfStress, upr,lwr) %>%
#   distinct()

s4 <- ggplot(binData, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
  scale_fill_manual(values = cols[c(3,5, 8)]) +
  scale_colour_manual(values = cols[c(3,5, 8)]) +
  annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
           alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Discharge", y=0.75), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=2.5, label="Current Discharge Range", y=0.72), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=5, label="Fall Current Discharge", y=0.72), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=-12, label="Spring Current Discharge", y=0.72), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        # plot.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s4

file.name1 <- "FinalFigures/Figure9_swp_prob_stress.jpg"
ggsave(s4, filename=file.name1, dpi=600, height=10, width=15)

# Figure 10 = Probability of stress species

write.csv(spDataPx, "FigureData/05_prob_stress_per_species.csv")

## plot

s1 <- ggplot(spDataPx, aes(x = CombQ, y = ProbabilityOfStress, group = Species, colour = Species)) +
  geom_line()  +
  # facet_wrap(~Season, scale = "free_x") +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress") +
  theme_bw()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    # plot.title = element_text(size = 13),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

s1

file.name1 <- "FinalFigures/Figure10_prob_stress_species.jpg"
ggsave(s1, filename=file.name1, dpi=600, height=5, width=8)


#  CV ---------------------------------------------------------------------

# Figure 7: Variable Importance

ImpMean <- read.csv("output_data/06_CV_var_imp.csv") 
ImpMean

ImpMean$VariableHuman <- factor(ImpMean$VariableHuman, levels=ImpMean[order(ImpMean$MeanImpPerc,decreasing=F),]$VariableHuman)
## add mod peformance
ImpMean <- ImpMean %>%
  mutate(RFVarExpl = mean(rf$rsq))
ImpMean

## plot

i1 <- ggplot(ImpMean, aes(x = MeanImpPerc, y = VariableHuman)) +
  geom_point() +
  scale_x_continuous("Relative Importance (%)") +
  scale_y_discrete("") +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
    # plot.title = element_text(size = 13),
    # legend.text = element_text(size = 9),
    # legend.title = element_text(size = 10)
  )
i1

file.name1 <- "FinalFigures/Figure7_CV_varImp.jpg"
ggsave(i1, filename=file.name1, dpi=600, height=5, width=8)

# Figure 8: Mixed model

rf.data.rescalex <- read.csv("FigureData/06_mixed_mod_perSpecies.csv")

m1 <- ggplot(data = rf.data.rescalex, aes(y = fit.c, x = SJC002_POM001Combined, group = Species, col = Species)) +
  geom_smooth(method = "lm") +
  facet_wrap(~Season) +
  scale_x_continuous(name = "Combined Q (MGD)") +
  scale_y_continuous(name = "CV: Mixed Model Fit (%)") +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    plot.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 15)
  )
# geom_line(aes(col = Species), linewidth = 1)  


m1

file.name1 <- "FinalFigures/Figure8_CV_MixedMod.jpg"
ggsave(m1, file = file.name1, dpi=600, height=10, width=12)

# Figure 11: Prob stress overall 


s3 <- ggplot(binDataP, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth()  +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
           alpha = .2, fill='red') +
  theme_classic()+
  geom_text(aes(x=0, label="Baseline Discharge", y=0.70), colour="gray30", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=2.5, label="Current Discharge Range", y=0.70), colour="red", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=5, label="Fall Current Discharge", y=0.70), colour="forestgreen", angle=90, vjust = 1.3, size=4)+
  geom_text(aes(x=-12, label="Spring Current Discharge", y=0.70), colour="forestgreen", angle=90, vjust = 1.3, size=4)+
  geom_vline(xintercept = 0, col = "gray30", lty = "dashed") +
  geom_vline(xintercept = 5, col = "forestgreen", lty = "dashed") + ## fall line
  geom_vline(xintercept = -12, col = "forestgreen", lty = "dashed") + ## spring line
  # facet_wrap(~Season, scale = "free_x") +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (CV)", limits = c(0.3,0.8)) +
  theme_classic()+
  theme(legend.position = "none",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          # plot.title = element_text(size = 13),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 12)) 
  

s3

file.name1 <- "FinalFigures/Figure11_probability_of_stress_overall_current_delta_range.jpg"
ggsave(s3, filename=file.name1, dpi=600, height=5, width=8)

# Figure 12: Prob stress per species

spDataPx <- read.csv("FigureData/06_CV_prob_stress_species.csv")
head(spDataPx)

## plot

s1 <- ggplot(spDataPx, aes(x = CombQ, y = ProbabilityOfStress, group = Species, colour = Species)) +
  geom_line()  +
  # facet_wrap(~Season, scale = "free_x") +
  geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress") +
  theme_classic()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    # plot.title = element_text(size = 13),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

s1

file.name1 <- "FinalFigures/Figure12_CV_probability_of_stress_per_species.jpg"
ggsave(s1, filename=file.name1, dpi=600, height=5, width=8)


# Power analysis ----------------------------------------------------------

# Figure 13

alldf <- read.csv( "FigureData/03b_powerAnalysis_Fig13.csv")

# alldf <- alldf %>%
#   filter(alpha == 0.05)

p1 <- ggplot(alldf, aes(x=Year, y = Power, group = Stressor, col = Stressor)) +
  geom_point(aes(group = Stressor, col = Stressor)) +
  facet_wrap(~Stressor)+
  theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_hline(yintercept=0.8, col = "red", linetype="dashed") +
  # geom_vline(xintercept="2028", col = "red", linetype="dashed") +
  theme(legend.position = "none") +
  theme_classic()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
    # plot.title = element_text(size = 13),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 12)
  )

p1

file.name1 <- "FinalFigures/Figure13_CombQ_sensitivity.jpg"
ggsave(p1, filename=file.name1, dpi=600, height=5, width=8)

# Figure 14

alldf2 <- read.csv("FigureData/03c_power_analysis_sens.csv")

p3 <- ggplot(alldf2, aes(y=as.factor(round(firstYear, digits=0)), x = SignificanceLevel, group = Stressor, col = Stressor)) +
  geom_point(aes(group = Stressor, col = Stressor)) +
  facet_wrap(~Type)+
  ylab("Year") + xlab("Significance Level") +
  theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_blank())  +
  # theme_bw()+
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15)
    # plot.title = element_text(size = 13),
    # legend.text = element_text(size = 10),
    # legend.title = element_text(size = 12)
  )
# geom_hline(yintercept=0.8, col = "red", linetype="dashed") +
# geom_vline(xintercept="2028", col = "red", linetype="dashed") +
# theme(legend.position = "none")


p3

file.name1 <- "FinalFigures/Figure14_poerAn_sens.jpg"
ggsave(p3, filename=file.name1, dpi=600, height=5, width=8)

