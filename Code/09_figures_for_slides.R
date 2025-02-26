## figures for slides

library(tidyverse)
library(tidylog)


## col blind pallette
cols <- palette.colors()
cols
## data
binData <- read.csv( "05_prob_stress_data_for_figure.csv")
head(binData)

unique(binData$Type)
## plot just mod curve. no Q labels

s1 <- ggplot(subset(binData, Type =="Overall"), aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, , colour = NA) +
  scale_fill_manual(values = cols[c(8)]) +
  scale_colour_manual(values = cols[c( 8)]) +
  # annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
  #          alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Q", y=0.8), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=2.5, label="Current Discharge Range", y=0.80), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=5, label="Fall Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=-12, label="Spring Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  # geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  # geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s1

file.name1 <- "Figures/09_probability_of_stress_1_curve.jpg"
ggsave(s1, filename=file.name1, dpi=300, height=10, width=15)

## plot just mod and wet year curve. no Q labels

s2 <- ggplot(subset(binData, Type %in% c("Overall", "Dry Year: 2022")), aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
  scale_fill_manual(values = cols[c(5, 8)]) +
  scale_colour_manual(values = cols[c(5, 8)]) +
  # annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
  #          alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Q", y=0.80), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=2.5, label="Current Discharge Range", y=0.80), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=5, label="Fall Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=-12, label="Spring Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  # geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  # geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s2

file.name1 <- "Figures/09_probability_of_stress_2_curve.jpg"
ggsave(s2, filename=file.name1, dpi=300, height=10, width=15)

## plot all curves. no Q labels

s3 <- ggplot(binData, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
  scale_fill_manual(values = cols[c(3, 5, 8)]) +
  scale_colour_manual(values = cols[c(3, 5, 8)]) +
  # annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
  #          alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Q", y=0.80), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=2.5, label="Current Discharge Range", y=0.80), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=5, label="Fall Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=-12, label="Spring Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  # geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  # geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s3

file.name1 <- "Figures/09_probability_of_stress_3_curve.jpg"
ggsave(s3, filename=file.name1, dpi=300, height=10, width=15)

## plot all curves 

s5 <- ggplot(binData, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
  scale_fill_manual(values = cols[c(3, 5, 8)]) +
  scale_colour_manual(values = cols[c(3, 5, 8)]) +
  # annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
  #          alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Q", y=0.80), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=2.5, label="Current Discharge Range", y=0.80), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  # geom_text(aes(x=5, label="Fall Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=-12, label="Spring Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  # geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s5

file.name1 <- "Figures/09_probability_of_stress_3_curve_thresholds.jpg"
ggsave(s5, filename=file.name1, dpi=300, height=10, width=15)


## plot all data on one plot

s4 <- ggplot(binData, aes(x = CombQ, y = ProbabilityOfStress)) +
  geom_smooth(aes(col = Type, fill = Type))  +
  geom_ribbon( aes(col = Type, fill = Type,ymin = lwr, ymax = upr), alpha = .15, colour = NA) +
  scale_fill_manual(values = cols[c(3,5, 8)]) +
  scale_colour_manual(values = cols[c(3,5, 8)]) +
  annotate("rect", xmin = 2, xmax = 2.5, ymin = -Inf, ymax = Inf,
           alpha = .2, fill=cols[7]) +
  geom_text(aes(x=0, label="Baseline Discharge", y=0.80), colour=cols[9], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=2.5, label="Current Discharge Range", y=0.80), colour=cols[7], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=5, label="Fall Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_text(aes(x=-12, label="Spring Current Discharge", y=0.80), colour=cols[4], angle=90, vjust = 1.3, size=6)+
  geom_vline(xintercept = 0, col = cols[9], lty = "dashed") +
  geom_vline(xintercept = 5, col = cols[4], lty = "dashed") + ## fall line
  geom_vline(xintercept = -12, col = cols[4], lty = "dashed") + ## spring line
  # facet_wrap(~Year, scale = "free_x") +
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 25)) +
  scale_x_continuous(name = "Delta Combined Q (MGD)") +
  scale_y_continuous(name = "Probability of Stress (SWP)")

s4

file.name1 <- "Figures/09_probability_of_stress_overall_current_delta_range_ALL.jpg"
ggsave(s4, filename=file.name1, dpi=300, height=10, width=15)
