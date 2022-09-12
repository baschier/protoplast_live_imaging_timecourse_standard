## Analysis of image data collected 6/29/22 of B73 GPS2 and GPS1-NR protoplasts treated with 0, 1, 10, and 100uM GA3
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(agricolae)
library(emmeans)
library(ggbeeswarm)
library(readr)
library(stringr)


A1_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A1_before_Results.csv", col_select = c(1:5))
A1_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A1_after_Results.csv", col_select = c(1:5))
A2_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A2_before_Results.csv", col_select = c(1:5))
A2_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A2_after_Results.csv", col_select = c(1:5))
A3_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A3_before_Results.csv", col_select = c(1:5))
A3_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A3_after_Results.csv", col_select = c(1:5))
A4_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A4_before_Results.csv", col_select = c(1:5))
A4_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/A4_after_Results.csv", col_select = c(1:5))
B1_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B1_before_Results.csv", col_select = c(1:5))
B1_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B1_after_Results.csv", col_select = c(1:5))
B2_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B2_before_Results.csv", col_select = c(1:5))
B2_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B2_after_Results.csv", col_select = c(1:5))
B3_before_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B3_before_Results.csv", col_select = c(1:5))
B3_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B3_after_Results.csv", col_select = c(1:5))
B4_before_Results <- read_csv("/Volumes/Seagate/2022_6_29_B73_GPS2_GPS1NR/data/B4_before_Results.csv", col_select = c(1:5))
B4_after_Results <- read_csv("~/Downloads/2022_6_29_B73_GPS2_GPS1NR_data/data/B4_after_Results.csv", col_select = c(1:5))

B73_6_30_22 <- rbind(A1_before_Results, A1_after_Results, A2_before_Results, A2_after_Results, A3_before_Results, A3_after_Results, A4_before_Results, A4_after_Results, B1_before_Results, B1_after_Results, B2_before_Results, B1_after_Results, B2_before_Results, B2_after_Results, B3_before_Results, B3_after_Results, B4_before_Results, B4_after_Results)
B73_6_30_22 <- subset(B73_6_30_22, Area <= 300 & Area > 0)

B73_6_30_22$TIME <- "value"
B73_6_30_22 <- B73_6_30_22 %>% mutate(TIME = case_when(
  grepl(pattern = "Result of BEFORE", x = Label) ~ "-10min",
  grepl(pattern = "slice0001_Result of After", x = Label) ~ "0h",
  grepl(pattern = "slice0002_Result of After", x = Label) ~ "2h",
  grepl(pattern = "slice0003_Result of After", x = Label) ~ "4h",
  grepl(pattern = "slice0004_Result of After", x = Label) ~ "6h",
  grepl(pattern = "slice0005_Result of After", x = Label) ~ "8h",
  grepl(pattern = "slice0006_Result of After", x = Label) ~ "10h",
  grepl(pattern = "slice0007_Result of After", x = Label) ~ "12h",
  grepl(pattern = "slice0008_Result of After", x = Label) ~ "14h",
  grepl(pattern = "slice0009_Result of After", x = Label) ~ "16h",
  grepl(pattern = "slice0010_Result of After", x = Label) ~ "18h",
  grepl(pattern = "slice0011_Result of After", x = Label) ~ "20h",
  grepl(pattern = "slice0012_Result of After", x = Label) ~ "22h",
  grepl(pattern = "slice0013_Result of After", x = Label) ~ "24h",
))

B73_6_30_22$TREATMENT <- "value"  
B73_6_30_22 <- B73_6_30_22 %>% mutate(TREATMENT = case_when(
  grepl(pattern = "GA3_A1", x = Label) ~ "0uM GA3",
  grepl(pattern = "GA3_A2", x = Label) ~ "1uM GA3",
  grepl(pattern = "GA3_A3", x = Label) ~ "10uM GA3",
  grepl(pattern = "GA3_A4", x = Label) ~ "100uM GA3",
  grepl(pattern = "GA3_B1", x = Label) ~ "0uM GA3",
  grepl(pattern = "GA3_B2", x = Label) ~ "1uM GA3",
  grepl(pattern = "GA3_B3", x = Label) ~ "10uM GA3",
  grepl(pattern = "GA3_B4", x = Label) ~ "100uM GA3",
))

B73_6_30_22$Construct <- "value"  
B73_6_30_22 <- B73_6_30_22 %>% mutate(Construct = case_when(
  grepl(pattern = "GA3_A", x = Label) ~ "GPS2",
  grepl(pattern = "GA3_B", x = Label) ~ "GPS1-NR",
))

colnames(B73_6_30_22)[4] <- "Ratio"
B73_6_30_22$Ratio <- as.numeric(B73_6_30_22$Ratio)

write_csv(B73_6_30_22, "2022_6_30_all_B73_GPS2_GPS1NR.csv")

B73_6_30_22 <- within(B73_6_30_22, {
  TIME <- factor(TIME, levels = c("-10min", "0h", "2h", "4h", "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h", "22h", "24h"))})

B73_6_30_22 <- within(B73_6_30_22, {
  Construct <- factor(B73_6_30_22$Construct)})

B73_6_30_22 <- within(B73_6_30_22, {
  TREATMENT <- factor(B73_6_30_22$TREATMENT)})

data_summary <- B73_6_30_22 %>% group_by(TIME, TREATMENT, Construct) %>% summarise(n = n(), mean = mean(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem) 

data_summary <- within(data_summary, {
  TIME <- factor(TIME, levels = c("-10min", "0h", "2h", "4h", "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h", "22h", "24h"))
  TREATMENT <- factor(TREATMENT, levels = c("0uM GA3", "1uM GA3", "10uM GA3", "100uM GA3"))})

plottitle <- "Maize B73 GPS2 and GPS1-NR FRET response to GA3"
ggplot(data_summary, aes(x = TIME, y = mean, color = Construct, group = TREATMENT, shape = TREATMENT)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, size = 1, position = position_dodge(width = 0.3)) + 
  geom_point(lwd = 4, position = position_dodge(width = 0.3)) + labs(x = "TIME(minutes)", y = "Ratio (FRET/CFP)", title = str_wrap(plottitle, 30)) +  theme(axis.text = element_text(size = 22)) + theme(axis.title = element_text(size = 28)) + theme(plot.title = element_text(size = 36)) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 22)) + scale_color_manual(values = c("GPS1-NR" = "Orange", "GPS2" = "black")) + theme(plot.title = element_text(hjust = 0.25)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2.75))

#wondering why GPS1-NR 100uM GA3 at -10min FRET is so low, so I want to look at that data
GPS1NR_100uM <- subset(GPS1NR, TREATMENT == "100uM GA3")
  
# Do statistical analysis for GPS2 construct first
GPS2 <- subset(B73_6_30_22, Construct == "GPS2")
B73_GPS2_2way_aov <- aov(Ratio ~ TREATMENT + TIME, data = GPS2)
B73_GPS2_2way_aov_treatment_x_time <- aov(Ratio ~ TREATMENT*TIME, data = GPS2)
summary(B73_GPS2_2way_aov)
summary(B73_GPS2_2way_aov_treatment_x_time)
HSD.test(B73_GPS2_2way_aov, "TREATMENT", group = T, console = T)

GPS1NR <- subset(B73_6_30_22, Construct == "GPS1-NR")
B73_GPS1NR_2way_aov <- aov(Ratio ~ TREATMENT + TIME, data = GPS1NR)
B73_GPS1NR_2way_aov_treatment_x_time <- aov(Ratio ~ TREATMENT*TIME, data = GPS1NR)
summary(B73_GPS1NR_2way_aov)
summary(B73_GPS1NR_2way_aov_treatment_x_time)
HSD.test(B73_GPS1NR_2way_aov, "TREATMENT", group = T, console = T)

#Both construct anova and post hoc, post hoc doesn't work for interaction 
GPS.aov <- aov(Ratio ~ Construct*TREATMENT + TIME, data = B73_6_30_22)
summary(GPS.aov)
HSD.test(GPS.aov, "Construct:TREATMENT", group = T, console = T)

#now I want to see the survival of protoplasts from 0h-24h for maize
ggplot(data_summary, aes(x = TIME, y = n, color = Construct, group = TREATMENT, shape = TREATMENT)) + geom_point() + labs(x = "Time", y = "N protoplasts", title = "Maize B73 protoplast nuclei numbers") + scale_color_manual(values = c("GPS1-NR" = "darkorange2", "GPS2" = "black")) +
  theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 26)) + theme(plot.title = element_text(size = 32)) + theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.25)) + theme(panel.background = element_rect(fill = "papayawhip"))




  