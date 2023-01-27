################################### PREAMBLE #######
library(tidyverse)
library(janitor)
library(hablar)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(ggstatsplot)
library(easystats)
library(MASS)
library(ggeffects)
library(patchwork)
library(magrittr)
library(MuMIn)
library(broom)
library(DHARMa)
library(tidylog)
library(dabestr)
library(pilot)
library(gtools)
library(gt)
library(broom.mixed)
library(tinytex)
library(huxtable)
library(gtsummary)
library(segmented)
library(sjPlot)
library(Hmisc)

##############################
# Fin Growth Growth
##############################

#load data
fin_len <- read_csv("FinGrowth.csv") %>%
  janitor::clean_names()

#growth/day fasting#####

#three-way
m1 <- glmmTMB(fin_len ~ treatment*sex*day_num  + (1|rep/id2), data = findat %>% filter(day_num == "3" | day_num == "7" |day_num == "15" ), REML = T)
summary(m1)
car::Anova(m1, type = "III")
emtrends(m1, specs = pairwise~treatment*sex, var = "day_num", type = "response", adjust = "none")
tab_model(m1,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#male-one
m2 <- glmmTMB(fin_len ~ treatment + day_num  + (1|rep/id2), data = findat %>% filter(sex == "male")  %>% filter(day_num == "3" | day_num == "7" |day_num == "15" ), REML = T)
summary(m2)
car::Anova(m2, type = "III")
emtrends(m2, specs = pairwise~treatment, var = "day_num", type = "response", adjust = "none")

#female-one
m3 <- glmmTMB(fin_len ~ treatment + day_num  + (1|rep/id2), data = findat %>% filter(sex == "female") %>% filter(day_num == "3" | day_num == "7" |day_num == "15" ), REML = T)
summary(m3)
car::Anova(m3, type = "III")
emtrends(m3, specs = pairwise~treatment, var = "day_num", type = "response", adjust = "none")

#growth/day refeeding#####

#two-way
m4 <- glmmTMB(fin_len ~ treatment + sex + day_num + treatment:sex + treatment:day_num + sex:day_num + (1|rep/id2), data = findat  %>% filter(day_num != "3" & day_num != "7 "& day_num != "15"), REML = T)
summary(m4)
car::Anova(m4, type = "III")
emtrends(m4, specs = pairwise~treatment, var = "day_num", type = "response", adjust = "none")
emtrends(m4, specs = pairwise~sex, var = "day_num", type = "response", adjust = "none")
tab_model(m4,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#male-one
m5 <- glmmTMB(fin_len ~ treatment + day_num  + (1|rep/id2), data = findat %>% filter(sex == "male") %>% filter(day_num != "3" & day_num != "7" & day_num != "15"), REML = T)
summary(m5)
car::Anova(m5, type = "III")
emtrends(m5, specs = pairwise~treatment, var = "day_num", type = "response", adjust = "none")

#female-two
m6 <- glmmTMB(fin_len ~ treatment*day_num  +  (1|rep/id2), data = findat %>% filter(sex == "female") %>% filter(day_num != "3" & day_num != "7" & day_num != "15"), REML = T)
summary(m6)
car::Anova(m6, type = "III")
emtrends(m6, specs = pairwise~treatment, var = "day_num", type = "response", adjust = "none")
tab_model(m6,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#GGEffect Graphs########
labs1 <- c("3", "7","15")
labs2 <- c("21", "28", "35")

#for graphing
m1a <- glmmTMB(fin_len ~ treatment*day_num  + (1|rep/id2), data = findat %>% filter(sex == "female") %>% filter(day_num == "3" | day_num == "7" |day_num == "15" ), REML = T)
m1b <- glmmTMB(fin_len ~ treatment*day_num  + (1|rep/id2), data = findat %>% filter(sex == "male")  %>% filter(day_num == "3" | day_num == "7" |day_num == "15" ), REML = T)
m1c <- glmmTMB(fin_len ~ treatment*day_num  + (1|rep/id2), data = findat %>% filter(sex == "male") %>% filter(day_num != "3" & day_num != "7" & day_num != "15"), REML = T)
m1d <- glmmTMB(fin_len ~ treatment*day_num  +  (1|rep/id2), data = findat %>% filter(sex == "female") %>% filter(day_num != "3" & day_num != "7" & day_num != "15"), REML = T)


FemFast <- as.tibble(ggeffects:: ggemmeans(m1a, terms = c("day_num", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m1b, terms = c("day_num", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MalePost <- as.tibble(ggeffects:: ggemmeans(m1c, terms = c("day_num", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))
FemPost <- as.tibble(ggeffects:: ggemmeans(m1d, terms = c("day_num", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)
TotPost <- rbind(MalePost,FemPost)

g1 <- ggplot(TotFast, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  ylab("Fin Length (mm)") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs1, breaks = c(3, 7, 15)) +
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))


g2 <- ggplot(TotPost, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab("Day") +
  ylab("Fin Length (mm)") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs2, breaks = c(21, 28, 35)) +
  theme(strip.text.x = element_blank()) +
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

g1/g2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom') 

ggplot2::ggsave("FinGrowth.png", width = 8, height = 8, device = "png", dpi= 600, path = NULL)


##############################
# Reproduction
##############################

#load data
repodat <- read_csv("ReproductionLong.csv") %>%
  janitor::clean_names()

#ASR fasting#####

#female-one
m7 <- glmmTMB(total_no ~ day_14 + treatment + (1|replicate/id), data = reprodat %>% filter(sex == "Female") %>% filter(day_14 == "7" | day_14 == "15"), family = "nbinom2", zi = ~1, REML = T)
summary(m7)
car::Anova(m7, type = "III")
emmeans(m7, specs = pairwise ~ treatment, type = "response")
tab_model(m7,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#ASR refeeding#####

#female-one
m8 <- glmmTMB(total_no ~ day_14 + treatment + (1|replicate/id), data = reprodat %>% filter(sex == "Female") %>% filter(day_14 != "7" & day_14 != "15"), family = "nbinom2", zi = ~1, REML = T )
summary(m8)
car::Anova(m8, type = "III")
emmeans(m8, specs = pairwise ~ treatment, type = "response")
tab_model(m8,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)


#GGEffect Graphs########

#for graphing
m2a <- glmmTMB(total_no ~ day_14*treatment + (1|replicate/id), data = reprodat %>% filter(sex == "Female") %>% filter(day_14 == "7" | day_14 == "15"), family = "nbinom2", zi = ~1, REML = T)
m2b <- glmmTMB(total_no ~ day_14*treatment + (1|replicate/id), data = reprodat %>% filter(sex == "Male") %>% filter(day_14 == "7" | day_14 == "15"), family = "nbinom2", zi = ~1, REML = T)
m2c <- glmmTMB(total_no ~ day_14*treatment + (1|replicate/id), data = reprodat %>% filter(sex == "Male") %>% filter(day_14 != "7" & day_14 != "15"), family = "nbinom2", zi = ~1, REML = T)
m2d <- glmmTMB(total_no ~ day_14*treatment +  (1|replicate/id), data = reprodat %>% filter(sex == "Female") %>% filter(day_14 != "7" & day_14 != "15"), family = 'nbinom2', zi = ~1, REML = T )


FemFast <- as.tibble(ggeffects:: ggemmeans(m2a, terms = c("day_14", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m2b, terms = c("day_14", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MalePost <- as.tibble(ggeffects:: ggemmeans(m2c, terms = c("day_14", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))
FemPost <- as.tibble(ggeffects:: ggemmeans(m2d, terms = c("day_14", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)
TotPost <- rbind(MalePost,FemPost)

labs1 <- c("7","15")
labs2 <- c("21", "28", "35")


g1 <- ggplot(TotFast, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  ylab("Offspring Count") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs1, breaks = c(7, 15))+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))


g2 <- ggplot(TotPost, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab("Day") +
  ylab("Offspring Count") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs2, breaks = c(21, 28, 35)) +
  theme(strip.text.x = element_blank())+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))


g1/g2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom') 

ggplot2::ggsave("ASR.png", width = 8, height = 8, device = "png", dpi= 600, path = NULL)

##############################
# LRS
##############################

#LRS fasting#####

#male-one
m9 <- glmmTMB(totrep ~ treatment + (1|replicate), zi = ~1, family = "nbinom2", data = reprodatwidefast %>% filter(sex == "Male"), REML = T)
summary(m9)
emmeans(m9, pairwise ~ treatment, type = "response")
tab_model(m9,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#female-one
m10 <- glmmTMB(totrep ~ treatment  + (1|replicate), family = "nbinom2", data = reprodatwidefast %>% filter(sex == "Female"), REML = T)
summary(m10)
emmeans(m10, pairwise ~ treatment, type = "response")
tab_model(m10,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#LRS refeeding#####

#male-one
m11 <- glmmTMB(totrep ~ treatment + (1|replicate), family = "nbinom2", data = reprodatwidepost %>% filter(sex == "Male"), REML = T)
summary(m11)
emmeans(m11, specs = pairwise ~ treatment, type = "response")
tab_model(m11,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#female-one
m12 <- glmmTMB(totrep ~ treatment + (1|replicate),family = "nbinom2", data = reprodatwidepost %>% filter(sex == "Female"), REML = T)
summary(m12)
emmeans(m12, specs = pairwise ~ treatment, type = "response")
tab_model(m12,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)


#GGEffect Graphs########
FemFast <- as.tibble(ggeffects:: ggemmeans(m10, terms = c("treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m9, terms = c("treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MalePost <- as.tibble(ggeffects:: ggemmeans(m11, terms = c("treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))
FemPost <- as.tibble(ggeffects:: ggemmeans(m12, terms = c("treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)
TotPost <- rbind(MalePost,FemPost)

g3 <- ggplot(TotFast, aes(x = x, y = predicted, colour = x)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = x), show.legend = F) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, show.legend = F) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  ylab("Total Offspring Number") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment")+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted")) +
  scale_x_discrete(labels=c("Fed","Fasted"))


g4 <-  ggplot(TotPost, aes(x = x, y = predicted, colour = x)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = x), show.legend = F) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, show.legend = F) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab("Treatment") +
  ylab("Total Offspring Number") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  theme(strip.text.x = element_blank()) +
scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted")) +
  scale_x_discrete(labels=c("Fed","Fasted"))


(g1|g3)/(g2|g4) + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom') 

ggplot2::ggsave("LRSRepr.png", width = 8, height = 8, device = "png", dpi= 600, path = NULL)

##############################
# Survival
##############################

#load data
ehDat1 <- read_csv("Survival.csv") %>%
  janitor::clean_names()

#Survival fasting######

#male-two
m13 <- glmmTMB(lifespan ~ treatment*day +  (1|replicate/id), family = "binomial", data = ehDat1 %>% filter(sex == "Male") %>% filter(day == "7" | day == "15"), REML = T)
summary(m13)
car::Anova(m13, type = "III")
emtrends(m13, specs = pairwise~treatment, var = "day", type = "response", adjust = "none") 
tab_model(m13 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#female-two
m14 <- glmmTMB(lifespan ~ treatment*day +  (1|replicate/id), family = "binomial", data = ehDat1 %>% filter(sex == "Female") %>% filter(day == "7" | day == "15"), REML = T)
summary(m14)
car::Anova(m14, type = "III")
emtrends(m14, specs = pairwise~treatment, var = "day", type = "response", adjust = "none")
tab_model(m14 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Survival refeeding######

#male-two
m15 <- glmmTMB(lifespan ~ treatment*day +  (1|replicate/id), family = "binomial", data = ehDat1 %>% filter(sex == "Male") %>% filter(day != "7" & day != "15"), REML = T)
summary(m15)
car::Anova(m15, type = "III")
emtrends(m15, specs = pairwise~treatment, var = "day", type = "response", adjust = "none") 
tab_model(m15 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#female-two
m16<- glmmTMB(lifespan ~ treatment*day +  (1|replicate/id), family = "binomial", data = ehDat1 %>% filter(sex == "Female") %>% filter(day != "7" & day != "15"), REML = T)
summary(m16)
car::Anova(m16, type = "III")
emtrends(m16, specs = pairwise~treatment, var = "day", type = "response", adjust = "none")
tab_model(m16 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)


#GGEffect Graphs########
FemFast <- as.tibble(ggeffects:: ggemmeans(m14, terms = c("day", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m13, terms = c("day", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MalePost <- as.tibble(ggeffects:: ggemmeans(m15, terms = c("day", "treatment"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))
FemPost <- as.tibble(ggeffects:: ggemmeans(m16, terms = c("day", "treatment"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Post", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)
TotPost <- rbind(MalePost,FemPost)


labs1 <- c("7","15")
labs2 <- c("21", "28", "35")


g1 <- ggplot(TotFast, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  ylab("Probability of an Egg Surviving") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs1, breaks = c(7, 15))+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

g2 <- ggplot(TotPost, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab("Day") +
  ylab("Probability of an Egg Surviving") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(labels = labs2, breaks = c(21, 28, 35)) +
  theme(strip.text.x = element_blank())+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

g1/g2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom') 

ggplot2::ggsave("Egg Prob.png", width = 8, height = 8, device = "png", dpi= 600, path = NULL)

##############################
# Sperm
##############################

#load data
spermdat <- read_csv("Sperm.csv") %>%
  janitor::clean_names()

#Sperm fasting######

#VLC-two
m17 <- glmmTMB(log(VCL) ~ Treat*Day + block + (1|fish_id), family = gaussian, data = spermdat %>% filter(Day == "7"|Day == "15"), REML = T)
summary(m17) 
emtrends(m17, specs = pairwise~Treat, var = "Day", type = "response", adjust = "none")
tab_model(m17 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#VAP-two
m18 <- glmmTMB(log(VAP) ~ Treat*Day + block + (1|fish_id), family = gaussian, data = spermdat %>% filter(Day == "7"|Day == "15"), REML = T)
summary(m18) 
emtrends(m18, specs = pairwise~Treat, var = "Day", type = "response", adjust = "none")
tab_model(m18 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#Sperm refeeding#####

#VLC-two
m19 <- glmmTMB(log(VCL) ~ Treat*Day + block + (1|fish_id), family = gaussian, data = spermdat %>% filter(Day != "7" & Day != "15"), REML = T)
summary(m19)  
emtrends(m19, specs = pairwise~Treat, var = "Day", type = "response", adjust = "none")
tab_model(m19 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#VAP-two
m20 <- glmmTMB(log(VAP) ~ Treat*Day + block + (1|fish_id), family = gaussian, data = spermdat %>% filter(Day != "7" & Day != "15"), REML = T)
summary(m20)
emtrends(m20, specs = pairwise~Treat, var = "Day", type = "response", adjust = "none")
tab_model(m20 ,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#GGEffect Graphs########
MaleFastVCL <- as.tibble(ggeffects:: ggemmeans(m17, terms = c("Day", "Treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fasting", nrow(.)))
MalePostVCL <- as.tibble(ggeffects:: ggemmeans(18, terms = c("Day", "Treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Refeeding", nrow(.)))

MaleFastVAP <- as.tibble(ggeffects:: ggemmeans(m19, terms = c("Day", "Treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fasting", nrow(.)))
MalePostVAP <- as.tibble(ggeffects:: ggemmeans(m20, terms = c("Day", "Treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Refeeding", nrow(.)))

TotFastVCL <- rbind(MaleFastVCL,MalePostVCL)
TotFastVAP <- rbind(MaleFastVAP,MalePostVAP)

labs1 <- c(7,15)
labs2 <- c(21, 28, 35)

my_breaks <- function(x) { if (max(x) < 20) labs2 else labs1 }

g1 <- ggplot(TotFastVCL, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_wrap(.~period, scales = "free_x") +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Velocity Curvilinear ( ", mu, "m/s )"))) +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(breaks = c(7,15, 21,28,35)) +
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

g2 <- ggplot(TotFastVAP, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_wrap(.~period, scales = "free_x") +
  theme_bw() +
  xlab(NULL) +
  ylab(expression(paste("Average Path Velocity ( ", mu, "m/s )"))) +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment") +
  scale_x_continuous(breaks = c(7,15, 21,28,35)) +
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

g1/g2 + plot_annotation(tag_levels = 'A') + plot_layout(guides = "collect") & theme(legend.position = 'bottom') 

ggplot2::ggsave("Sperm.png", width = 8, height = 8, device = "png", dpi= 600, path = NULL)


##############################
# Growth
##############################

#load data
dat <- read_csv("FryGrowth.csv") %>%
  janitor::clean_names()

#FryGrowth Per Day#######

#two-way
m21 <- glmmTMB(len_mm ~ treat + dpf + sex + treat:dpf + treat:sex + sex:dpf + (1|rep/id), family = "gaussian", data = dat, REML = T)
summary(m21)
car::Anova(m21, type = "III")
emtrends(m21, specs = pairwise~sex, var = "dpf", type = "response")
emtrends(m21, specs =pairwise ~treat, var = "dpf", type = "response")
tab_model(m21,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#male-two
m22 <- glmmTMB(len_mm ~ treat + dpf + treat:dpf + (1|rep/id), family = "gaussian", data = dat %>% filter(sex == "male"), REML = T)
summary(m22)
car::Anova(m22, type = "III")
emtrends(m22, specs = pairwise~treat, var = "dpf", type = "response")
tab_model(m22,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)

#female-two
m23 <- glmmTMB(len_mm ~ treat + dpf + treat:dpf + (1|rep/sample), family = "gaussian", data = dat %>% filter(sex == "female"), REML = T)
summary(m23)
car::Anova(m23, type = "III")
emtrends(m23, specs = pairwise~treat, var = "dpf", type = "response")
tab_model(m23,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F)


#GGEffect Graphs########
FemFast <- as.tibble(ggeffects:: ggemmeans(m23, terms = c("dpf", "treat"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m22, terms = c("dpf", "treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)

g1 <- ggplot(TotFast, aes(x = x, y = predicted, colour = group, linetype = group)) +
  geom_line() +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.15) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  ylab("Fry Length (mm)") +
  labs(colour = "Treatment", linetype = "Treatment", fill = "Treatment")+
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted"))

#FryGrowth Per Hour#######

#One-way
m24 <- glmmTMB(growth_day ~ treat + sex2 + rep, family = "gaussian", data = datwide, REML = T)
summary(m24)
car::Anova(m24, type = "III")
emmeans(m24, specs = pairwise~sex2, type = "response")
emmeans(m24, specs = pairwise~treat, type = "response")
tab_model(m24,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F, digits = 5)

#male-one
m25 <- glmmTMB(growth_day ~ treat + rep, family = "gaussian", data = datwide %>% filter(sex2 == "Male"), REML = T)
summary(m25)
emmeans(m25, specs = pairwise~treat, type = "response")
tab_model(m25,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F, digits = 5)

#female-one
m26 <- glmmTMB(growth_day ~ treat + rep, family = "gaussian", data = datwide %>% filter(sex2 == "Female"), REML = T)
summary(m26)
emmeans(m26, specs = pairwise~treat, type = "response")
tab_model(m26,transform = NULL, show.ci = F, show.se = T, show.r2 = F, show.stat = T, show.icc = F, digits = 5)

#GGEffect Graphs########
FemFast <- as.tibble(ggeffects:: ggemmeans(m26, terms = c("treat"))) %>%
  mutate(sex = rep("Female", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))
MaleFast <- as.tibble(ggeffects:: ggemmeans(m125, terms = c("treat"))) %>%
  mutate(sex = rep("Male", nrow(.))) %>%
  mutate(period = rep("Fast", nrow(.)))

TotFast <- rbind(FemFast,MaleFast)

g3 <- ggplot(TotFast, aes(x = x, y = predicted, colour = x)) +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high, colour = x), show.legend = F) +
  stat_summary(fun.data = "mean_cl_boot", geom = "point", size = 2, show.legend = F) +
  facet_grid(.~sex) +
  theme_bw() +
  xlab(NULL) +
  xlab("Treatment") +
  ylab("Fry Growth per Hour") +
  labs(colour = "Treatment") +
  scale_colour_discrete(labels=c("Fed", "Fasted")) +
  scale_linetype_discrete(labels=c("Fed", "Fasted")) +
  scale_fill_discrete(labels=c("Fed", "Fasted")) +
  scale_x_discrete(labels=c("Fed","Fasted"))


g1/g3 + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = 'bottom') 

ggplot2::ggsave("FRGrowth.png", width = 8, height = 10, device = "png", dpi= 600, path = NULL)
