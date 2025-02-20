"
THAMP_R_EEG.R
Arun Asthagiri
01-21-2025
"

library(R.matlab)

mat_filepath = "./../mat_files/"
filename = "PLV_to_R.mat"
filename_plv = "PLV_to_R.mat"
data <- readMat(paste(mat_filepath, filename, sep=""))
data_PLV <- readMat(paste(mat_filepath, filename_plv, sep=""))

library(tidyr)
library(reshape2)
library(dplyr)

mod_RT <- melt(data$mod.SART.RT) %>% 
  setNames(c("sub.idx","song.num", "RT"))

mod_RTCV <- melt(data$mod.SART.RTCV) %>% 
  setNames(c("sub.idx","song.num", "RTCV"))

mod_PLV <- melt(data_PLV$sart.mod.PLV.modfreq) %>% 
  setNames(c("sub.idx","song.num", "electrode", "PLV"))

mod_data <- mod_RT %>% 
  merge(mod_RTCV, by=c("sub.idx", "song.num")) %>% 
  merge(mod_PLV, by=c("sub.idx", "song.num")) %>% 
  mutate(condition="Mod")

unmod_RT <- melt(data$unmod.SART.RT) %>% 
  setNames(c("sub.idx","song.num", "RT")) 

unmod_RTCV <- melt(data$unmod.SART.RTCV) %>% 
  setNames(c("sub.idx","song.num", "RTCV"))

unmod_PLV <- melt(data_PLV$sart.unmod.PLV.modfreq) %>% 
  setNames(c("sub.idx","song.num", "electrode", "PLV"))

unmod_data <- unmod_RT %>% 
  merge(unmod_RTCV, by=c("sub.idx", "song.num")) %>% 
  merge(unmod_PLV, by=c("sub.idx", "song.num")) %>% 
  mutate(condition="Unmod")

all_data<-rbind(mod_data, unmod_data)

mean_electrode_data<-all_data %>% 
  # mutate(RTCV = log1p(RTCV)) %>% 
  group_by(sub.idx, song.num, condition) %>% 
  mutate(m_PLV = (mean(PLV))) %>% 
  select(-c(electrode, PLV)) %>% 
  distinct() %>% ungroup() %>%
  group_by(sub.idx) 
# mutate(m_PLV = m_PLV-mean(m_PLV)) %>%
# mutate(RTCV = RTCV-mean(RTCV))
metadata_filepath <- "./../metadata/"
qualtrics_filename<- "qualtrics.csv"
data_qualtrics <- read.csv(paste(metadata_filepath, qualtrics_filename,sep=""), header=FALSE)
colnames(data_qualtrics) <- c("participantID", "ASRS", "eBMRQ")
data_qualtrics<-data_qualtrics %>% mutate(sub.idx = row_number())
mean_electrode_data<-mean_electrode_data%>%merge(data_qualtrics, by=c("sub.idx"))%>%
  drop_na() 
mean_electrode_data$ASRS <- as.factor(mean_electrode_data$ASRS)
mean_electrode_data$condition <- relevel(as.factor(mean_electrode_data$condition), ref="Unmod")
mean_electrode_data<-mean_electrode_data %>% mutate(eBMRQ.tertile = factor(ntile(eBMRQ,3)))
levels(mean_electrode_data$eBMRQ.tertile) <- c("Low", "Med", "High")
levels(mean_electrode_data$ASRS) <- c("ASRS.Negative","ASRS.Positive")
# rescale variables
# mean_electrode_data$RTCV <- scale(mean_electrode_data$RTCV, center = FALSE, scale = max(mean_electrode_data$RTCV))
# mean_electrode_data$m_PLV <- scale(mean_electrode_data$m_PLV, center = FALSE, scale = max(mean_electrode_data$m_PLV))
mean_electrode_data$eBMRQ <- scale(mean_electrode_data$eBMRQ, center = TRUE, scale = max(mean_electrode_data$eBMRQ))




library(ggplot2)
library(ggridges)
ggplot(mean_electrode_data, aes(x=m_PLV, fill=condition)) +
  geom_density(alpha=.5)

ggplot(mean_electrode_data, aes(x=condition, y=m_PLV, color=condition)) +
  geom_violin()+
  # geom_point() +
  geom_jitter() +
  geom_line(aes(group = sub.idx), color = "black", alpha = 0.1) 

library(see)
ggplot(mean_electrode_data, aes(x=condition, y=m_PLV, color=condition)) +
  geom_violinhalf(aes(group = interaction(condition, ASRS))) +
  theme_classic()

ggplot(mean_electrode_data, aes(x=condition, y=m_PLV, fill=ASRS)) +
  geom_violinhalf(aes(group = interaction(condition, ASRS)), 
                  position = position_dodge(width = 0.3), 
                  alpha = 0.5) +
  theme_classic() +
  labs(x = "Condition", y = "m_PLV", fill = "ASRS Group")
  geom_jitter() 
  
library(gghalves)
ggplot(data=mean_electrode_data, 
       aes(x=condition, y=m_PLV, split = ASRS, fill=ASRS)) +
    geom_half_violin(
      position = "identity", 
      alpha = 1) +
    theme_classic() +
    labs(x = "Condition", y = "m_PLV", fill = "ASRS Group")
  
# IMPORTANT PLOT
ggplot(data=mean_electrode_data%>%filter(condition=="Mod"), 
  aes(x=eBMRQ.tertile, y=m_PLV, split = ASRS, fill=ASRS)) +
  geom_half_violin(
    position=position_nudge(x=0),
    alpha = 1) +
  geom_half_boxplot(
    position = position_dodge2(preserve="single"),
    fill="white", side = "r", 
    width = 0.4, 
    alpha = 0.7,
    center = TRUE, errorbar.draw = FALSE, outlier.color = NA) +
  geom_half_point(aes(color=ASRS), size=.4)+
  theme_classic() +
  labs(x = "eBMRQ", y = "m_PLV", color = "ASRS Group") +
  guides(fill="none")


t.test(filter(mean_electrode_data, condition=="Mod")$m_PLV, filter(mean_electrode_data, condition=="Unmod")$m_PLV, alternative="greater",paired=TRUE)

ggplot(mean_electrode_data, aes(x=m_PLV, y=RTCV, color=ASRS)) +
  geom_point() +
  stat_smooth(method = "lm", geom="line",alpha=.8, size=1) +
  theme_minimal() +
  facet_wrap(vars(condition)) +
  scale_color_brewer(palette = "Paired")


ggplot(mean_electrode_data%>%filter(condition=="Mod"), aes(x=m_PLV, y=RTCV, color=eBMRQ.tertile)) +
  geom_point() +
  stat_smooth(method = "lm", geom="line",alpha=.8, size=1) +
  theme_minimal() +
  facet_wrap(vars(ASRS)) +
  scale_color_brewer(palette = "Paired")


ggplot(mean_electrode_data%>%filter(condition=="Mod"), aes(x=m_PLV, y=RTCV)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_minimal() +
  scale_color_brewer(palette = "Paired")


library(lme4)
library(lmerTest)
library(emuR)
library(pracma)



t.test(filter(mean_electrode_data, condition=="Mod")$RTCV, filter(mean_electrode_data, condition=="Unmod")$RTCV, alternative="less",paired=TRUE)
model.out <- lm(RTCV~m_PLV, data=mean_electrode_data)


model.out <- lmer(RTCV~m_PLV*condition + (1|sub.idx), data=mean_electrode_data)
summary(model.out)

# Important GLM Model
model.out <- glm(RTCV~m_PLV,family=Gamma(link = "log"), data=mean_electrode_data%>%filter(condition=="Mod"))

# Not significant indicating that PLV variability over subjects is a predictor
model.out <- glmer(RTCV~m_PLV + (1|sub.idx),family=Gamma(link = "log"), data=mean_electrode_data%>%filter(condition=="Mod"))

# IMPORTANT LMER Model
model.out<-lmer(m_PLV~condition*ASRS*eBMRQ + (1|sub.idx), data=mean_electrode_data)

summary(model.out)
qqnorm(resid(model.out))
qqline(resid(model.out))


library(ggplot2)
library(GGally) # extends ggplot2
ggplot(mean_electrode_data, aes(x=m_PLV,y=RTCV,color=condition)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)
