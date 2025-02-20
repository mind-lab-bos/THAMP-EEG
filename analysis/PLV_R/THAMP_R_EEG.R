"
THAMP_R_EEG.R
Arun Asthagiri
01-21-2025
"


# library(eegUtils)
# plot_butterfly(demo_epochs)
# topoplot(demo_epochs, 
#          time_lim = c(.22, .25 ))

library(R.matlab)


filepath = "/Users/arun/Documents/MINDLab/THAMP/EEG_Data/PLV_R/"
filename = "PLV_to_R.mat"
filename_plv = "PLV_to_R.mat"
data <- readMat(paste(filepath, filename, sep=""))
data_PLV <- readMat(paste(filepath, filename_plv, sep=""))

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

qualtrics_filename<- "qualtrics.csv"
data_qualtrics <- read.csv(paste(filepath, qualtrics_filename,sep=""), header=FALSE)
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

model.out<-lmer(RTCV~m_PLV + (m_PLV|sub.idx), data=mean_electrode_data%>%filter(condition=="Mod"))

model.out <- glmer(RTCV~m_PLV*ASRS+ (1|sub.idx), 
                   family=Gamma(link = "log"), data=mean_electrode_data%>%filter(condition=="Mod"))
summary(model.out)



model.out <- lmer(RTCV~m_PLV + (1|sub.idx), data=mean_electrode_data%>%filter(condition=="Mod"))
qqnorm(resid(model.out))
qqline(resid(model.out))
summary(model.out)

model.out <- lmer(RTCV~m_PLV + (1|sub.idx), data=mean_electrode_data%>%filter(condition=="Unmod"))
summary(model.out)

model.out <- lmer(RT~condition*m_PLV + (condition|sub.idx), data=mean_electrode_data)
summary(model.out)


library(ggplot2)
library(GGally) # extends ggplot2
ggplot(mean_electrode_data, aes(x=m_PLV,y=RTCV,color=condition)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)




######################
filepath = "/Users/arun/Documents/MINDLab/THAMP/EEG_Data/PLV_R/"
filename = "THAMP_EEG_RTCV_RT.mat"
filename_plv = "PLV_by_Song.mat"
data <- readMat(paste(filepath, filename, sep=""))
data_PLV <- readMat(paste(filepath, filename_plv, sep=""))

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

library(ggplot2)
library(ggridges)
ggplot(mean_electrode_data, aes(x=m_PLV, fill=condition)) +
  geom_density(alpha=.5)
  
ggplot(mean_electrode_data, aes(x=condition, y=m_PLV, color=condition)) +
  geom_violin()+
  # geom_point() +
  geom_jitter() +
  geom_line(aes(group = sub.idx), color = "black", alpha = 0.1) 

t.test(filter(mean_electrode_data, condition=="Mod")$m_PLV, filter(mean_electrode_data, condition=="Unmod")$m_PLV, alternative="greater",paired=TRUE)


ggplot(mean_electrode_data, aes(x=m_PLV, y=RTCV, color=condition)) +
  geom_point() +
  stat_smooth(method = "lm") +
  theme_minimal() +
  scale_color_brewer(palette = "Paired")


library(lme4)
library(lmerTest)
library(emuR)
library(pracma)

mean_electrode_data$condition <- relevel(as.factor(mean_electrode_data$condition), ref="Unmod")
# mean_electrode_data$condition <- as.factor(mean_electrode_data$condition)

t.test(filter(mean_electrode_data, condition=="Mod")$RTCV, filter(mean_electrode_data, condition=="Unmod")$RTCV, alternative="less",paired=TRUE)
model.out <- lm(RTCV~m_PLV*condition, data=mean_electrode_data)


model.out <- lmer(RTCV~m_PLV*condition + (1|sub.idx), data=mean_electrode_data)
summary(model.out)
model.out <- glmer(RTCV~m_PLV+ (1|sub.idx), 
                   family=Gamma(link = "log"), data=mean_electrode_data)
summary(model.out)


model.out <- lmer(RTCV~m_PLV + (1|sub.idx), data=mean_electrode_data%>%filter(condition=="Mod"))
qqnorm(resid(model.out))
qqline(resid(model.out))
summary(model.out)

model.out <- lmer(RTCV~m_PLV + (1|sub.idx), data=mean_electrode_data%>%filter(condition=="Unmod"))
summary(model.out)

model.out <- lmer(RT~condition*m_PLV + (condition|sub.idx), data=mean_electrode_data)
summary(model.out)


library(ggplot2)
library(GGally) # extends ggplot2
ggplot(mean_electrode_data, aes(x=m_PLV,y=RTCV,color=condition)) + 
  geom_point()+
  geom_smooth(method='lm', formula= y~x)


library(R.matlab)
filepath = "/Users/arun/Documents/MINDLab/THAMP/EEG_Data/PLV_R/"
filename = "PLV_over_time_to_R.mat"
data <- readMat(paste(filepath, filename, sep=""))

library(tidyr)
library(reshape2)
library(dplyr)

mod_PLV <- melt(data$mod.PLV.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "electrode", "PLV")) %>%
  pivot_wider(names_prefix = "PLV_electrode_",names_from=electrode, values_from=PLV)
mod_RTCV <- melt(data$mod.SART.RTCV.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "RTCV")) 
mod_RT <- melt(data$mod.SART.RT.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "RT")) 
mod_data <- mod_PLV %>% merge(mod_RTCV, by=c("sub.idx", "time", "song.num")) %>%
  merge(mod_RT, by=c("sub.idx", "time", "song.num")) %>%
  mutate(condition = "Mod")

unmod_PLV <- melt(data$unmod.PLV.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "electrode", "PLV")) %>%
  pivot_wider(names_prefix = "PLV_electrode_",names_from=electrode, values_from=PLV) # %>% 
unmod_RTCV <- melt(data$unmod.SART.RTCV.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "RTCV"))
unmod_RT <- melt(data$unmod.SART.RT.over.time) %>% 
  setNames(c("sub.idx","time", "song.num", "RT"))
unmod_data <- unmod_PLV %>% merge(unmod_RTCV, by=c("sub.idx", "time", "song.num")) %>%
  merge(unmod_RT, by=c("sub.idx", "time", "song.num")) %>%
  mutate(condition = "Unmod")

all_data<-rbind(mod_data, unmod_data)

library(lme4)
library(lmerTest)
library(emuR)
library(pracma)

new_data <- all_data %>% pivot_longer(cols = starts_with("PLV_electrode_"), names_prefix = "PLV_electrode_", names_to="electrode", values_to="PLV") %>%
  group_by(sub.idx, song.num, electrode) %>%
  mutate(RTCV.corr = cor(detrend(log(RTCV)), detrend(log(PLV)))^2)


model.out <- lmer(RTCV.corr~condition  + (1|sub.idx:song.num), data=new_data%>%filter(electrode==1))
summary(model.out)

library(ggplot2)


library(GGally) # extends ggplot2
# pairwise exploration of your data columns

ggpairs(data=new_data%>%filter(sub.idx==1), columns = c("PLV", "RTCV", "RT"))


# test_data <- data.frame(t(rbind(c(rep(1, 10),rep(2, 10)), c(1:10,1:10), c(1:10, (11:20)^2))))
# colnames(test_data)<-c("var0", "var1", "var2")
# test_data<-test_data %>% mutate_at(c("var0"), as.factor)
# test_data<-test_data %>% group_by(var0) %>% mutate(var3 = resid(lm(var2~var1)))
# test_data<-test_data %>% group_by(var0) %>% mutate(var3 = detrend(var2))

mean_PLV_data <- all_data %>% 
  pivot_longer(names_to="electrode",cols = starts_with("PLV_electrode_"), names_prefix= "PLV_electrode_",values_to="PLV") %>%
  group_by(sub.idx, song.num, time) %>% 
  mutate(m_PLV = mean(PLV)) %>% 
  select(-c(electrode, PLV)) %>% 
  distinct() %>%
  group_by(sub.idx, song.num) %>% arrange(time) %>%
  mutate(m_PLV = detrend(m_PLV), RTCV =(detrend(RTCV)), RT = detrend(RT)) %>% # detrend all variables
  mutate(m_PLV = m_PLV-min(m_PLV)+10e-5, RTCV =RTCV-min(RTCV)+10e-5) # detrend all variables



model.out <- lmer(RTCV~m_PLV+(m_PLV|sub.idx:song.num), data=mean_PLV_data%>%filter(condition=="Mod"))
model.out <- lmer(RTCV~m_PLV+(m_PLV|sub.idx), data=mean_PLV_data%>%filter(condition=="Mod"))

shapiro.test(resid(model.out))
library(nortest)
ad.test(resid(model.out))
ks.test(resid(model.out), "pnorm")

model.out <- glmer(RTCV~m_PLV+(m_PLV|sub.idx), data=mean_PLV_data%>%filter(condition=="Mod"), family = Gamma(link = "log"))
summary(model.out)
qqnorm(resid(model.out))
qqline(resid(model.out), col = "red")



library(brms)
mean_PLV_data$condition <- relevel(as.factor(mean_PLV_data$condition), ref="Unmod")
full.model.bayes <- brm(RTCV ~ m_PLV*condition + (m_PLV | sub.idx:song.num), 
                   data = mean_PLV_data, family = student())
model.bayes <- brm(RTCV ~ m_PLV+ (m_PLV | sub.idx:song.num), 
                   data = mean_PLV_data%>%filter(condition=="Mod"), family = student())

summary(model.bayes)
qqnorm(residuals(model.bayes))
qqline(residuals(model.bayes))
ad.test(residuals(model.bayes))



library(forcats)
topos <- array(dim=61)
surrogate_size <- 50
data <- mod_data
for (electrode_idx in 1:length(topos)){
  print(electrode_idx)
  lmer_fn <- paste("RTCV~PLV_electrode_", electrode_idx,"+(time|sub.idx)+(PLV_electrode_", electrode_idx,"|sub.idx)", sep="")
  model.out <- lmer(lmer_fn, data=data)
  true_value <-summary(model.out)$coefficients[2,4]
  topos[electrode_idx] <- true_value
  # surr_values <- array(dim=surrogate_size)
  # for (surr_idx in 1:surrogate_size){
  #   # print(surr_idx)
  #   shuffled_data <- data%>% mutate(sub.idx = sample(sub.idx))# randomly shuffle subjects
  #   model.out <- lmer(lmer_fn, data=shuffled_data)
  #   surr_values[surr_idx] <- summary(model.out)$coefficients[2,4]
  # }
  # topos[electrode_idx] <- (topos[electrode_idx] -mean(surr_values))/sd(surr_values)
  # 
  
  
  # p <- summary(model.out)$coefficients[2, 5]
  # true_z_value <- ifelse(summary(model.out)$coefficients[2,4] >0, qnorm(1 - p / 2), -qnorm(1 - p / 2))
  # random_lag = sample(data$time, surrogate_size, replace=TRUE)
  # 
  # surr_values <- array(dim=surrogate_size)
  # for (surr_idx in 1:surrogate_size){
  #   shuffled_data <- data %>%
  #     group_by(sub.idx, song.num) %>%
  #     arrange(sub.idx, song.num, time) %>%
  #     mutate(shuffled_RT = shift(RT,random_lag[surr_idx]-1)) %>%
  #     mutate(shuffled_RTCV = shift(RTCV,random_lag[surr_idx]-1))
  #   lmer_fn <- paste("shuffled_RTCV~PLV_electrode_", electrode_idx,"+time+(PLV_electrode_", electrode_idx,"|sub.idx:song.num)", sep="")
  #   model.out <- lmer(lmer_fn, data=shuffled_data)
  #   p <- summary(model.out)$coefficients[2, 5]
  #   surr_z_value <- ifelse(summary(model.out)$coefficients[2,4] >0, qnorm(1 - p / 2), -qnorm(1 - p / 2))
  #   surr_values[surr_idx] <- surr_z_value
  # }
  # topos[electrode_idx] <- (true_z_value - mean(surr_values))/sd(surr_values)
}

writeMat(paste(filepath,"topos_out.mat", sep=""), topos=topos)
# writeMat(paste(filepath,"topos_unmod_RTCV.mat", sep=""), topos=topos)
# 
# for (idx in 1:length(topos)){
#   lmer_fn <- paste("RT~PLV_electrode_", idx,"+time+(1|sub.idx/song.num)", sep="")
#   model.out <- lmer(lmer_fn, data=mod_data)
#   topos[idx] <- summary(model.out)$coefficients[2, 4]
# }
# 
# summary(model.out)
# summary(model.out)$coefficients[2, 4]

