# tone.R
# Ola Ozernov-Palchik 
# oozernov@mit.edu
# Updated-11 2 2020
#
# MIT Speech Perception study, Tone Anchoring Task
#

#### Setup ####
dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd("dir")

Packages <- c("dplyr", "readr", "magrittr", "tidyr", "ggplot2", "lme4", "lmerTest",
              "emmeans", "sjstats","dabestr","gridExtra")
lapply(Packages, library, character.only = TRUE)

#### Load and organize the data ####
# adult data
adult <- read_csv("data/tonethres_data_a_010719.csv")
groups_a <- read_csv("data/adult_groups_012418.csv") %>%
    mutate(group = as.factor(group)) %>%
    rename(PartID = Subject)    

adult %<>%
    rename(PartID = Subject) %>%
    inner_join(groups_a, by = c("PartID")) %>% # add DD columns to the data
    select(-X1)

# child data
child <- read_csv("data/child_tonethres_data_122718.csv")
groups_c <- read.csv("data/groups_041817.csv") %>%
    dplyr::mutate(group = ifelse(Typical == 1, "Typ", "Dys")) %>%
    dplyr::select(-DD, -Typical)

child %<>%
    inner_join(groups_c, by = "PartID") %>% # add DD columns to the data
    mutate(group = as.factor(group)) %>%
    dplyr::select(-Subject)

# combine the datasets
d <- bind_rows("Adult" = adult,
               "Child" = child,
               .id = "age") %>%
    mutate(com_cond = paste(group, "-", age),
           # rename condition to standard/no-standard
           cond = ifelse(cond == "a", "Standard", "No-standard"))

groups <- dplyr::bind_rows("Adult" = groups_a,
                    "Child" = groups_c,
                    .id = "age")

# exclude outlier participant
#d = d %>% filter(PartID != 5174)

# convert hz to cent function
hz_to_cents <- function(a, b) {
    abs(1200 * log2(b/a))
}

# add freqDiff in cents
d$cents <- hz_to_cents(d$lowFreq, d$highFreq)

#### Get sample sizes by group ####
counts <- d %>%
    group_by(PartID, age) %>%
    summarize(m = mean(cents)) %>%
    ungroup() %>%
    left_join(groups, by = c("PartID", "age"))
count(counts, group, age)

#### Analysis ####

#### Percent Accuracy ####

#### Plot ####
d_correct <- d %>%
    group_by(cond, com_cond, group) %>%
    summarize(m_correct = mean(correct, na.rm = TRUE),
              se = plotrix::std.error(correct, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(upper = m_correct + 1.96*se,
           lower = m_correct - 1.96*se)

ggplot(d_correct, aes(x = cond, y = m_correct*100, fill = cond)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    scale_fill_grey() +
    geom_errorbar(aes(ymin = lower*100, ymax = upper*100, width = .2), position = position_dodge(.9)) +
    labs(title = "Tone Threshold Accuracy", x = "Condition", y = "% Correct") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ com_cond, ncol = 2)

#### Analysis ####

d_2 <- d %>%
    group_by(PartID, cond, group, age) %>%
    summarize(m_correct = mean(correct, na.rm = TRUE))
LME_model1 <- lmer(m_correct ~ cond*group*age + (1|PartID), data = d_2, REML = FALSE)
anova(LME_model1) #main effect for accuracy by cond, group, age

#### JND ####
# Calculate
# mean frequency difference in last seven reversals (9-16)

d_reversals <- d %>%
    group_by(PartID)  %>% 
    summarize(max_rev=max(reversals),rev_num=max_rev-7) %>% 
    left_join(d, by = c("PartID"))%>% 
    filter(reversals >= rev_num) 

    #mutate(cents = hz_to_cents(lowFreq, highFreq))
d_idthres <- d_reversals %>%
    group_by(cond, com_cond, group) %>%
    summarize(m_fd = mean(freqDiff))

d_jnd_hz <- d_reversals %>%
    mutate(cond = ifelse(cond == "No-standard", "NS", "S"))%>%
    group_by(cond, com_cond, group) %>%
    summarize(m_fd = mean(freqDiff),
              se = plotrix::std.error(freqDiff)) %>%
    mutate(upper = m_fd + 1.96 * se,
           lower = m_fd - 1.96 * se)

##Individual performance
d_ind_jnd_hz <- d_reversals %>%
    mutate(cond = ifelse(cond == "No-standard", "NS", "S"))%>%
    group_by(PartID, cond, group, age) %>%
    summarize(jnd = mean(freqDiff, na.rm = T))

#plot individual slopes by age
ggplot(data = d_ind_jnd_hz, aes(x = cond, y = jnd, group = PartID, color = group)) +
    geom_point() +
    geom_line() +
    labs(x = 'Condition', y = "JND (Hz)", title = "Individual JND in Hz") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~age, ncol = 2)

####  Jnd effects ####

#Plot effects
ggplot(d_jnd_hz, aes(x = cond, y = m_fd, fill = cond)) +
    geom_bar(stat = 'identity', position = position_dodge()) +
    scale_fill_grey() +
    geom_errorbar(aes(ymin = lower, ymax = upper, width = .2), position = position_dodge(.9)) +
    labs(title = "JND in Hz", x = "Condition", y = "JND in Hertz (Hz)") +
    theme(plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
    facet_wrap(~com_cond, ncol = 2)

d_ind_jnd_hz$com_cond <- paste(d_ind_jnd_hz$group, "-", d_ind_jnd_hz$age)
d_ind_jnd_hz$group_cond = paste(d_ind_jnd_hz$group, "-", d_ind_jnd_hz$cond)

multi.two.group.unpaired <- 
    d_ind_jnd_hz %>%
    dabest(group_cond,jnd,
           idx = list(c("Dys - S", "Dys - NS"), 
                      c("Typ - S", "Typ - NS")),
           paired = FALSE
    )
p1<-plot(multi.two.group.unpaired,
     rawplot.ylabel = "JND",color.column = age,
     effsize.ylabel = "Effect Size")

#Analyze effects
jnd_model <- lmer(jnd ~ cond*group*age + (1|PartID), data = d_ind_jnd_hz,REML = TRUE)
rand(jnd_model) #test the random effect in the model
anova(jnd_model)

# ###Post hoc comparisons
lsmeans(jnd_model, list(pairwise ~ group|age), adjust = "tukey")
#lsmeans(jnd_model, list(pairwise ~ group|age), adjust = "tukey")

#### NTD ####
# Calculations

d_s_td = d_ind_jnd_hz %>% 
    filter(cond == 'S') %>% 
    group_by(PartID) %>% 
    summarise(m_fd_hz = mean(jnd))

d_ns_td = d_ind_jnd_hz %>% 
    filter(cond == 'NS')%>% 
    group_by(PartID) %>% 
    summarise(m_fd_hz = mean(jnd))

#remove participants who don't have both conditions
d_extra = d_ns_td %>% filter(!(PartID %in% d_s_td$PartID))
d_ns_td = d_ns_td %>% filter(!(PartID %in% d_extra$PartID))

d_ntd = as.data.frame((d_s_td$m_fd_hz - d_ns_td$m_fd_hz )/(d_s_td$m_fd_hz + d_ns_td$m_fd_hz))
names(d_ntd)[names(d_ntd) == '(d_s_td$m_fd_hz - d_ns_td$m_fd_hz)/(d_s_td$m_fd_hz + d_ns_td$m_fd_hz)'] <- 'm_fd_hz'
d_ntd$PartID = d_s_td$PartID

# now add risk columns and split into groups, then put back together for plot
d_ntd = merge(d_ntd, groups, "PartID")
d_ntd$com_cond = paste(d_ntd$group, "-", d_ntd$age)

#Plot effect size
#https://cran.r-project.org/web/packages/dabestr/vignettes/using-dabestr.html

multi.two.group.unpaired <- 
    d_ntd %>%
    dabest(com_cond,m_fd_hz,
           idx = list(c("Dys - Adult", "Typ - Adult"), 
                      c("Dys - Child", "Typ - Child")),
           paired = FALSE
    )
p2<-plot(multi.two.group.unpaired,
     rawplot.ylim = c(-1, 0.7),color.column = group,
     rawplot.ylabel = "NTD",
     effsize.ylabel = "Effect Size")

pdf("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT/jnd_ntd_plots.pdf") # Open a new pdf file
grid.arrange(p1, p2, nrow=2) # Write the grid.arrange in the file
dev.off() # Close the file

#### Analysis ####
model_ntd <- lm(m_fd_hz ~ group*age, data = d_ntd)
anova(model_ntd)
eta_sq(model_ntd)
lsmeans(model_ntd, list(pairwise ~ group|age), adjust = "tukey")

#### Frequency difference by trial ####
d_fd_trial_hz <- d %>%
    group_by(trialNum, group, com_cond, cond) %>%
    summarize(m_fd = mean(freqDiff),
              se = plotrix::std.error(freqDiff)) %>%
    ungroup() %>%
    mutate(upper = m_fd + 1.96 * se,
           lower = m_fd - 1.96 * se)

p1<-ggplot(d_fd_trial_hz, aes(colour = cond, y = m_fd, x = trialNum)) +
    scale_color_manual(values = c('red', 'blue')) +
    geom_line() + geom_errorbar(aes(ymin = lower, ymax = upper, width = 1)) +
    labs(x = "Trial number", y = "Frequency Difference (Hz)", title = "Frequency Difference by Trial") +
    theme(legend.title = element_blank(),plot.title = element_blank()+
              scale_color_brewer(palette="Dark2"))+
    facet_wrap(~com_cond, ncol = 2)

new_fd_hz <- d_fd_trial_hz %>%
    select(trialNum, cond, com_cond, m_fd) %>%
    mutate(cond = ifelse(cond == "No-standard", "ns", "s")) %>%
    spread(cond, m_fd) %>%
    mutate(ns_diff_s = ns - s)

p2<-ggplot(new_fd_hz, aes(colour = com_cond, y = ns_diff_s, x = trialNum)) +
    geom_point() + stat_smooth(method = "loess") +
    labs(x = "Trial number", y = "Frequency Difference (NS - S) in Hz") +
    theme(legend.title=element_blank(),plot.title = element_blank()+scale_color_brewer(palette="Dark2"))

pdf("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT/freq_by_trial.pdf") # Open a new pdf file
grid.arrange(p1,p2, nrow=2)
grid.arrange(
    grobs = gl,
    widths = c(2, 1),
    layout_matrix = rbind(c(1, 2))
)
dev.off() # Close the file
#########################
#fit first degree polynomial equation:
fit  <- lm(new_fd_hz$ns_diff_s~new_fd_hz$trialNum)
#second degree
fit2 <- lm(new_fd_hz$ns_diff_s~poly(new_fd_hz$trialNum,2,raw=TRUE))
#third degree
fit3 <- lm(new_fd_hz$ns_diff_s~poly(new_fd_hz$trialNum,3,raw=TRUE))
#fourth degree
fit4 <- lm(new_fd_hz$ns_diff_s~poly(new_fd_hz$trialNum,raw=TRUE))
#generate range of 50 numbers starting from 30 and ending at 160
xx <- seq(0,70, length=280)
plot(new_fd_hz$trialNum,new_fd_hz$ns_diff_s,pch=19,ylim=c(0,150))
lines(xx, predict(fit, data.frame(x=xx)), col="red")
lines(xx, predict(fit2, data.frame(x=xx)), col="green")
lines(xx, predict(fit3, data.frame(x=xx)), col="blue")
lines(xx, predict(fit4, data.frame(x=xx)), col="purple")

#########
names(d)
library(reshape2)
data_wide <- dcast(d, PartID+trialNum + age+group ~ cond, value.var="freqDiff",fun.aggregate = mean)
data_wide$diff<-(data_wide$Standard-data_wide$`No-standard`)/(data_wide$Standard+data_wide$`No-standard`)

#### Analysis ####
ids=unique(data_wide$PartID)
#data_wide<-data_wide%>%dplyr::filter(trialNum<30)

d_glm_fit_list <- vector(mode = "list", length = nrow(data_wide)*2)
index <- 0
for (subj in ids) {
        index <- index + 1
        print(paste(index, subj))
        d_subj <- data_wide %>%
            filter(PartID == subj)
        subj_model<-lm(d_subj$diff~poly(d_subj$trialNum,3,raw=TRUE),na.action=na.exclude)
        b <- as.numeric(subj_model$coefficients[1]) # coefficient intercept
        m <- as.numeric(subj_model$coefficients[2]) # coefficient slope
        dev_fit <- as.numeric(deviance(subj_model))
        d_glm_fit_list[[index]] <- list(
            PartID = subj,
            group = filter(groups, PartID == subj)$group,
            age = filter(groups, PartID == subj)$age,
            slope_coef = abs(m),
            deviance = dev_fit
            )}

d_glm_fit <- bind_rows(d_glm_fit_list)
d_glm_fit$age<-as.factor(d_glm_fit$age)
d_glm_fit$group<-as.factor(d_glm_fit$group)

m1<-lm(slope_coef~group+age+group:age,weights = deviance,data=d_glm_fit)
anova(m1)
lsmeans(m1, list(pairwise ~ group), adjust = "tukey")
lsmeans(m1, list(pairwise ~ age), adjust = "tukey")

#### Summary of results ####
# #create d_export (for individual differences analysis)
# d_jnd <- spread(d_ind_jnd_hz, cond, jnd)
# d_export <- d_ntd %>%
#     select(PartID, com_cond, m_fd_hz, m_ntd_hz)
# write.csv(d_export, "~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/dys_comp_tonethres_030419.csv")
