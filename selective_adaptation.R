# tone.R
# Ola Ozernov-Palchik 
# oozernov@mit.edu
# Updated-11 2 2020
#
# MIT Speech Perception study, Selective Adaptation Task
#

#### Setup ####
dir<-dirname(rstudioapi::getActiveDocumentContext()$path)
setwd("dir")
#
Packages <- c("dplyr", "readr", "magrittr", "tidyr", "ggplot2", 
              "lme4", "lmerTest","emmeans", "sjstats", "plotrix","dabestr","lmerTest",
              "lmPerm","gridExtra", "grid")
lapply(Packages, library, character.only = TRUE)

#### Load and organize the data ####
# adult data
d_a <- read_csv("data/phonselec_data_a_010719.csv")
groups_a <- read_csv("data/adult_groups_012418.csv") %>%
    dplyr::rename(PartID = Subject) # rename Subject column

d_a %<>%
    dplyr::rename(PartID = Subject) %>% # rename Subject column
    inner_join(groups_a, by = "PartID") %>% # add DD groups
    dplyr::select(-X1) # drop X1 column

# child data
d_c <- read_csv("data/phonselec_data_122718.csv")
groups_c <- read_csv("data/groups_041817.csv") %>%
    mutate(group = ifelse(Typical == 1, "Typ", "Dys")) # label groups

d_c %<>% inner_join(groups_c, by = "PartID") %>% # add DD groups
    dplyr::select(-Subject, -DD, -Typical) # drop Subject, DD, Typical columns

# combined data
d <- bind_rows("Adult" = d_a,  # combine adult and child datasets
               "Child" = d_c,
               .id = "age") %>%
    # exclude participants who had <0 for end d prime
    filter(!PartID %in% c('READER_002', 'READ_6103','ABCD_1776')) %>%

    mutate_at(vars(group, age), as.factor) %>% # convert group and age to factors

    # get step column from soundfile column
    separate(soundfile, into = c("x", "step", "y"), sep = "_", remove = FALSE) %>%
    separate(step, into = c("x", "step"), sep = "d", remove = TRUE) %>%
    dplyr::select(-x, -y) %>% # only keep step column from these split columns
    mutate(step = as.numeric(step), # convert step to numeric
           # add column for whether response was /b/ (ignores NAs)
           b = ifelse(is.na(response), NA, ifelse(response == 'b', 1, 0)))

groups <- bind_rows("Adult" = groups_a,
                    "Child" = groups_c,
                    .id = "age") %>%
    dplyr::select(PartID, group, age) %>%
    mutate_at(vars(group, age), as.factor) # convert group and age to factors

#### Display individual plots ####
ind_d <- d %>%
    group_by(PartID, step, adaptor, age) %>%
    dplyr::summarize(prop_b = sum(b)/4) %>%
    left_join(groups, by = c("PartID", "age"))

ind_d_subs <- ind_d %>%
    group_by(PartID, adaptor, group, age) %>%
    dplyr::summarize(m = mean(prop_b)) 

# plot
ind_d_subs_plot <- ggplot(data = ind_d_subs,
                          aes(x = adaptor, y = m, color = group, group = PartID)) +
    geom_point() +
    geom_line() +
    facet_wrap( ~ age) +
    labs(x = 'Adaptor', y = "mean /b/ response",
         title = "Individual Performance (Mean /b/ response by Condition)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.title = element_text(hjust = 0.5))
ind_d_subs_plot

#### Subject count by group ####
counts <- d %>%
    group_by(PartID, age) %>%
    summarize(m = mean(as.double(step))) %>%
    ungroup() %>%
    left_join(groups, by = c("PartID", "age")) %>%
    count(group, age)
counts
#### Create d export ####
ind_d_subs_b <- ind_d_subs %>%
    filter(adaptor == "b") %>%
    rename(mean_b_under_b_adaptor = m)
d_export <- select(ungroup(ind_d_subs_b), PartID)

#### Run logistic regression models to extract slope and inflection ####
d_glm_fit_list <- vector(mode = "list", length = nrow(d_export)*2)
index <- 0
for (subj in d_export$PartID) {
    for (a in c("d", "b")) {
        index <- index + 1
        #print(paste(index, subj, a))
        d_subj <- ind_d %>%
            filter(adaptor == a, PartID == subj)
        subj_model <- glm(prop_b ~ step,
                          data = d_subj)
        b <- as.numeric(subj_model$coefficients[1]) # coefficient intercept
        m <- as.numeric(subj_model$coefficients[2]) # coefficient slope
        i <- ifelse(a == "d", -b/m, i) # only recalculate step for "d"
        dev_fit <- as.numeric(deviance(subj_model))
        value <- as.numeric(predict(subj_model, list(step = i), type = "response"))

        d_glm_fit_list[[index]] <- list(
            PartID = subj,
            adaptor = a,
            group = filter(groups, PartID == subj)$group,
            age = filter(groups, PartID == subj)$age,
            slope_coef = m,
            inflection_step = -b/m, # recalculate for "d" and "b"
            value_at_da_i = value,
            deviance = dev_fit
        )
    }
}
d_glm_fit <- bind_rows(d_glm_fit_list)
str(d_glm_fit)
d_glm_fit$adaptor<-as.factor(d_glm_fit$adaptor)


#### Slope transformation ####
d_glm_fit$slope_coef_tran <- d_glm_fit$slope_coef * -1
min_slope <- min(d_glm_fit$slope_coef_tran)
d_glm_fit$slope_coef_tran <- (d_glm_fit$slope_coef_tran + abs(min_slope)) + 1
d_glm_fit$slope_coef_tran <- log10(d_glm_fit$slope_coef_tran)

#### Plot distribution ####
ggplot(d_glm_fit, aes(x = slope_coef_tran, fill = group)) +
    geom_density(alpha = .5) +
    labs(title = "Transformed Slope Distribution") +
    theme(plot.title = element_text(hjust = 0.5)) +
    facet_wrap(~ age)

#### Main effects of slope ####
#with adaptor effect
LME_model1 <- lmer(slope_coef_tran ~ adaptor*group*age + (1|PartID),
                   weights = deviance,
                   data = d_glm_fit,
                   REML = TRUE)
anova(LME_model1)
eta_sq(LME_model1)
#### Post hoc comparisons ####
lsmeans(LME_model1, list(pairwise ~ group), adjust = "tukey")
lsmeans(LME_model1, list(pairwise ~ group|age), adjust = "tukey")
lsmeans(LME_model1, list(pairwise ~ adaptor), adjust = "tukey")

#across both adaptors
d_glm_fit2 = d_glm_fit %>%
    dplyr::group_by(PartID)%>%
    dplyr::summarise(mean_slope = mean(slope_coef_tran, na.rm = T))%>%
    left_join(groups, by = c("PartID"))

m3 <- lm(mean_slope ~ group*age, data = d_glm_fit2)
anova(m3)
lsmeans(m3, list(pairwise ~ group|age), adjust = "tukey")
eta_sq(m3)


#Plot slope effects
##All together
# d_glm_fit_grp<-d_glm_fit%>%
#     mutate(com_cond = paste(d_glm_fit$group, "-", d_glm_fit$age))%>%
#     select(PartID, com_cond,adaptor,slope_coef_tran) %>%
#     group_by(PartID,com_cond) %>%
#     summarise(mean_slope = mean(as.double(slope_coef_tran)))
# 
# multi.two.group.unpaired <- 
#     d_glm_fit_grp %>%
#     dabest(com_cond,mean_slope,
#            idx = list(c("Dys - Adult", "Typ - Adult"), 
#                       c("Dys - Child", "Typ - Child")),
#            paired = FALSE
#     )
# plot(multi.two.group.unpaired,
#               #rawplot.ylim = c(-1, 0.7),
#               rawplot.ylabel = "Slope",
#               effsize.ylabel = "Effect Size")
##By age seperately
d_glm_fit_grp<-d_glm_fit%>%
    mutate(com_cond = paste(d_glm_fit$group, "-", d_glm_fit$age))%>%
    select(PartID, age,group,adaptor,slope_coef_tran) %>%
    group_by(PartID,age, group) %>%
    summarise(mean_slope = mean(as.double(slope_coef_tran)))

multi.two.group.unpaired_a <- 
    d_glm_fit_grp %>%
    filter(age=="Adult")%>%
    dabest(group,mean_slope,
           idx = list(c("Dys", "Typ")),
           paired = FALSE
    )
a<-plot(multi.two.group.unpaired_a,
    rawplot.ylim = c(0, 0.04),
     rawplot.ylabel = "Adult Slope",
     effsize.ylabel = "Effect Size")

multi.two.group.unpaired_c <- 
    d_glm_fit_grp %>%
    filter(age=="Child")%>%
    dabest(group,mean_slope,
           idx = list(c("Dys", "Typ")),
           paired = FALSE
    )
c<-plot(multi.two.group.unpaired_c,
        rawplot.ylim = c(0, 0.04),
        rawplot.ylabel = "Child Slope",
        effsize.ylabel = "Effect Size")

grid.arrange(a,c)
#### Inflection ####????WHAT IS HAPPENING?
d_glm_fit_na<-na.omit(d_glm_fit)
LME_model_i <- lmer(inflection_step ~ group*age*adaptor + (1|PartID),
                    weights = deviance,
                    data = d_glm_fit_na,
                    REML = TRUE)
rand(LME_model_i)
anova(LME_model_i)
lsmeans(LME_model_i, list(pairwise ~ adaptor), adjust = "tukey")
lsmeans(LME_model_i, list(pairwise ~ adaptor|age|group), adjust = "tukey")
eta_sq(LME_model_i)

#### Calculate inflection differences ####
adaptor_diff <- d_glm_fit %>%
    dplyr::select(PartID, adaptor, inflection_step) %>%
    spread(adaptor, inflection_step) %>%
    dplyr::mutate(diff = d - b)

group_inflection <- d_glm_fit %>%
    filter(adaptor == "d") %>%
    left_join(adaptor_diff, by = "PartID") %>%
    select(PartID,
           da_i = slope_coef,
           adaptor_diff = diff,
           group,
           age,
           deviance)

#### Main Effects for Inflection ####

#### Plot group effects ####
selec_models <- read_csv("data/selec_models.csv")
raw_df_list = vector("list", length = nrow(selec_models))
predict_df_list = vector("list", length = nrow(selec_models))
selec_models$age[selec_models$age=="adult"]<-'Adult' #rename variables
selec_models$age[selec_models$age=="child"]<-"Child"


for (i in 1:nrow(selec_models)) {
    temp_params <- selec_models[i,]
    temp_pos <- d %>%
        filter(adaptor == temp_params$adaptor,
               group == temp_params$group,
               age == temp_params$age) %>%
        group_by(step) %>%
        summarize(prop_b = sum(b)/4)
    temp_model <- glm(prop_b ~ step, family = binomial(link = "logit"), data = temp_pos)
    new_step <- seq(min(temp_pos$step), max(temp_pos$step), len = 100)
    temp_data <- data.frame(step = new_step,
                            value = predict(temp_model, list(step = new_step),
                                            type = "response"),
                            adaptor = temp_params$adaptor,
                            group = temp_params$group,
                            age = temp_params$age) %>%
        
        mutate(group_age = paste(group, "-", age),
               interact = interaction(adaptor, group, age)) %>%
        mutate_if(is.factor, as.character)
    temp_pos %<>% rename(value = prop_b) %>%
        mutate(adaptor = temp_params$adaptor,
               group = temp_params$group,
               age = temp_params$age,
               group_age = paste(group,"-",age),
               interact = interaction(adaptor, group, age)) %>%
        mutate_if(is.factor, as.character)

    raw_df_list[[i]] <- temp_pos
    predict_df_list[[i]] <- temp_data
}
raw_df <- bind_rows(raw_df_list)
predict_df <- bind_rows(predict_df_list)

group_effects_plot <- ggplot() +
    geom_line(data = predict_df,
              aes(x = step, y = value, group = adaptor, color = adaptor),
              size = 1.05)  +
    geom_hline(yintercept = 0.5, color = 'black') +
    geom_line(data = raw_df,
        aes(x = step, y = value, color = adaptor),
        linetype = "dotted",
        size = 1) +
    geom_point(data = raw_df,
        aes(x = step, y = value, color = adaptor),
        size = 1.5,
        shape = 1) +
    facet_wrap(~ group_age, nrow = 2, scales = "free") + # free scales to repeat axes
    scale_y_continuous(limits = c(0, 1.0)) +
    scale_x_continuous(breaks = seq.int(2, 18, 2)) +
    labs(x = 'Step', y = "Proportion /b/ Response") +
    theme_bw() +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        #legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12)
    )
group_effects_plot +scale_color_brewer(palette="Dark2")+ 
    theme(legend.position="bottom",legend.text = element_text(size=14),legend.title = element_blank())
                                                                                                                  

#### Adaptation effect ####
d$com_cond <- paste(d$group, "-", d$age)
d_conds <- d %>%
   dplyr:: select(PartID, group, age, com_cond)

adaptor_diff <- d %>%
    group_by(adaptor, PartID, step, com_cond, group, age) %>%
    summarize(mean_b = mean(b, na.rm = TRUE)) %>%
    group_by(PartID) %>%
    spread(adaptor, mean_b) %>%
    ungroup() %>%
    mutate(diff = d - b) %>%
    select(PartID, step, diff)


adaptor_diff_pos <- adaptor_diff %>%
    group_by(PartID) %>%
    summarize(mean_diff = mean(diff, na.rm = TRUE)) %>%
    left_join(d_conds, by = "PartID") %>%
    unique()

#### Plot adaptation effects ####
adaptor_diff %<>% left_join(d_conds, by = "PartID") %>%
    unique()

adaptor_diff_total <- adaptor_diff %>%
    dplyr::group_by(step, group, age, com_cond) %>%
    dplyr::summarize(mean_diff = mean(diff, na.rm = TRUE))

adaptor_diff_total$group2 <- factor(adaptor_diff_total$group, levels = c("Dys", "Typ"))
adaptor_diff_total$step <- factor(adaptor_diff_total$step, levels = seq.int(2, 18, 2))

#plot with confidence interval 
adaptor_diff_total_2 <- adaptor_diff %>%
    group_by(step, group, age, com_cond) %>%
    summarize(mean.diff = mean(diff, na.rm = TRUE),
              sd.diff = sd(diff, na.rm = TRUE),
              n.diff = n()) %>%
    ungroup() %>%
    mutate(se.diff = sd.diff / sqrt(n.diff),
           lower.ci.diff = mean.diff - qt(1 - (0.05 / 2), n.diff - 1) * se.diff,
           upper.ci.diff = mean.diff + qt(1 - (0.05 / 2), n.diff - 1) * se.diff) %>%
    mutate_at(vars(step, group), as.factor)

ggplot(data = adaptor_diff_total_2,
       aes(x = step, y = mean.diff, ymin = lower.ci.diff, ymax = upper.ci.diff,
           group = group, fill = group)) +
    geom_ribbon(alpha = 0.5) +
    geom_point() +
    geom_line() +
    labs(x = 'Step', y = "Adaptor Difference", title = element_blank()) +
    theme_bw() +
    theme(text = element_text(size = 18),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.background = element_rect(fill = 'white'),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    facet_wrap( ~ as.factor(age), nrow = 2)+theme(
        strip.text = element_text(size=14),
        legend.text = element_text(size=18),
        legend.title = element_blank())
    


#### Analysis by step ####
str(adaptor_diff)
model_step<-lm(diff ~ step*group*age, data = adaptor_diff)
#model_step <- lmer(diff ~ group*age*step + (1|PartID), data = adaptor_diff, REML = FALSE) #this model is overfitted
summary(model_step)
# post hocs
lsmeans(model_step, list(pairwise ~ step), adjust = "tukey")
#lsmeans(model_step, list(pairwise ~ age|step), adjust = "tukey")


###Difference at step 8
d_8 <- adaptor_diff %>%
    filter(step == "8") %>%
    select(PartID, diff,group,age)
model_8<-lm(diff~group*age, data=d_8)
anova(model_8)
eta_sq(model_8)


#### Test for overall accuracy (d-primes) ####
d_prime <- read.csv("~/Dropbox (MIT)/Com_Dys_2016_data/final_for_sharing/final_code_JT/data/d_prime_123018.csv") %>%
    group_by(PartID, group, age) %>%
    summarize(d = mean(d, na.rm = TRUE))
d_prime2<-d_prime%>%filter(d>0)
model1 <- lm(d ~ group*age, data = d_prime2)
anova(model1)
lsmeans(model1, list(pairwise ~ group|age), adjust = "tukey")
lsmeans(model1, list(pairwise ~ age|group), adjust = "tukey")

aggregate(d ~ group + age, data = d_prime, FUN = "mean")
eta_sq(model1)

#### d_export ####
#adaptor_diff_pos
#d_prime
#step10
d_8 <- adaptor_diff %>%
    filter(step == "8") %>%
    select(PartID, diff)
slope <- d_glm_fit %>%
    group_by(PartID) %>%
    summarize(slope = mean(slope_coef, na.rm = TRUE))
d_export <- left_join(d_10, slope, by = "PartID") %>%
    left_join(adaptor_diff_pos, by = "PartID") %>%
    select(PartID,
           diff_10 = diff,
           slope,
           adapt = mean_diff,
           condition = com_cond)
write.csv(d_export, "dys_comp_selec_adapt_011819.csv", row.names = FALSE)
#ABCD_1776 included?
d_ends<-d%>%filter(step==2|step==18)
d_ends$correct<-ifelse(d_ends$step=='2' & d_ends$response=='b',1,ifelse(d_ends$step=='18' & d_ends$response=='d',1,0))
d_ends2<-d_ends%>% group_by(PartID) %>%
    summarize(accuracy = sum(correct, na.rm = TRUE))%>%
    left_join(groups, by = "PartID") 
m2<-lm(accuracy~age*group,data=d_ends2)
anova(m2)
eta_sq(m2)

