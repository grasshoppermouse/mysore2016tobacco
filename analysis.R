#+ message=F, warning=F

# Load libraries --------------------------------------------------------

# Note: To preserve anonymity, identifying info, e.g., age, 
# has been aggregated in placektobaccopublic
# Results using placektobaccopublic will therefore differ
# slightly from placektobacco

library(placektobaccopublic)
# library(placektobacco) # Not publicly available

library(tidyverse)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(MASS)
library(mixtools)
library(effects)
library(ggalt)
library(ggmosaic)
library(viridis)
library(car)
library(boot)
library(visreg)
library(glmnetUtils)
library(broom)
library(broom.mixed)
library(hagenutils) # for scale_colour_binary, scale_fill_binary

set.seed(5678765)

# Estimate a cotinine cutoff from the data -----------------------------------------------

# Fit Gaussian mixture model
m_mix <- normalmixEM(log10(samples$cotinine))
samples$comp1 <- m_mix$posterior[,1]
samples$comp2 <- m_mix$posterior[,2]

mu <- m_mix$mu
sigman <- m_mix$sigma

cutoff <- 6.5 # Point at which there is an equal probability of belonging to either component
p_cot_mix <-
  ggplot(samples, aes(cotinine, comp1)) + 
  geom_line(colour='red') + 
  geom_line(aes(y=comp2), colour='green') +
  geom_vline(xintercept = cutoff, linetype='dotted') +
  scale_x_log10() +
  labs(x = '\nCotinine (ng/ml)', y = 'Probability\n') +
  theme_bw(15)
p_cot_mix

# Prepare data for regression modeling ------------------------------------

# Combine baseline vars with baseline and followup cotinine values

women_vars <- c('ID', 'GroupID', 'age', 'Presentation', 'pregnancy_status', 'trimester', 'nicotine1', 'nicotine2')

# ID 202 provided followup saliva but did not attend presentation; remove
df <-
  samples %>%
  dplyr::filter(ID != '202') %>% 
  dplyr::select(ID, pre_post, cotinine) %>%
  spread(key = pre_post, value = cotinine) %>%
  left_join(women[women_vars], by = 'ID') %>%
  mutate(
    pregnancy_status = factor(pregnancy_status), 
    Pregnant = ifelse(pregnancy_status == 'P', 'Yes', 'No'),
    Presentation = factor(Presentation),
    baseline_user = pre > 3,
    followup_user = post > 3,
    tobacco_user = pre > 3 & post > 3,
    tobacco_user2 = pre > 3 | post > 3,
    non_user = pre <= 3 & post <= 3,
    nicotine1 = ifelse(is.na(nicotine1), 0, nicotine1),
    nicotine2 = ifelse(is.na(nicotine2), 0, nicotine2),
    "Tobacco use 12hrs (baseline)" = ifelse(nicotine1 == 0, 'No', 'Yes'),
    "Tobacco use 12hrs (followup)" = ifelse(nicotine2 == 0, 'No', 'Yes'),
    Underreporting = (nicotine1 == 0 & pre >= 3) | (nicotine2 == 0 & post >= 3),
    pre_misreport = nicotine1 == 0 & pre >= 3,
    post_misreport = nicotine2 == 0 & post >= 3,
    Overreporting = (nicotine1 == 1 & pre < 3) | (nicotine2 == 1 & post < 3)
    )

# ggplot likes data in long format
df_long <- 
  df %>%
  dplyr::select(ID, pre, post, `Tobacco use 12hrs (baseline)`, `Tobacco use 12hrs (followup)`) %>% 
  rename(
    cotinine_pre = pre,
    cotinine_post = post,
    use_pre = `Tobacco use 12hrs (baseline)`,
    use_post = `Tobacco use 12hrs (followup)`
  ) %>% 
  pivot_longer(-ID, names_to = c('.value', 'Time'), names_sep = '_') %>% 
  rename(
    `Self-reported use` = use
  )

# binwidth Freedman-Diaconis
bw <- function(x){
  2 * IQR(x) / length(x)^(1/3)
}

time_lbl <- c(
  pre = 'Baseline',
  post = 'Follow-up'
)

# Pre/post cotinine by self-reported use
plot_cot_self_pre <-
  ggplot(df_long, aes(cotinine, fill = cotinine>=3)) + 
  geom_histogram(binwidth = 0.35) + 
  geom_vline(xintercept = 3, linetype = 'dashed', color='red') +
  scale_x_log10(breaks = c(0.1, 1, 3, 10, 100, 500)) +
  scale_y_continuous(breaks = c(2,4,6,8,10,12,14)) +
  scale_fill_binary() +
  guides(fill=guide_legend(title = "Probable recent use")) +
  facet_grid(Time~`Self-reported use`, labeller = labeller(`Self-reported use` = label_both, Time = time_lbl)) +
  labs(x = '\nCotinine ng/ml', y = 'Count\n') +
  theme_bw(15)
plot_cot_self_pre

percent_nonusers <- signif(100*sum(df$non_user)/nrow(df), 3)
percent_nonusers

percent_nonpreg_tobacco <- signif(100*sum(df$Pregnant == 'No' & df$tobacco_user2, na.rm = T)/sum(df$Pregnant == 'No', na.rm = T), 2)
percent_nonpreg_tobacco

percent_preg_tobacco <- signif(100*sum(df$Pregnant == 'Yes' & df$tobacco_user2, na.rm = T)/sum(df$Pregnant == 'Yes', na.rm = T), 2)
percent_preg_tobacco

# Misreporting tobacco use ----------------------------------------------------------

misreportPreOrPost <- signif(100*sum(df$Underreporting)/nrow(df), 2)
misreportPreOrPost2 <- signif(100*sum((df$nicotine1 == 0 & df$pre>=6.5) | (df$nicotine2 == 0 & df$post>=6.5))/nrow(df), 2)
overreportPreOrPost <- signif(100*sum(df$Overreporting)/nrow(df), 2)

# Explore predictors of tobacco use  --------------------------------------

explorevars <- 
  c('ID',
    'age',
    'education',
    'pregnancy_status',
    'monthly_income', 
    'arranged_marriage', 
    'number_children',
    'family_tobacco_use2',
    'friend_tobacco_use2',
    'mother_tobacco_use',
    'motherinlaw_tobacco_use',
    # 'tobacco_ageofonset2', # Too many missing values
    'marital_status',
    'total_baseline_harms',
    'domestic_work_hrs',
    'nondomestic_work_hrs'
  )

# Nice names
explore_rename <- c(
  "(Intercept)" = "Intercept",
  "age" = "Age",
  "pregnancy_statusNP" = "Not pregnant",
  "pregnancy_statusP" = "Pregnant",
  "pregnant" = "Pregnant",
  'family_tobacco_use2' = "Family tobacco use",
  'friend_tobacco_use2' = "Friend tobacco use",
  'mother_tobacco_use' = "Mother tobacco use",
  'motherinlaw_tobacco_use' = "Mother-in-law tobacco use",
  "peer_tobacco_useNo" = "No family tobacco use",
  "peer_tobacco_useYes" = "Family tobacco use",
  "education" = "Education (years)",
  "monthly_income" = "Monthy income (rupees)",
  "arranged_marriage" = "Arranged marriage",
  "arranged_marriageYes" = "Arranged marriage",
  "arranged_marriageNo" = "Love marriage",
  "number_children" = "Number of children",
  "married" = "Married",
  "marital_statusMarried" =  "Married",
  "marital_statusNot married" = "Not married",
  "total_baseline_harms" = "Number of freelisted harms (baseline)",
  'domestic_work_hrs' = "Domestic work hours",
  'nondomestic_work_hrs' = "Non-domestic work hours"
)

# Prepare data for elastic net regression
# Standardize continuous vars by 2SD, per Gelman
dfw0 <- 
  df %>% 
  left_join(women[explorevars]) %>%
  mutate(
    pregnant = ifelse(pregnancy_status == 'P', T, F),
    arranged_marriage = ifelse(arranged_marriage == 'Yes', T, F),
    married = ifelse(marital_status == 'Married', T, F),
    user = !non_user
  ) %>% 
  dplyr::select(
    GroupID,
    user,
    Underreporting,
    pregnant,
    married,
    all_of(explorevars),
    -ID, -pregnancy_status, -marital_status
  ) %>% 
  # Standardize numeric vars to 2SD to resemble dichotomous vars, per Gelman
  map_if(is.numeric, ~ (. - mean(.))/(2*sd(.))) %>% 
  as_tibble

dfw1 <- dplyr::select(dfw0, -Underreporting, -GroupID)

# Cross-validate glmnet many times to
# obtain stable estimate of lambda.min
cvm <- 0
for (i in 1:80){
  cv.user <- cv.glmnet(user  ~., family = 'binomial', data = dfw1, alpha = 0.2, standardize = F)
  cvm <- cvm + cv.user$cvm
}
lambda.min <- cv.user$lambda[which.min(cvm)]

coef_user <- coef(cv.user, s = lambda.min)
rownames(coef_user) <- explore_rename[rownames(coef_user)]

# Predictors of under-reporting -------------------------------------------

# dfw2 <- 
#   dfw0 %>%
#   mutate(Underreporting= factor(Underreporting)) %>% 
#   dplyr::filter(!non_user) %>% 
#   dplyr::select(-non_user)
# 
# cv.under <- cv.glmnet(Underreporting  ~., family = 'binomial', data = dfw2, alpha = 0.5, standardize = F)
# coef_under <- coef(cv.under, s = cv.under$lambda.min)
# rownames(coef_under) <- explore_rename[rownames(coef_under)]
# ggdotchart(exp(coef_under[-1, 1])) + scale_x_log10()

# Psychological outcomes analysis -----------------------------------------

followup <- 
  followup %>%
  left_join(women[c('ID', 'GroupID', 'Presentation', 'tobacco_use')]) %>% 
  left_join(df[c('ID', 'tobacco_user', 'baseline_user', 'followup_user', 'non_user')]) %>% 
  mutate(
    cuttingback2 = factor(cuttingback2, levels = c('Never', 'Sometimes')), 
    thoughtquitting_binary = ifelse(thoughtquitting == 'sometimes', 'often', as.character(thoughtquitting)),
    consequence_binary = ifelse(consequence == 'sometimes', 'often', as.character(consequence)),
    craved_binary = ifelse(craved == 'sometimes', 'often', as.character(craved))
  )
  
followup_usersonly <-
  followup %>% 
  dplyr::filter(
    !non_user
  )

# Prepare data for ggplot mosaic plot
cog_table <-
  followup_usersonly %>% 
  dplyr::select(ID, cuttingback2, thoughtquitting, consequence, craved) %>% 
  gather(key = Variable, value = Frequency, -ID) %>% 
  mutate(
    Frequency = str_to_title(Frequency),
    Frequency = factor(Frequency, levels = c('Never', 'Sometimes', 'Often')),
    Variable = factor(
      Variable, 
      levels = c('cuttingback2', 'thoughtquitting', 'consequence', 'craved'),
      labels = c('Cutting back', 'Thought of quitting', 'Consequences', 'Craved')
    )
  ) %>% 
  na.omit

p_cog <- 
  ggplot(data=cog_table) + 
  geom_mosaic(
    aes(
      x=product(Frequency, Variable), 
      fill = Frequency
    ),
    show.legend = F
  ) +
  scale_fill_viridis(discrete = T) +
  scale_x_productlist('') +
  scale_y_productlist('') +
  theme_minimal(15)


# PCA of psychological vars as requested by editor ----------------------------------------------

# Prepare psychological data for PCA
pca_cog <-
  followup_usersonly %>% 
  dplyr::select(ID, Presentation, cuttingback2, thoughtquitting, consequence, craved) %>% 
  mutate(
    across(-c(ID, Presentation), ~str_to_lower(.)),
    across(
      -c(ID, Presentation),
      ~case_when(
        . == 'never' ~ 0,
        . == 'sometimes' ~ 1,
        . == 'often' ~ 2
      )
    ),
    Presentation = ifelse(Presentation == 'GHP', 0, 1)
  ) %>% 
  na.omit

m_pca <- prcomp(pca_cog[-c(1,2)], scale. = F)
pca_cog$PC1 <- m_pca$x[,1]

ggplot(pca_cog, aes(PC1, factor(Presentation))) + geom_violin() + geom_jitter()

cor.test(~Presentation+PC1, pca_cog)

# Tobacco use frequency ---------------------------------------------------

# Tobacco use X Presentation table
tobacco_use_presentation <- xtabs(~tobacco_use+Presentation, followup)

# Prepare data for model of followup frequency
followup <-
  followup %>% 
  left_join(women[c('ID', 'tobacco_24hr_freq2', 'trimester')])

followup$Presentation <- factor(followup$Presentation)
followup$pregnancy_status <- factor(followup$trimester)

m_tobacco_freq <-
  glmmTMB(
    tobaccotoday2 ~ 
      trimester +
      tobacco_24hr_freq2 *
      Presentation +
      (1|GroupID),
    family = nbinom1,
    data = followup
  )

# Salience ----------------------------------------------------------------

df_baseline_salience0 <-
  women %>% 
  dplyr::select(ID, Presentation, starts_with('specific_harm')) %>% 
  gather(key = order, value = harm, -ID, -Presentation) %>% 
  mutate(harm = factor(harm, levels = unique(specific_harms))) %>% 
  mutate(order = as.integer(str_sub(order, 14,-1))) %>% 
  na.omit %>% 
  group_by(ID) %>% 
  mutate(salience = (max(order) - order + 1)/max(order)) %>% 
  ungroup 

# There are some dupes because some harms distinguished
# by participants we group together
#
# Remove

df_baseline_salience0 <- df_baseline_salience0[!duplicated(df_baseline_salience0[c(1,2,4)]),]

df_baseline_wide <- 
  df_baseline_salience0 %>% 
  dplyr::select(-Presentation, -order) %>% # work around tidyr bug
  spread(key = harm, value = salience, fill = 0, drop = F) %>% 
  left_join(women[c('ID', 'Presentation')]) %>%  # work around tidyr bug
  group_by(ID, Presentation) %>% # Need to keep Presentation
  summarise_all(sum) %>% 
  ungroup

df_baseline_long <-
  df_baseline_wide %>%
  gather(key = harm, value = salience, -Presentation, -ID) %>% 
  mutate(Type = harms[harm]) %>% 
  ungroup

# Simplify Figure per reviewer request
harms2 <- str_replace(harms, ' presentations', '')
harms2 <- str_replace(harms2, ' presentation', '')
harms2 <- str_replace(harms2, ' emic', '')
names(harms2) <- names(harms)

df_baseline_mean <-
  df_baseline_long %>% 
  group_by(harm) %>% 
  summarise(mean = mean(salience)) %>% 
  mutate(
    Type = factor(harms[harm]),
    Type2 = case_when(
      str_detect(Type, 'general')  ~ 'General',
      str_detect(Type, 'nothing')  ~ 'General',
      str_detect(Type, 'both')  ~ 'General',
      str_detect(Type, 'reproductive')  ~ 'Reproductive',
      TRUE ~ 'Other'
    ),
    Presentation = ifelse(str_detect(Type, 'presentation'), 'Yes', 'No'),
    # Type = factor(harms2[harm], levels = c('general', 'reproductive', 'both', 'nothing')),
    Mentioned = ifelse(mean > 0, 'Yes', 'No')
  ) %>% 
  dplyr::filter(Type2 != 'Other') %>% 
  bind_rows(tibble(harm='high blood pressure', mean=0, Type='Both', Type2='Reproductive', Presentation='Yes', Mentioned='No'))

plot_baseline_salience <-
  ggplot(
    df_baseline_mean, 
    aes(mean, forcats::fct_reorder(harm, mean), colour = Presentation)
    ) + 
  geom_point(size=3) +
  scale_colour_binary(direction = 1, guide=guide_legend(reverse=T)) +
  facet_wrap(~Type2, scales = 'free_y') +
  labs(
    title = "Salience of harms at baseline",
    x = '\nMean salience',
    y = '',
    colour = 'Included in presentation?') +
  theme_minimal(15)

plot_baseline_salience

informant_num <- length(unique(followup$ID))

df_followup_salience0 <-
  followup %>% 
  dplyr::select(ID, Presentation, starts_with('specific_harm')) %>% 
  gather(key = order, value = harm, -ID, -Presentation) %>% 
  mutate(harm = factor(harm, levels = unique(specific_harms))) %>% 
  mutate(order = as.integer(str_sub(order, 14,-1))) %>% 
  na.omit %>%
  group_by(ID) %>% 
  mutate(salience = (max(order) - order + 1)/max(order)) 

# Remove dupes here for the same reason as above.
df_followup_salience0 <- df_followup_salience0[!duplicated(df_followup_salience0[c(1,2,4)]),]

df_followup_wide <- 
  df_followup_salience0 %>% 
  dplyr::select(-Presentation, -order) %>% # work around tidyr bug
  spread(key = harm, value = salience, fill = 0, drop = F) %>% 
  left_join(women[c('ID', 'Presentation')]) %>%  # work around tidyr bug
  group_by(ID, Presentation) %>% # Need to keep Presentation
  summarise_all(sum) %>% 
  ungroup

df_followup_long <-
  df_followup_wide %>%
  gather(key = harm, value = salience, -Presentation, -ID) %>% 
  mutate(Type = harms[harm]) %>% 
  ungroup

df_followup_salience <-
  df_followup_salience0 %>% 
  group_by(Presentation, harm) %>% 
  summarise(mean = sum(salience)/informant_num)

repro <- 
  c(
    'pregnancy loss', 
    'premature delivery', 
    'low birth weight', 
    'fetal brain development', 
    'early menopause', 
    'reproductive cancer'
    )

df_followup_salience_combine <-
  df_followup_salience0 %>% 
  group_by(harm) %>% 
  summarise(mean = sum(salience)/informant_num) %>% 
  mutate(
    Type = ifelse(harm %in% repro, 'Reproductive', 'Other'),
    Type = factor(Type, levels = c('Reproductive', 'Other'))
  )

# Reproductive only

df_followup_salience_repro <-
  df_followup_salience0 %>% 
  dplyr::filter(Presentation == 'RHP') %>% 
  group_by(harm) %>% 
  summarise(mean = sum(salience)/informant_num) %>% 
  mutate(
    Type = ifelse(harm %in% repro, 'Reproductive', 'Other'),
    Type = factor(Type, levels = c('Reproductive', 'Other'))
  )

# Compare change

df_full <- 
  df_followup_long %>% 
  dplyr::select(-Presentation, -Type) %>% 
  left_join(df_baseline_long, by = c('ID', 'harm')) %>% 
  rename(salience1 = salience.y, salience2 = salience.x) %>% 
  mutate(diff = salience2 - salience1)

# Define function
mean_salience <-
  . %>% 
  group_by(Presentation, harm) %>% 
  summarise(s1mean = mean(salience1), s2mean = mean(salience2), meandiff = mean(diff))

df_full_mean <-
  df_full %>% mean_salience

mean_diff <- function(d, w){
  df <- mean_salience(d[w,])
  df$meandiff
}

b <- boot(df_full, mean_diff, 1000, strata = factor(df_full$harm))
b_ci <- confint(b, type = 'perc')
df_full_mean$lower <- b_ci[,1]
df_full_mean$upper <- b_ci[,2]

dumbbell_plot <- function(pres) {
    
    type <- ifelse(pres == 'GHP', 'General', 'Reproductive')
    cols <- magma(11)[c(5,8)]
    names(cols) <- c('Decrease', 'Increase')
    
    df <-
      df_full_mean %>% 
      dplyr::filter(Presentation == pres) %>% 
      mutate(
        harm2 = factor(harm, levels = harm[order(ifelse(s2mean > 0, s2mean, -1-meandiff))]),
        Change = ifelse(meandiff < 0, 'Decrease', 'Increase'), 
        Sig = ifelse(sign(lower) == sign(upper) & lower != 0, 1, 0.3)
      ) %>% 
      ggplot(
        aes(x = s1mean, xend = s2mean, y = harm2, colour = Change, alpha = Sig)
      ) +
      geom_dumbbell(
        colour_x = 'black', 
        size = 3
        ) +
      geom_point(size = 3, alpha=1, colour='black') +
      scale_alpha(range = c(0.3, 1)) +
      scale_colour_binary(guide=guide_legend(reverse=T)) +
      scale_x_continuous(limits = c(0, 0.6)) +
      labs(
        title = {type}
        ) +
      guides(alpha = F) +
      theme_minimal(15)
  }

df2 <-
  followup %>%
  dplyr::select(
    ID, 
    quit2, 
    thoughtquitting2, 
    tobaccotoday2, 
    cuttingback2,
    consequence,
    craved,
    presentation_share2,
    number_general,
    number_reproductive,
    number_emic,
    delta_days) %>% 
  mutate(ID = str_replace(ID, 'NP-', '')) %>% 
  mutate(ID = str_replace(ID, 'P-', '')) %>% 
  gather(key = harm_type, value = followup_number, number_general, number_reproductive, number_emic) %>% 
  mutate(harm_type = str_replace(harm_type, 'number_', ''))

df2 <-
  women %>%
  dplyr::select(
    ID,
    GroupID,
    Presentation, 
    pregnancy_status,
    trimester,
    number_general,
    number_reproductive,
    number_emic
    ) %>%
  gather(key = harm_type, value = baseline_number, number_general, number_reproductive, number_emic) %>% 
  mutate(harm_type = str_replace(harm_type, 'number_', '')) %>% 
  left_join(df2, by = c('ID', 'harm_type')) %>% 
  mutate(
    harm_type = factor(harm_type), 
    Presentation = factor(Presentation, levels = c('GHP', 'RHP')),
    pregnancy_status = factor(pregnancy_status)
    )

df2 <-
  df2 %>% left_join(df[c('ID', 'tobacco_user')]) %>% 
  mutate(tobacco_user = factor(tobacco_user))

# total harms
df_fharms <-
  df2 %>% 
  group_by(ID) %>% 
  summarise(total_followup_harms = sum(followup_number), total_baseline_harms = sum(baseline_number))

m_wilcox <- wilcox.test(df_fharms$total_baseline_harms, df_fharms$total_followup_harms, paired=T)

# Regression models of harm types -----------------------------------------

# remove emic harms because no predictions about them

df2b <-
  df2 %>% 
  dplyr::filter(harm_type != 'emic', !is.na(followup_number)) %>% 
  mutate(
    harm_type = case_when(
      harm_type == 'general' ~ 'General harms',
      harm_type == 'reproductive' ~ 'Reproductive harms'
    )
  )

m_harms <- glmmTMB(
  followup_number ~ 
    baseline_number +
    trimester +
    Presentation * harm_type +
    (1|GroupID/ID),
  family = poisson,
  data = df2b
)

mean_baseline_general <- mean(df2b$baseline_number[df2b$harm_type == 'General harms'])
mean_baseline_repro <- mean(df2b$baseline_number[df2b$harm_type == 'Reproductive harms'])

df_annotate <- tibble(
  Presentation = c('RHP', 'RHP'), 
  label = c('General harms baseline', 'Reproductive harms baseline'),
  x = c(1.1,1.1),
  y = c(mean_baseline_general, mean_baseline_repro),
  harm_type = c('')
  )

p_m_harms <- 
  visreg(m_harms, by ="Presentation", xvar = 'harm_type', scale = 'response', gg = 'T', rug=F) +
  geom_hline(yintercept = mean_baseline_general, linetype = 'dotted') +
  geom_hline(yintercept = mean_baseline_repro, linetype = 'dashed') +
  geom_text(data = df_annotate, aes(x, y, label = label), hjust=0) +
  scale_y_continuous(limits = c(0, 4)) +
  coord_cartesian(clip='off', xlim=c(0,1)) +
  labs(x = '\nHarm type', y = 'Number of harms mentioned at follow-up\n') +
  theme_bw(15) +
  theme(
    plot.margin = unit(c(1,11,1,1), "lines"), 
    panel.background = element_rect(fill='white', colour='gray'),
    panel.grid.major = element_line(colour = '#eeeeee')
    )
p_m_harms

# Cotinine regression models ----------------------------------------------

#' Effect of Presentation type (G = General harm; R = Reproductive harm) on
#' post-intervention cotinine, controlling for pre-intervention cotinine and
#' trimester. A significant main effect for Presentation, or significant interaction between
#' Presentation and pre-intervention smoking indicates the intervention had an effect (the latter on
#' the argument that the intervention should have a bigger effect on heavier tobacco_users).

#' All data, Gaussian errors

m_post1 <- lmer(
  post ~ 
    pre + 
    Presentation + 
    trimester +
    (1|GroupID),
  data = df
)

m_post2 <- lmer(
  post ~ 
    pre * Presentation + 
    trimester +
    (1|GroupID),
  data = df
  )

# The predict method of this model produces SE's for visreg
m_post2b <- glmmTMB(
  post ~ 
    pre * Presentation + 
    trimester +
    (1|GroupID),
  family = gaussian,
  data = df
)

p_m_post <- 
  visreg(m_post2b, xvar = 'pre', by = 'Presentation', gg = T) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dotted') +
  labs(x = '\nBaseline salivary cotinine (ng/ml)', y = 'Followup salivary cotinine (ng/ml)\n') +
  theme_bw(15)
p_m_post

# Sharing presentation ----------------------------------------------------

df3 <-
  df2 %>% 
  group_by(ID, Presentation, presentation_share2, GroupID) %>% 
  summarise(baseline_number = sum(baseline_number), followup_number = sum(followup_number)) %>% 
  mutate(Share = presentation_share2 == 'yes') %>% 
  left_join(women[c('ID', 'trimester')])

m_share <-
  glmer(
    Share ~
      followup_number +
      Presentation +
      trimester +
      (1|GroupID),
      family = binomial,
    data = df3
    )
  
summary(m_share)


# Effect of delta harms on cotinine ---------------------------------------

df4 <- 
  left_join(df2b[c('ID', 'GroupID', 'trimester', 'harm_type', 'baseline_number', 'followup_number')], df[c('ID', 'Presentation', 'pre', 'post')]) %>% 
  mutate(delta_number = followup_number - baseline_number) %>%
  dplyr::select(-followup_number, -baseline_number) %>% 
  pivot_wider(names_from = 'harm_type', values_from = delta_number) %>% 
  rename(General_harms = `General harms`, Reproductive_harms = `Reproductive harms`) %>% 
  mutate(Total_harms = General_harms + Reproductive_harms)

m_cot_deltargharms <- glmmTMB(log10(post) ~ log10(pre) + trimester + Presentation + General_harms + Reproductive_harms + (1|GroupID), family = gaussian, data = df4)
summary(m_cot_deltargharms)

p_m_cot_deltargharmsGH <- 
  visreg(m_cot_deltargharms, xvar='General_harms', gg=T) +
  labs(x='\nΔNumber of general harms', y='Followup log10 salivary cotinine\n') +
  theme_bw()
p_m_cot_deltargharmsGH

p_m_cot_deltargharmsRH <- 
  visreg(m_cot_deltargharms, xvar='Reproductive_harms', gg=T) +
  labs(x='\nΔNumber of reproductive harms', y='') +
  theme_bw()
p_m_cot_deltargharmsRH


# Alternative regression table --------------------------------------------

var_dict <- c(
  '(Intercept)' = '(Intercept)',
  # 'sd__(Intercept)' = 'SD (random intercept)',
  # 'sd__Observation' = 'SD (observation)',
  'trimester2' = 'Trimester 2',
  'trimester3' = 'Trimester 3',
  'trimesterNot pregnant' = 'Not pregnant',
  'PresentationRHP' = 'RHP',
  'tobacco_24hr_freq2' = 'Baseline Tobacco Frequency',
  'tobacco_24hr_freq2:PresentationRHP' = 'Baseline Tobacco Frequency * RHP',
  'pre' = 'Baseline cotinine',
  'pre:PresentationRHP' = 'Baseline cotinine * RHP',
  'baseline_number' = 'Baseline number of harms',
  # 'Reproductive_harms' = 'Number of reproductive harms',
  # 'General_harms' = 'Number of general harms',
  'harm_typeReproductive harms' = 'Number of reproductive harms',
  'PresentationRHP:harm_typeReproductive harms' = 'Number of reproductive harms * RHP',
  'followup_number' = 'Total number of harms at followup'
)

varnames <- c('term', 'estimate', 'std.error', 'statistic', 'p.value', 'conf.low', 'conf.high')
# varnames_round <- c('estimate', 'std.error', 'statistic', 'conf.low', 'conf.high')

models <- list(
  'Tobacco self-report' = m_tobacco_freq, 
  'Cotinine 1' = m_post1, 
  'Cotinine 2' = m_post2b,
  'Number of harms' = m_harms, 
  'Share presentation' = m_share
  )

model_stats <- 
  map_df(models, ~tidy(., conf.int=T)) %>% 
  dplyr::select(all_of(varnames)) %>%
  dplyr::filter(!str_detect(term, 'sd')) %>% 
  mutate(
    term = var_dict[term],
    across(estimate:statistic, ~round(., 2)),
    p.value = signif(p.value, 2),
    across(conf.low:conf.high, ~signif(., 3)),
    across(everything(), as.character),
    across(everything(), ~ ifelse(is.na(.), '', .))
    ) %>% 
  rename(
    Variable = term,
    Estimate = estimate,
    SE = std.error,
    'Z-value' = statistic,
    'P-value' = p.value,
    'Lower 95% CI' = conf.low,
    'Upper 95% CI' = conf.high
  )

