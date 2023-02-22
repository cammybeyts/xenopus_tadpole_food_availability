rm(list = ls())
if(Sys.info()[["user"]] == "s1437006"){
  setwd("/Users/s1437006/Documents/Edinburgh_PhD_documents/Papers/Xenopus_2020/Github/")
} 
if(Sys.info()[["user"]] == "cbeyts"){
  setwd("/Users/cbeyts/Documents/Edinburgh_PhD_documents/Papers/Xenopus_2020/Github/")
}else{
  setwd("/Users/cammybeyts/Documents/Edinburgh_PhD_documents/Papers/Xenopus_2020/Github/")
}



#### libraries ####
library(rstan)
library(brms)
library(loo)
library(bayesplot)
library(tidyverse)
library(kableExtra)
library(tidybayes)

#### loading models and data ####
load("tadpole_noegg_model_yrep.rda")

read.data <- read.delim("Xenopus_ACT_EXP.txt", header = TRUE)

#some data formating
full.data1 <- read.data %>%
  mutate(NameTreatment = Treatment) %>%
  mutate(Treatment = recode(Treatment, Low = 1, High = 2)) %>%
  # low high converted into 1, 2
  mutate(Treatment = as.factor(Treatment)) %>%
  mutate(exp = as.numeric(scale(EXP_Total_Dist))) %>%
  mutate(act = as.numeric(scale(ACT_Total_Dist))) %>%
  mutate(Scale_SVL = scale(SVL)) %>%
  mutate(TadpoleID = as.factor(TadpoleID)) %>%
  mutate(Set = as.factor(Set)) %>%
  mutate(Egg_MassID = as.factor(Egg_MassID)) %>%
  mutate(Rep = recode(Rep, 0, 1, 2, 3, 4, 5, 6, 7)) %>%
  # reps 12345678 recoded to 01234567 %>%
  mutate(RepCon = as.numeric(Rep))

#### param list ####

##list of parameter names from stan model 
#helps to select parameters for plotting ect
params <- c(
  "b_ScaleACTTotalDist_Intercept",
  "b_ScaleACTTotalDist",
  "b_ScaleEXPTotalDist10_Intercept",
  "b_ScaleEXPTotalDist10",
  "b_sigma_ScaleACTTotalDist_Intercept",
  "b_sigma_ScaleACTTotalDist",
  "b_sigma_ScaleEXPTotalDist10_Intercept",
  "b_sigma_ScaleEXPTotalDist10",
  "V_tadpole_Low",
  "V_tadpole_High",
  "cor_tadpole_Low",
  "cor_tadpole_High",
  "rescor"
)

##rename parameter names to make easier to work with data
params_names <- c(
  "M_act_m_int",
  "M_act_m_treatCont", "act_m_SVL", "act_m_rep", "act_m_set",
  "M_exp_m_int",
  "M_exp_m_treatCont", "exp_m_SVL", "exp_m_rep", "exp_m_set",
  "S_act_v_int",
  "S_act_v_treatCont",
  "S_exp_v_int",
  "S_exp_v_treatCont",
  "V_id_em_Low",
  "V_id_es_Low",
  "V_id_ev_Low",
  "V_id_am_Low",
  "V_id_as_Low",
  "V_id_av_Low",
  "V_id_em_High",
  "V_id_es_High",
  "V_id_ev_High",
  "V_id_am_High",
  "V_id_as_High",
  "V_id_av_High",
  paste("cor_id", c("em"), rep("es", 1), "Low", sep = "_"),
  paste("cor_id", c("em", "es"), rep("ev", 2), "Low", sep = "_"),
  paste("cor_id", c("em", "es", "ev"), rep("am", 3), "Low", sep = "_"),
  paste("cor_id", c("em", "es", "ev", "am"), rep("as", 4), "Low", sep = "_"),
  paste("cor_id", c("em", "es", "ev", "am", "as"), rep("av", 5), "Low",
        sep = "_"
  ),
  paste("cor_id", c("em"), rep("es", 1), "High", sep = "_"),
  paste("cor_id", c("em", "es"), rep("ev", 2), "High", sep = "_"),
  paste("cor_id", c("em", "es", "ev"), rep("am", 3), "High", sep = "_"),
  paste("cor_id", c("em", "es", "ev", "am"), rep("as", 4), "High", sep = "_"),
  paste("cor_id", c("em", "es", "ev", "am", "as"), rep("av", 5), "High",
        sep = "_"
  ),
  "rescor_Low",
  "rescor_High"
)

params_names_full <- c(
  "M_act_m_int",
  "M_act_m_treatCont", "act_m_SVL", "act_m_rep", "act_m_set",
  "M_exp_m_int",
  "M_exp_m_treatCont", "exp_m_SVL", "exp_m_rep", "exp_m_set",
  "S_act_v_int",
  "S_act_v_treatCont",
  "S_exp_v_int",
  "S_exp_v_treatCont",
  "V_id_em_Low",
  "V_id_es_Low",
  "V_id_ev_Low",
  "V_id_am_Low",
  "V_id_as_Low",
  "V_id_av_Low",
  "V_id_em_High",
  "V_id_es_High",
  "V_id_ev_High",
  "V_id_am_High",
  "V_id_as_High",
  "V_id_av_High",
  paste(c("ExpVI(m)"), rep("ExpVIxE(m)", 1), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)"), rep("ExpVI(v)", 2), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)"), rep("ActVI(m)", 3), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)", "ActVI(m)"), rep("ActVIxE(m)", 4), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)", "ActVI(m)", "ActVIxE(m)"), rep("ActVI(v)", 5),
        sep = "_"
  ),
  paste(c("ExpVI(m)"), rep("ExpVIxE(m)", 1), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)"), rep("ExpVI(v)", 2), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)"), rep("ActVI(m)", 3), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)", "ActVI(m)"), rep("ActVIxE(m)", 4), sep = "_"),
  paste(c("ExpVI(m)", "ExpVIxE(m)", "ExpVI(v)", "ActVI(m)", "ActVIxE(m)"), rep("ActVI(v)", 5),
        sep = "_"
  ),
  "Exp Act residual correlation",
  "Exp Act residual correlation"
)
#### ggplot theme for plots ####

##some ggplot themes for easy plotting
theme_catplots <- function() {
  theme_classic() + # use a white background
    theme(
      axis.title.x = element_text(face = "bold", size = 14),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 14),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      strip.text.x = element_text(face = "bold", size=12),
      strip.background = element_rect(colour=FALSE, fill=FALSE),
      legend.position="none"
    )
}

theme_boxplots <- function(){
  theme_classic() + # use a white background
    theme(axis.title.x = element_text(face="bold", size=12),
          axis.title.y = element_text(face="bold", size=12),
          axis.text.y = element_text(face="bold", size=12),
          axis.text.x = element_text(face="bold", size=12),
          plot.title = element_text(hjust = 0.5, face = "bold", size=14),
          legend.position="none")
}
#### functions to calculate post difference ####

#function to calculate the difference between the high vs the low feeding treatments
post_diff <- function(data_pd, params, type) {
  if (type == "M") { #mean model
    vcHigh <- data_pd[, paste("M", params[1], "treatHigh", sep = "_")]
    vcLow <- data_pd[, paste("M", params[1], "treatLow", sep = "_")]
    vcDiff <- data.frame(vcHigh - vcLow)
  } 
  if (type == "S") { #resiuals
    vcHigh <- data_pd[, paste("S", params[1], "treatHigh", sep = "_")]
    vcLow <- data_pd[, paste("S", params[1], "treatLow", sep = "_")]
    vcDiff <- data.frame(vcHigh - vcLow)
  } 
  if (type == "V") { #variance among individuals
    vcHigh <- data_pd[, paste("V_id", params[1], "High", sep = "_")]
    vcLow <- data_pd[, paste("V_id", params[1], "Low", sep = "_")]
    vcDiff <- data.frame(vcHigh - vcLow)
  } 
  if (type == "rescor") { #residual variance
    vcHigh <- data_pd[, paste(params[1], "High", sep = "_")]
    vcLow <- data_pd[, paste(params[1], "Low", sep = "_")]
    vcDiff <- data.frame(vcHigh - vcLow)
  }
  else if (type == "corr") { #correlations among individuals
    vcLow <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "Low",
      sep = "_"
    )]
    vcHigh <- data_pd[, paste(
      "cor_id", paste(params, collapse = "_"), "High",
      sep = "_"
    )]
    vcDiff <- data.frame(vcHigh - vcLow)
  }
  return(vcDiff)
}

#### Get parameters from model and rename #####
noegg_pd <- as.data.frame(noegg_model, pars = params)
names(noegg_pd) <- params_names
#### remove treatment differences from fixed effects#####
#removes the intercept from the mean and dispersion model, allowing us to see estimate of ACT and EXP parameters, not comparisons
noegg_pd <- noegg_pd %>%
  mutate(
    M_act_m_treatLow = M_act_m_int,
    M_act_m_treatHigh = M_act_m_int + M_act_m_treatCont,
    S_act_v_treatLow = S_act_v_int,
    S_act_v_treatHigh = S_act_v_int + S_act_v_treatCont,
    M_exp_m_treatLow = M_exp_m_int,
    M_exp_m_treatHigh = M_exp_m_int + M_exp_m_treatCont,
    S_exp_v_treatLow = S_exp_v_int,
    S_exp_v_treatHigh = S_exp_v_int + S_exp_v_treatCont
  )

#### treatment differences #####
## data for differences used for creating summary plots

#calculate treatment differences between the high vs low treatments
noegg_pd <- noegg_pd %>%
  mutate(
    d_M_act = unlist(post_diff(., params = "act_m", type = "M")),
    d_M_exp = unlist(post_diff(., params = "exp_m", type = "M")),
    d_S_act = unlist(post_diff(., params = "act_v", type = "S")),
    d_S_exp = unlist(post_diff(., params = "exp_v", type = "S")),
    d_V_id_em = unlist(post_diff(., params = "em", type = "V")),
    d_V_id_am = unlist(post_diff(., params = "am", type = "V")),
    d_V_id_es = unlist(post_diff(., params = "es", type = "V")),
    d_V_id_as = unlist(post_diff(., params = "as", type = "V")),
    d_V_id_ev = unlist(post_diff(., params = "ev", type = "V")),
    d_V_id_av = unlist(post_diff(., params = "av", type = "V")),
    d_cor_id_em_es = unlist(post_diff(., params = c("em", "es"), type = "corr")), 
    d_cor_id_em_ev = unlist(post_diff(., params = c("em", "ev"), type = "corr")),
    d_cor_id_es_ev = unlist(post_diff(., params = c("es", "ev"), type = "corr")),
    d_cor_id_em_am = unlist(post_diff(., params = c("em", "am"), type = "corr")),
    d_cor_id_es_am = unlist(post_diff(., params = c("es", "am"), type = "corr")),
    d_cor_id_ev_am = unlist(post_diff(., params = c("ev", "am"), type = "corr")),
    d_cor_id_em_as = unlist(post_diff(., params = c("em", "as"), type = "corr")),
    d_cor_id_es_as = unlist(post_diff(., params = c("es", "as"), type = "corr")),
    d_cor_id_ev_as = unlist(post_diff(., params = c("ev", "as"), type = "corr")),
    d_cor_id_am_as = unlist(post_diff(., params = c("am", "as"), type = "corr")),
    d_cor_id_em_av = unlist(post_diff(., params = c("em", "av"), type = "corr")),
    d_cor_id_es_av = unlist(post_diff(., params = c("es", "av"), type = "corr")),
    d_cor_id_ev_av = unlist(post_diff(., params = c("ev", "av"), type = "corr")),
    d_cor_id_am_av = unlist(post_diff(., params = c("am", "av"), type = "corr")),
    d_cor_id_as_av = unlist(post_diff(., params = c("as", "av"), type = "corr")), 
    d_rescor = unlist(post_diff(., params = "rescor", type = "rescor"))
  )

#####preparing data for summary tables#####


Variable_names_table <- c("SVL",
                          "Trial",
                          "Set",
                          "Low",
                          "High",
                          "High - Low",
                          "Low ",
                          "High ",
                          "High - Low ",
                          "Low  ",
                          "High  ",
                          "High - Low  ",
                          "Low   ",
                          "High   ",
                          "High - Low   ",
                          "Low     ",
                          "High     ",
                          "High - Low     "
                          )

noegg_pd_mean_ACT <- noegg_pd %>%
  select(act_m_SVL,
         act_m_rep,
         act_m_set,
         M_act_m_treatLow,
         M_act_m_treatHigh,
         d_M_act,
         S_act_v_treatLow,
         S_act_v_treatHigh,
         d_S_act,
         V_id_am_Low,
         V_id_am_High,
         d_V_id_am,
         V_id_as_Low,
         V_id_as_High,
         d_V_id_as,
         V_id_av_Low,
         V_id_av_High,
         d_V_id_av
         ) %>%
  summarise_all(list(mean)) %>%
  round(.,3) %>%
  gather()


noegg_pd_CI_ACT <- noegg_pd %>%
  select(act_m_SVL,
         act_m_rep,
         act_m_set,
         M_act_m_treatLow,
         M_act_m_treatHigh,
         d_M_act,
         S_act_v_treatLow,
         S_act_v_treatHigh,
         d_S_act,
         V_id_am_Low,
         V_id_am_High,
         d_V_id_am,
         V_id_as_Low,
         V_id_as_High,
         d_V_id_as,
         V_id_av_Low,
         V_id_av_High,
         d_V_id_av
  ) %>%
  as.mcmc(.) %>%
  HPDinterval(.) %>%
  round(.,3) %>%
  as.data.frame()
colnames(noegg_pd_CI_ACT) <- c("actlower", "actupper")
rownames(noegg_pd_CI_ACT) <- NULL

noegg_pd_mean_EXP <- noegg_pd %>%
  select(exp_m_SVL,
         exp_m_rep,
         exp_m_set,
         M_exp_m_treatLow,
         M_exp_m_treatHigh,
         d_M_exp,
         S_exp_v_treatLow,
         S_exp_v_treatHigh,
         d_S_exp,
         V_id_em_Low,
         V_id_em_High,
         d_V_id_em,
         V_id_es_Low,
         V_id_es_High,
         d_V_id_es,
         V_id_ev_Low,
         V_id_ev_High,
         d_V_id_ev,
  ) %>%
  summarise_all(list(mean)
  ) %>%
  round(.,3) %>%
  gather(.)

noegg_pd_CI_EXP <- noegg_pd %>%
  select(exp_m_SVL,
         exp_m_rep,
         exp_m_set,
         M_exp_m_treatLow,
         M_exp_m_treatHigh,
         d_M_exp,
         S_exp_v_treatLow,
         S_exp_v_treatHigh,
         d_S_exp,
         V_id_em_Low,
         V_id_em_High,
         d_V_id_em,
         V_id_es_Low,
         V_id_es_High,
         d_V_id_es,
         V_id_ev_Low,
         V_id_ev_High,
         d_V_id_ev,
  ) %>%
  as.mcmc(.) %>%
  HPDinterval(.) %>%
  round(.,3) %>%
  as.data.frame() 
colnames(noegg_pd_CI_EXP) <- c("explower", "exupper")
rownames(noegg_pd_CI_EXP) <- NULL


noegg_pd_mean_ACT_EXP <- cbind(noegg_pd_mean_ACT[2], 
                               noegg_pd_CI_ACT, 
                               noegg_pd_mean_EXP[2], 
                               noegg_pd_CI_EXP)
colnames(noegg_pd_mean_ACT_EXP) <- c("Mean", "Lower", "Upper", "Mean", "Lower", "Upper")
rownames(noegg_pd_mean_ACT_EXP) <- Variable_names_table


#### summary tables ####
ACT_EXP_results_table <- noegg_pd_mean_ACT_EXP  %>%
  knitr::kable(
    format = "html",
    escape = F
  ) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"), full_width = T) %>% 
  add_header_above(c(" " = 2, "95% CI" = 2, " " = 1, "95% CI" = 2)) %>%
  add_header_above(c(" " = 1, "Activity" = 3, "Exploration" = 3)) %>%
  pack_rows("Fixed Effects", 1, 3) %>%
  pack_rows("Population: Mean behavior", 4, 6) %>%
  pack_rows("Population: Predictability", 7, 9) %>%
  pack_rows("Variance among individuals: Personality", 10, 12) %>%
  pack_rows("Variance among individuals: Plasticity", 13, 15) %>%
  pack_rows("Variance among individuals: Predictability", 16, 18)
ACT_EXP_results_table
#save_kable(ACT_EXP_results_table, "ACT_EXP_results_table.png")

####new plots ####
pop_act_lowCI <- noegg_pd %>%
  select(M_act_m_treatLow) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_act_lowmean <- noegg_pd %>%
  select(M_act_m_treatLow) %>%
  as.data.frame() %>%
  summarise(
    mean_param = mean(M_act_m_treatLow)
  )

pop_act_highCI <- noegg_pd %>%
  select(M_act_m_treatHigh) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_act_highmean <- noegg_pd %>%
  select(M_act_m_treatHigh) %>%
  as.data.frame() %>%
  summarise(
    M_act_m_treatHigh = mean(M_act_m_treatHigh)
  )

pop_act_low <- cbind(pop_act_lowmean, pop_act_lowCI) %>%
  mutate(Treatment = "Low")
colnames(pop_act_low) <- c("mean", "lower", "upper", "Treatment")
pop_act_high <- cbind(pop_act_highmean, pop_act_highCI) %>%
  mutate(Treatment = "High")
colnames(pop_act_high) <- c("mean", "lower", "upper", "Treatment")
pop_act <- rbind(pop_act_low , pop_act_high)

pop_act_plot <- ggplot() + 
  geom_point(data = pop_act, aes(x = Treatment, y =mean)) +
  geom_linerange(data = pop_act, aes(x = Treatment, ymin = lower, ymax = upper)) +
  xlab("\nTreatment") +
  #ylab("Scaled Distance Travelled (pixels)\n") +
  theme_boxplots()
pop_act_plot
ggsave(filename = "stan_plots/Act_mean_treat_effect.png", pop_act_plot)

#act - pop predict
pop_V_act_lowCI <- noegg_pd %>%
  select(S_act_v_treatLow) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_V_act_lowmean <- noegg_pd %>%
  select(S_act_v_treatLow) %>%
  as.data.frame() %>%
  summarise(
    mean_param = mean(S_act_v_treatLow)
  )

pop_V_act_highCI <- noegg_pd %>%
  select(S_act_v_treatHigh) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_V_act_highmean <- noegg_pd %>%
  select(S_act_v_treatHigh) %>%
  as.data.frame() %>%
  summarise(
    S_act_v_treatHigh = mean(S_act_v_treatHigh)
  )

pop_V_act_low <- cbind(pop_V_act_lowmean, pop_V_act_lowCI) %>%
  mutate(Treatment = "Low")
colnames(pop_V_act_low) <- c("mean", "lower", "upper", "Treatment")
pop_V_act_high <- cbind(pop_V_act_highmean, pop_V_act_highCI) %>%
  mutate(Treatment = "High")
colnames(pop_V_act_high) <- c("mean", "lower", "upper", "Treatment")
pop_V_act <- rbind(pop_V_act_low , pop_V_act_high)

pop_V_act_plot <- ggplot() + 
  geom_point(data = pop_V_act, aes(x = Treatment, y =mean)) +
  geom_linerange(data = pop_V_act, aes(x = Treatment, ymin = lower, ymax = upper)) +
  #xlab("\nTreatment") +
  #ylab("Scaled Distance Travelled (pixels)\n") +
  theme_boxplots()
pop_V_act_plot
ggsave(filename = "stan_plots/Act_mean_predict_treat_effect.png", pop_V_act_plot)


#mean exploration
pop_exp_lowCI <- noegg_pd %>%
  select(M_exp_m_treatLow) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_exp_lowmean <- noegg_pd %>%
  select(M_exp_m_treatLow) %>%
  as.data.frame() %>%
  summarise(
    mean_param = mean(M_exp_m_treatLow)
  )

pop_exp_highCI <- noegg_pd %>%
  select(M_exp_m_treatHigh) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_exp_highmean <- noegg_pd %>%
  select(M_exp_m_treatHigh) %>%
  as.data.frame() %>%
  summarise(
    M_exp_m_treatHigh = mean(M_exp_m_treatHigh)
  )

pop_exp_low <- cbind(pop_exp_lowmean, pop_exp_lowCI) %>%
  mutate(Treatment = "Low")
colnames(pop_exp_low) <- c("mean", "lower", "upper", "Treatment")
pop_exp_high <- cbind(pop_exp_highmean, pop_exp_highCI) %>%
  mutate(Treatment = "High")
colnames(pop_exp_high) <- c("mean", "lower", "upper", "Treatment")
pop_exp <- rbind(pop_exp_low , pop_exp_high)

pop_exp_plot <- ggplot() + 
  geom_point(data = pop_exp, aes(x = Treatment, y =mean)) +
  geom_linerange(data = pop_exp, aes(x = Treatment, ymin = lower, ymax = upper)) +
  xlab("\nTreatment") +
  ylab("Scaled Distance Travelled (pixels)\n") +
  theme_boxplots()
pop_exp_plot
ggsave(filename = "stan_plots/Exp_mean_treat_effect.png", pop_exp_plot)

#exp - pop predict
pop_V_exp_lowCI <- noegg_pd %>%
  select(S_exp_v_treatLow) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_V_exp_lowmean <- noegg_pd %>%
  select(S_exp_v_treatLow) %>%
  as.data.frame() %>%
  summarise(
    mean_param = mean(S_exp_v_treatLow)
  )

pop_V_exp_highCI <- noegg_pd %>%
  select(S_exp_v_treatHigh) %>%
  as.mcmc() %>%
  HPDinterval() %>%
  as.data.frame()

pop_V_exp_highmean <- noegg_pd %>%
  select(S_exp_v_treatHigh) %>%
  as.data.frame() %>%
  summarise(
    S_exp_v_treatHigh = mean(S_exp_v_treatHigh)
  )

pop_V_exp_low <- cbind(pop_V_exp_lowmean, pop_V_exp_lowCI) %>%
  mutate(Treatment = "Low")
colnames(pop_V_exp_low) <- c("mean", "lower", "upper", "Treatment")
pop_V_exp_high <- cbind(pop_V_exp_highmean, pop_V_exp_highCI) %>%
  mutate(Treatment = "High")
colnames(pop_V_exp_high) <- c("mean", "lower", "upper", "Treatment")
pop_V_exp <- rbind(pop_V_exp_low , pop_V_exp_high)

pop_V_exp_plot <- ggplot() + 
  geom_point(data = pop_V_exp, aes(x = Treatment, y =mean)) +
  geom_linerange(data = pop_V_exp, aes(x = Treatment, ymin = lower, ymax = upper)) +
  xlab("\nTreatment") +
  ylab("Scaled Distance Travelled (pixels)\n") +
  theme_boxplots()
pop_V_exp_plot
ggsave(filename = "stan_plots/Exp_mean_predict_treat_effect.png", pop_V_exp_plot)


