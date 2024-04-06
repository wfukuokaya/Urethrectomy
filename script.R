pacman::p_load(missForest, gdata, survminer, gtsummary, WeightIt, cobalt, smd, survival, ggsurvfit, ggsci, patchwork, ggeffects)
library(tidyverse)
set.cobalt.options(binary = "std")

d # the study data

dfull <- d %>% 
  select(fu, death_all, death_cancer, nutrfs, nutr, comp_g3, ebl,
         urethrectomy, urethrectomy_bin, age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, neck_tur, prostate_tur, cis_tur, hist_rc, ct, cn,
         type, lnd, div, pt, pn, grade, nac, ac, hgb, lognlr, risk_urethra, id) %>% mutate_if(is.character, as.factor) %>% as.data.frame() # full dataset

# define ggplot2 theme
theme_juog <- function(base_size = 10, 
                       dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.4)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

# trend analysis
dt <- d %>% filter(quarter >= 2013 & quarter < 2020) %>% 
  group_by(quarter, urethrectomy) %>% tally() %>% 
  mutate(percentage = n / sum(n)) %>% 
  filter(urethrectomy == "yes") 

trend <- dt %>% 
  ggplot2::ggplot(ggplot2::aes(x = quarter, y = percentage)) + 
  ggplot2::geom_line(linewidth = 1, color = "grey") + 
  geom_point(size = 3, color = "grey") + 
  geom_smooth(method = "loess", level = 0.95, color = "#39568CFF", fill = "#39568CFF") + 
  theme_juog() + 
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) + 
  zoo::scale_x_yearqtr(format = "%Y", n = 10) + 
  labs(x = "", 
       y = "Proportion of patients\nreceiving prophylactic urethrectomy"
  )

# imputation.based on random forest...
set.seed(3)
df_imp <- missForest(dfull, maxiter = 1, verbose = TRUE) # random forest imputation
dimp <- df_imp$ximp # extract the imputed set

# define formula for PS calculation
iptw_formula <- as.formula(urethrectomy_bin ~ 
                             age_rc + sex + bmi + ecog_ps_cat + smoking + utuc + nmibc + bcg + multi_tur + neck_tur + prostate_tur + cis_tur + hist_rc + 
                             ct + cn + type + lnd + div + pt + pn + grade + nac + ac + hgb + lognlr + risk_urethra)

# IPTW calculation using multivariable logistic regression model
wgt_out <- weightit(
  iptw_formula,
  data = dimp, stabilize = TRUE, estimand = "ATT", # evaluating average treatment effect in the treated
  method = "ps", verbose = TRUE
)

balance_imp <- bal.tab(wgt_out, abs = FALSE, un = TRUE, thresholds = c(m = 0.1)) # the summary of min/max SMDs...
balance_imp # check balance between groups

label <- data.frame(
  old = c("age_rc", "sex_male", "bmi", "ecog_ps_cat_2 or more", 
          "smoking_never", "smoking_former", "smoking_current", 
          "utuc_yes", "nmibc_yes", "bcg_yes", "multi_tur_multifocal", "neck_tur_yes", "prostate_tur_yes", 
          "cis_tur_yes", "hist_rc_Urothelial carcinoma", "hist_rc_UC with variant histology", "hist_rc_Pure non-UC", 
          "ct_cTa or cT1", "ct_cT2", "ct_cT3 or cT4", 
          "cn_cN1", 
          "type_Open", "type_Laparoscopic", "type_Robotic", 
          "lnd_yes", "div_Other", 
          "pt_pTa or pT1", "pt_pT2", "pt_pT3 or pT4", 
          "pn_pN1 to pN3", 
          "grade_High grade", "nac_yes", "ac_yes", "hgb", "lognlr", "risk_urethra_High risk"),
  new = c("Age", "Sex: Male", "BMI", "ECOG PS: 2 or more", 
          "Smoking status: Never", "Smoking status: Former", "Smoking status: Current",
          "Previous UTUC: Yes", "Previous NMIBC: Yes", "Previous BCG: Yes",
          "Tumor multifocality: Multifocal", "Bladder neck tumor: Yes", 
          "Prostatic urethral tumor involvement: Yes", "Concomitant carcinoma in situ: Yes",
          "Tumor histology: UC", "Tumor histology: UC with variant", "Tumor histology: Pure non-UC", 
          "Clinical T stage: cTa-cT1", "Clinical T stage: cT2", "Clinical T stage: cT3-cT4", "Clinical N stage: cN1", 
          "Procedure: Open", "Procedure: Laparoscopic", "Procedure: Robotic", 
          "LN dissection: Yes", "Urinary diversion: Other",
          "Pathological T stage: pTa-pT1", "Pathological T stage: pT2", "Pathological T stage: pT3-pT4", 
          "Pathological N stage: pN1-pN3", 
          "Tumor grade: High grade", 
          "Neoadjuvant chemotherapy: Yes", "Adjuvant chemotherapy: Yes", 
          "Hemoglobin", "NLR (log-transformed)",
          "Risk of urethral recurrence: High risk")
)

# standardized mean difference plot
smdplot <- love.plot(wgt_out, drop.distance = TRUE, abs = FALSE, 
                     title = "Covariate balance", colors = c("#EE0000FF", "#39568CFF"), shapes = c("diamond filled", "circle filled"), sample.names = c("Crude", "IPTW-adjusted"), 
                     thresholds = c(m = 0.1), size = 3.3, themes = theme_grey(), var.names = label, binary = "std") + 
  theme(legend.position = "top") + 
  labs(
    title = "Covariate balance (RC with PU versus RC alone)", 
    x = "Standardized mean differences"
  )

# dataset plus weights and pscores
df <- tibble(dimp, weight = wgt_out$weights, pscores = wgt_out$ps) %>% 
  mutate_at(vars(nutr, death_all, death_cancer),
            funs(case_when(. == "yes" ~ 1, TRUE ~ 0)))

# propensity score distributions
# define ggplot2 theme...
theme_juog <- function(base_size = 11, 
                       dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.4)),
          plot.subtitle = element_text(size = rel(1.1)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

psdist <- ggplot() + 
  geom_density(data = filter(df, urethrectomy == "yes"),
               aes(x = pscores, weight = weight, fill = urethrectomy), alpha = 0.2) + 
  geom_density(data = filter(df, urethrectomy == "no"),
               aes(x = pscores, weight = weight, fill = urethrectomy, y = -..density..), alpha = 0.2) + 
  geom_density(data = filter(df, urethrectomy == "yes"),
               aes(x = pscores, fill = urethrectomy), alpha = 0.6) + 
  geom_density(data = filter(df, urethrectomy == "no"),
               aes(x = pscores, fill = urethrectomy, y = -..density..), alpha = 0.6) +
  geom_hline(yintercept = 0, color = "white", size = 0.25) + 
  annotate(geom = "label", x = 0.22, y = 1.5, label = "RC with prophylactic urethrectomy\n(Actual and IPTW pseudo-population)",
           fill = "#EE0000FF", alpha = 0.6, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.42, y = -2.6, label = "RC alone\n(IPTW pseudo-population)",
           fill = "#3B4992FF", alpha = 0.3, color = "white", hjust = 0) +
  annotate(geom = "label", x = 0.3, y = -1.5, label = "RC alone\n(Actual)",
           fill = "#3B4992FF", alpha = 0.6, color = "white", hjust = 0) +
  scale_fill_aaas() + 
  theme_juog() + 
  scale_y_continuous(label = abs) + 
  labs(x = "Estimated probability of receiving prophylactic urethrectomy", 
       y = "Kernel density") # PS density for all patients

# table1 for all patients (crude data)
reset_gtsummary_theme()
tbl_theme <- list(
  "tbl_summary-str:missing_stat" = "{N_miss} ({p_miss})" # display the number of missing values and its percentage
)
set_gtsummary_theme(tbl_theme)

data_tbl <- dfull %>% mutate(imp = "Crude")
data_imp <- df_imp$ximp %>% mutate(imp = "Imputed")
df_tblimp <- rbind(data_tbl, data_imp)

tbl1_crudeimp <- df_tblimp %>%
  dplyr::select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, neck_tur, prostate_tur, cis_tur, hist_rc, 
                ct, cn, type, lnd, div, pt, pn, grade, nac, ac, hgb, lognlr, risk_urethra, urethrectomy, imp) %>%
  tbl_summary(
    by = imp,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urinary cancer",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      neck_tur ~ "Bladder neck involvement",
      prostate_tur ~ "Prostatic urethral tumor involvement",
      cis_tur ~ "Carcinoma in situ",
      hist_rc ~ "Tumor histology",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical procedure",
      lnd ~ "Lymph node dissection",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      grade ~ "Tumor grade",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      risk_urethra ~ "Risk of urethral recurrence", 
      urethrectomy ~ "prophylactic urethrecromy"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) 

tbl1_crude <- dfull %>%
  dplyr::select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, neck_tur, prostate_tur, cis_tur, hist_rc, 
                ct, cn, type, lnd, div, pt, pn, grade, nac, ac, hgb, lognlr, risk_urethra, urethrectomy) %>%
  tbl_summary(
    by = urethrectomy,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urinary cancer",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      neck_tur ~ "Bladder neck involvement",
      prostate_tur ~ "Prostatic urethral tumor involvement",
      cis_tur ~ "Carcinoma in situ",
      hist_rc ~ "Tumor histology",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical procedure",
      lnd ~ "Lymph node dissection",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      grade ~ "Tumor grade",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      risk_urethra ~ "Risk of urethral recurrence"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) %>%
  add_overall()

# table1 for all patients (imputed data)
tbl1_imp <- df %>%
  select(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, neck_tur, prostate_tur, cis_tur, hist_rc, 
         ct, cn, type, lnd, div, pt, pn, grade, nac, ac, hgb, lognlr, risk_urethra, urethrectomy) %>%
  tbl_summary(
    by = urethrectomy,
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urinary cancer",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      neck_tur ~ "Bladder neck involvement",
      prostate_tur ~ "Prostatic urethral tumor involvement",
      cis_tur ~ "Carcinoma in situ",
      hist_rc ~ "Tumor histology",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical procedure",
      lnd ~ "Lymph node dissection",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      grade ~ "Tumor grade",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      risk_urethra ~ "Risk of urethral recurrence"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1))) %>%
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) %>%
  add_overall()

# weighted imputed data
df_weight <- survey::svydesign(ids = ~ id, weights = ~ weight, data = df)
wtbl <- df_weight %>%
  tbl_svysummary(
    by = urethrectomy,
    include = c(age_rc, sex, bmi, ecog_ps_cat, smoking, utuc, nmibc, bcg, multi_tur, neck_tur, prostate_tur, cis_tur, hist_rc, 
                ct, cn, type, lnd, div, pt, pn, grade, nac, ac, hgb, lognlr, risk_urethra),
    label = list(
      age_rc ~ "Age at radical cystectomy, year",
      sex ~ "Sex",
      bmi ~ "Body mass index",
      ecog_ps_cat ~ "ECOG performance status",
      smoking ~ "Smoking status",
      utuc ~ "Previous upper-tract urinary cancer",
      nmibc ~ "Previous non-muscle invasive urothelial cancer",
      bcg ~ "Prior intravesical BCG therapy",
      multi_tur ~ "Multifocality",
      neck_tur ~ "Bladder neck involvement",
      prostate_tur ~ "Prostatic urethral tumor involvement",
      cis_tur ~ "Carcinoma in situ",
      hist_rc ~ "Tumor histology",
      ct ~ "Clinical T stage",
      cn ~ "Clinical N stage",
      type ~ "Surgical procedure",
      lnd ~ "Lymph node dissection",
      div ~ "Urinary tract diversion",
      pt ~ "Pathological T stage",
      pn ~ "Pathological N stage",
      grade ~ "Tumor grade",
      nac ~ "Neoadjuvant chemotherapy",
      ac ~ "Adjuvant chemotherapy",
      hgb ~ "Hemoglobin concentration",
      lognlr ~ "Neutrophil-to-lymphocyte ratio (log-transformed)",
      risk_urethra ~ "Risk of urethral recurrence"
    ),
    type = all_dichotomous() ~ "categorical",
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                     all_categorical() ~ "{p}"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(1, 1))
  ) %>%
  add_difference(test = list(everything() ~ "smd")) %>%
  modify_column_hide(columns = ci) 

tbl1_imputation <- tbl_merge(tbls = list(tbl1_crude, tbl1_imp), 
                             tab_spanner = c("Crude", "After imputation"))

tbl1 <- tbl_merge(tbls = list(tbl1_crude, tbl1_imp, wtbl),
                  tab_spanner = c("Unweighted population (before imputation)", "Unweighted population (after imputation)", "Weighted population"))

# survival metrics
# calculate the median follow-up
library(pec)
fu_time <- quantile(prodlim::prodlim(Hist(fu, death_all) ~ 1, data = df, reverse = TRUE))

# define ggplot2 theme...
theme_km <- function(base_size = 10, 
                     dark_text = "#1A242F") {
  mid_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[2]
  
  light_text <-  monochromeR::generate_palette(
    dark_text, "go_lighter",
    n_colors = 5)[3]
  
  theme_grey(base_size = base_size) +
    theme(text = element_text(color = mid_text, lineheight = 1.1),
          plot.title = element_text(color = dark_text, size = rel(1.2)),
          plot.subtitle = element_text(size = rel(1.0)),
          axis.text.y = element_text(color = light_text, size = rel(1)),
          axis.title.y = element_text(size = rel(1)),
          axis.text.x = element_text(color = mid_text, size = rel(1)),
          axis.title.x = element_text(size = rel(1)),
          legend.position = "top",
          legend.title = element_blank(),
          panel.grid = element_line(color = "#F3F4F5"),
          plot.caption = element_text(size = rel(1))) 
}

# Kaplan-Meier estimates
# stratified Kaplan-Meier curves
rfs_model_fit_crude <- survfit(data = df, Surv(nutrfs, nutr) ~ urethrectomy)
rfs_model_crude <- survminer::ggsurvplot(rfs_model_fit_crude, 
                                         risk.table = FALSE, font.tickslab = 10, font.x = 10, font.y = 10,
                                         title = "NUTRFS", 
                                         subtitle = "Crude",
                                         legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                         legend.title = "", 
                                         censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                         xlab = "Months since radical cystectomy", 
                                         ylab = "NUTRFS probability",
                                         palette = "aaas", size = 0.5,  break.time.by = 12,
                                         ggtheme = theme_grey(), conf.int = TRUE) 

rfs_model_crude$plot <- rfs_model_crude$plot + 
  theme_km() 

# overall survival
css_model_fit_crude <- survival::survfit(data = df, Surv(fu, death_cancer) ~ urethrectomy)
css_model_crude <- survminer::ggsurvplot(css_model_fit_crude, 
                                         risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                         title = "CSS", 
                                         subtitle = "Crude",
                                         legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                         legend.title = "", 
                                         censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                         xlab = "Months since radical cystectomy", 
                                         ylab = "CSS probability",
                                         palette = "aaas", size = 0.5,  break.time.by = 12,
                                         ggtheme = theme_grey(), risk.table.title = "Number at risk",
                                         risk.table.col = "strata",
                                         tables.height = 0.15, risk.table.fontsize = 3,
                                         conf.int = TRUE) 

css_model_crude$plot <- css_model_crude$plot + 
  theme_km() 

# overall survival
os_model_fit_crude <- survival::survfit(data = df, Surv(fu, death_all) ~ urethrectomy)
os_model_crude <- survminer::ggsurvplot(os_model_fit_crude, 
                                        risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                        title = "OS", 
                                        subtitle = "Crude",
                                        legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                        legend.title = "", 
                                        censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                        xlab = "Months since radical cystectomy", 
                                        ylab = "OS probability",
                                        palette = "aaas", size = 0.5,  break.time.by = 12,
                                        ggtheme = theme_grey(), conf.int = TRUE) 

os_model_crude$plot <- os_model_crude$plot + 
  theme_km() 

rfs_model_fit_iptw <- survfit2(data = df, Surv(nutrfs, nutr) ~ urethrectomy, weights = weight)
rfs_model_iptw <- survminer::ggsurvplot(rfs_model_fit_iptw, 
                                        risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                        title = "NUTRFS", 
                                        subtitle = "IPTW-adjusted",
                                        legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                        legend.title = "", 
                                        censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                        xlab = "Months since radical cystectomy", 
                                        ylab = "NUTRFS probability",
                                        palette = "aaas", size = 0.5,  break.time.by = 12,
                                        ggtheme = theme_grey(), conf.int = TRUE) 

rfs_model_iptw$plot <- rfs_model_iptw$plot + 
  theme_km() 

# overall survival
css_model_fit_iptw <- survival::survfit(data = df, Surv(fu, death_cancer) ~ urethrectomy, weights = weight)
css_model_iptw <- survminer::ggsurvplot(css_model_fit_iptw, 
                                        risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                        title = "CSS", 
                                        subtitle = "IPTW-adjusted",
                                        legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                        legend.title = "", 
                                        censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                        xlab = "Months since radical cystectomy", 
                                        ylab = "CSS probability",
                                        palette = "aaas", size = 0.5,  break.time.by = 12,
                                        ggtheme = theme_grey(), conf.int = TRUE) 

css_model_iptw$plot <- css_model_iptw$plot + 
  theme_km() 

# overall survival
os_model_fit_iptw <- survival::survfit(data = df, Surv(fu, death_all) ~ urethrectomy, weights = weight)
os_model_iptw <- survminer::ggsurvplot(os_model_fit_iptw, 
                                       risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                       title = "OS", 
                                       subtitle = "IPTW-adjusted",
                                       legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                       legend.title = "", 
                                       censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                       xlab = "Months since radical cystectomy", 
                                       ylab = "OS probability",
                                       palette = "aaas", size = 0.5,  break.time.by = 12,
                                       ggtheme = theme_grey(), conf.int = TRUE) 

os_model_iptw$plot <- os_model_iptw$plot + 
  theme_km() 

merge <- (rfs_model_crude$plot + css_model_crude$plot + os_model_crude$plot) / (rfs_model_iptw$plot + css_model_iptw$plot + os_model_iptw$plot) + 
  plot_annotation(caption = "Data are weighted proportions and not absolute numbers.")

# weighted log-rank test using RISCA package
df <- df %>% 
  mutate(urethrectomy = case_when(urethrectomy == "yes" ~ 1, TRUE ~ 0))

pval_rfs <- RISCA::ipw.log.rank(times = df$nutrfs, 
                                failures = df$nutr,
                                variable = df$urethrectomy,
                                weights = df$weight)
pval_css <- RISCA::ipw.log.rank(times = df$fu, 
                                failures = df$death_cancer,
                                variable = df$urethrectomy,
                                weights = df$weight)
pval_os <- RISCA::ipw.log.rank(times = df$fu, 
                               failures = df$death_all,
                               variable = df$urethrectomy,
                               weights = df$weight)

# crude Cox regression model
theme_gtsummary_journal(journal = "jama")
cox_rfs <- coxph(Surv(nutrfs, nutr) ~ urethrectomy, data = df) %>% 
  tbl_regression(exponentiate = TRUE)
cox_css <- coxph(Surv(fu, death_cancer) ~ urethrectomy, data = df) %>% 
  tbl_regression(exponentiate = TRUE)
cox_os <- coxph(Surv(fu, death_all) ~ urethrectomy, data = df) %>% 
  tbl_regression(exponentiate = TRUE)

# merge table
cox_crude_merge <- tbl_merge(tbls = list(cox_rfs, cox_css, cox_os),
                             tab_spanner = c("NUTRFS", "CSS", "OS")
)

# IPTW-adjusted Cox regression model
theme_gtsummary_journal(journal = "jama")
cox_rfs <- coxph(Surv(nutrfs, nutr) ~ urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)
cox_css <- coxph(Surv(fu, death_cancer) ~ urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)
cox_os <- coxph(Surv(fu, death_all) ~ urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)

# merge table
cox_iptw_merge <- tbl_merge(tbls = list(cox_rfs, cox_css, cox_os),
                            tab_spanner = c("NUTRFS", "CSS", "OS")
)

# check the proportional hazards assumption
# crude
coxph(Surv(nutrfs, nutr) ~ urethrectomy, data = df) %>% cox.zph()
coxph(Surv(fu, death_cancer) ~ urethrectomy, data = df) %>% cox.zph()
coxph(Surv(fu, death_all) ~ urethrectomy, data = df) %>% cox.zph()

# IPTW-adjusted
coxph(Surv(nutrfs, nutr) ~ urethrectomy, data = df, weights = weight) %>% cox.zph()
coxph(Surv(fu, death_cancer) ~ urethrectomy, data = df, weights = weight) %>% cox.zph()
coxph(Surv(fu, death_all) ~ urethrectomy, data = df, weights = weight) %>% cox.zph()

# visualized RMSTs...
# calculate last observed or censored survival time
os <- survfit(Surv(fu, death_all) ~ urethrectomy, data = df, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

css <- survfit(Surv(fu, death_cancer) ~ urethrectomy, data = df, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

rfs <- survfit(Surv(nutrfs, nutr) ~ urethrectomy, data = df, weights = weight) %>%
  broom::tidy() %>% dplyr::group_by(strata) %>% summarize(max = max(time)) %>% summarize(time = min(max))

# calculate IPTW-adjusted RMSTs
rmst_rfs <- akm_rmst(time = df$nutrfs, 
                     status = df$nutr, 
                     group = as.factor(df$urethrectomy), 
                     weight = df$weight, 
                     tau = 107) 

rmst_css <- akm_rmst(time = df$fu, 
                     status = df$death_cancer, 
                     group = as.factor(df$urethrectomy), 
                     weight = df$weight, 
                     tau = 107) 

rmst_os <- akm_rmst(time = df$fu, 
                    status = df$death_all, 
                    group = as.factor(df$urethrectomy), 
                    weight = df$weight, 
                    tau = 107) 

# calculate crude RMSTs
rmst_rfs_crude <- survRM2::rmst2(time = df$nutrfs,
                                 status = df$nutr,
                                 arm = as.factor(df$urethrectomy),
                                 tau = 106)

rmst_css_crude <- survRM2::rmst2(time = df$fu,
                                 status = df$death_cancer,
                                 arm = as.factor(df$urethrectomy),
                                 tau = 106)

rmst_os_crude <- survRM2::rmst2(time = df$fu,
                                status = df$death_all,
                                arm = as.factor(df$urethrectomy),
                                tau = 106)

# IPTW-adjusted logistic regression analysis for predicting Clavien-Dindo grade 3 complication
alog <- df %>% glm(comp_g3 ~ urethrectomy, data = ., weights = weight, family = binomial) %>% 
  tbl_regression(exponentiate = TRUE,
                 label = list(urethrectomy ~ "prophylactic urethrectomy"))

# estimated blood loss
debl <- d %>% select(id, ebl)
d_ebl <- df %>%
  left_join(debl, by = "id") %>% 
  filter(!is.na(ebl)) 

ebl_crude <- d_ebl %>% 
  dplyr::select(ebl, comp_g3, urethrectomy) %>% 
  tbl_summary(
    by = urethrectomy,
    label = list(
      ebl ~ "Estimated blood loss",
      comp_g3 ~ "Severe surgical complication"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_p() %>% 
  add_overall()

# weighted imputed data
svyebl <- survey::svydesign(ids = ~ id, weights = ~ weight, data = d_ebl)
ebl_iptw <- svyebl %>%
  tbl_svysummary(
    by = urethrectomy,
    include = c(ebl, comp_g3),
    label = list(
      ebl ~ "Estimated blood loss",
      comp_g3 ~ "Severe surgical complication"
    ),
    statistic = list(all_continuous() ~ "{median} ({p25}-{p75})",
                     all_categorical() ~ "{n} ({p})"),
    digits = list(all_continuous() ~ c(1, 1),
                  all_categorical() ~ c(0, 1)),
    missing_text = "Missing") %>%
  add_p() 

tbl_ae <- tbl_merge(tbls = list(ebl_crude, ebl_iptw),
                    tab_spanner = c("Crude", "IPTW-adjusted"))


# Heterogeneity of treatment effect...
theme_gtsummary_journal(journal = "jama")
# urethral risk group
hte_rfs_risk <- coxph(Surv(nutrfs, nutr) ~ risk_urethra * urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)

# sex
hte_rfs_sex <- coxph(Surv(nutrfs, nutr) ~ sex * urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)
hte_rfs_neck <- coxph(Surv(nutrfs, nutr) ~ neck_tur * urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)
hte_rfs_prostate <- coxph(Surv(nutrfs, nutr) ~ prostate_tur * urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)
hte_rfs_multi <- coxph(Surv(nutrfs, nutr) ~ multi_tur * urethrectomy, data = df, weights = weight) %>% 
  tbl_regression(exponentiate = TRUE)

# merge table
hte_merge <- tbl_stack(list(hte_rfs_sex, hte_rfs_prostate, hte_rfs_risk, hte_rfs_nac),
                       group_header = c("Sex - Urethrectomy interaction", "Prostatic involvement - Urethrectomy interaction", "Risk - Urethrectomy interaction", "NAC - Urethrectomy interaction")
)

# subgroup analysis
rfs_model_fit_iptw_low <- survival::survfit(data = df %>% filter(risk_urethra == "Low risk"), Surv(nutrfs, nutr) ~ urethrectomy, weights = weight)
rfs_model_iptw_low <- survminer::ggsurvplot(rfs_model_fit_iptw_low, 
                                            risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                            title = "NUTRFS", 
                                            subtitle = "Patients at low-risk of urethral recurrence",
                                            legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                            legend.title = "", 
                                            censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                            xlab = "Months since radical cystectomy", 
                                            ylab = "NUTRFS probability",
                                            palette = "aaas", size = 0.5,  break.time.by = 12,
                                            ggtheme = theme_grey(), conf.int = TRUE) 

rfs_model_iptw_low$plot <- rfs_model_iptw_low$plot + 
  theme_km() 

rfs_model_fit_iptw_low <- survival::survfit(data = df %>% filter(risk_urethra == "Low risk"), Surv(nutrfs, nutr) ~ urethrectomy, weights = weight)
rfs_model_iptw_low <- survminer::ggsurvplot(rfs_model_fit_iptw_low, 
                                            risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                            title = "NUTRFS", 
                                            subtitle = "Patients at low-risk of urethral recurrence",
                                            legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                            legend.title = "", 
                                            censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                            xlab = "Months since radical cystectomy", 
                                            ylab = "NUTRFS probability",
                                            palette = "aaas", size = 0.5,  break.time.by = 12,
                                            ggtheme = theme_grey(), conf.int = TRUE) 

rfs_model_iptw_low$plot <- rfs_model_iptw_low$plot + 
  theme_km() 

# overall survival
css_model_fit_iptw_low <- survival::survfit(data = df %>% filter(risk_urethra == "Low risk"), Surv(fu, death_cancer) ~ urethrectomy, weights = weight)
css_model_iptw_low <- survminer::ggsurvplot(css_model_fit_iptw_low, 
                                            risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                            title = "CSS", 
                                            subtitle = "Patients at low-risk of urethral recurrence",
                                            legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                            legend.title = "", 
                                            censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                            xlab = "Months since radical cystectomy", 
                                            ylab = "CSS probability",
                                            palette = "aaas", size = 0.5,  break.time.by = 12,
                                            ggtheme = theme_grey(), conf.int = TRUE) 

css_model_iptw_low$plot <- css_model_iptw_low$plot + 
  theme_km() 

# overall survival
os_model_fit_iptw_low <- survival::survfit(data = df %>% filter(risk_urethra == "Low risk"), Surv(fu, death_all) ~ urethrectomy, weights = weight)
os_model_iptw_low <- survminer::ggsurvplot(os_model_fit_iptw_low, 
                                           risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                           title = "OS", 
                                           subtitle = "Patients at low-risk of urethral recurrence",
                                           legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                           legend.title = "", 
                                           censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                           xlab = "Months since radical cystectomy", 
                                           ylab = "OS probability",
                                           palette = "aaas", size = 0.5,  break.time.by = 12,
                                           ggtheme = theme_grey(), conf.int = TRUE) 

os_model_iptw_low$plot <- os_model_iptw_low$plot + 
  theme_km() 

# subgroup analysis
rfs_model_fit_iptw_high <- survival::survfit(data = df %>% filter(risk_urethra == "High risk"), Surv(nutrfs, nutr) ~ urethrectomy, weights = weight)
rfs_model_iptw_high <- survminer::ggsurvplot(rfs_model_fit_iptw_high, 
                                             risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                             title = "NUTRFS", 
                                             subtitle = "Patients at high-risk of urethral recurrence",
                                             legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                             legend.title = "", 
                                             censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                             xlab = "Months since radical cystectomy", 
                                             ylab = "NUTRFS probability",
                                             palette = "aaas", size = 0.5,  break.time.by = 12,
                                             ggtheme = theme_grey(), conf.int = TRUE) 

rfs_model_iptw_high$plot <- rfs_model_iptw_high$plot + 
  theme_km() 

# overall survival
css_model_fit_iptw_high <- survival::survfit(data = df %>% filter(risk_urethra == "High risk"), Surv(fu, death_cancer) ~ urethrectomy, weights = weight)
css_model_iptw_high <- survminer::ggsurvplot(css_model_fit_iptw_high, 
                                             risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                             title = "CSS", 
                                             subtitle = "Patients at high-risk of urethral recurrence",
                                             legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                             legend.title = "", 
                                             censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                             xlab = "Months since radical cystectomy", 
                                             ylab = "CSS probability",
                                             palette = "aaas", size = 0.5,  break.time.by = 12,
                                             ggtheme = theme_grey(), conf.int = TRUE) 

css_model_iptw_high$plot <- css_model_iptw_high$plot + 
  theme_km() 

# overall survival
os_model_fit_iptw_high <- survival::survfit(data = df %>% filter(risk_urethra == "High risk"), Surv(fu, death_all) ~ urethrectomy, weights = weight)
os_model_iptw_high <- survminer::ggsurvplot(os_model_fit_iptw_high, 
                                            risk.table = TRUE, font.tickslab = 10, font.x = 10, font.y = 10,
                                            title = "OS", 
                                            subtitle = "Patients at high-risk of urethral recurrence",
                                            legend = "top", legend.labs = c("RC alone", "RC with PU"),
                                            legend.title = "", 
                                            censor = TRUE, censor.shape = "|", censor.size = 1.5,
                                            xlab = "Months since radical cystectomy", 
                                            ylab = "OS probability",
                                            palette = "aaas", size = 0.5,  break.time.by = 12,
                                            ggtheme = theme_grey(), conf.int = TRUE) 

os_model_iptw_high$plot <- os_model_iptw_high$plot + 
  theme_km() 

merge_subgroup <- (rfs_model_iptw_low$plot + css_model_iptw_low$plot + os_model_iptw_low$plot) / (rfs_model_iptw_high$plot + css_model_iptw_high$plot + os_model_iptw_high$plot) + 
  plot_annotation(caption = "IPTW-adjusted Kaplan-Meier curves. Data are weighted proportions and not absolute numbers.")

theme_gtsummary_journal(journal = "jama")

cox_rfs <- dimp %>% 
  mutate(nutr = case_when(nutr == "yes" ~ 1, TRUE ~ 0)) %>% 
  coxph(data = ., 
        Surv(nutrfs, nutr) ~ age_rc + sex + bmi + ecog_ps_cat + smoking + utuc + nmibc + bcg + multi_tur + neck_tur + prostate_tur + cis_tur + hist_rc + 
          ct + cn + type + lnd + div + pt + pn + grade + nac + ac + hgb + lognlr + risk_urethra + urethrectomy) %>% 
  cox.zph()

cox_css <- dimp %>% 
  mutate(death_cancer = case_when(death_cancer == "yes" ~ 1, TRUE ~ 0)) %>% 
  coxph(data = ., 
        Surv(fu, death_cancer) ~ age_rc + sex + bmi + ecog_ps_cat + smoking + utuc + nmibc + bcg + multi_tur + neck_tur + prostate_tur + cis_tur + hist_rc + 
          ct + cn + type + lnd + div + pt + pn + grade + nac + ac + hgb + lognlr + risk_urethra + urethrectomy) %>% 
  cox.zph()

cox_os <- dimp %>% 
  mutate(death_all = case_when(death_all == "yes" ~ 1, TRUE ~ 0)) %>% 
  coxph(data = ., 
        Surv(fu, death_all) ~ age_rc + sex + bmi + ecog_ps_cat + smoking + utuc + nmibc + bcg + multi_tur + neck_tur + prostate_tur + cis_tur + hist_rc + 
          ct + cn + type + lnd + div + pt + pn + grade + nac + ac + hgb + lognlr + risk_urethra + urethrectomy) %>% 
  cox.zph()

# subgroup analysis
subos <- function(var, data) {
  theme_gtsummary_journal(journal = "jama")
  data %>%
    select(fu, death_all, urethrectomy, weight, var) %>%
    tbl_strata(
      strata = var,
      .tbl_fun =
        ~ coxph(Surv(fu, death_all) ~ urethrectomy, data = ., weights = weight) %>% 
        tbl_regression(exponentiate = TRUE) %>% 
        add_n(),
      .combine_with = "tbl_stack"
    )
} 

# sex
os_sex <- subos(var = "sex", data = df)
os_neck <- subos(var = "neck_tur", data = df)
os_multi <- subos(var = "multi_tur", data = df)
os_prostate <- subos(var = "prostate_tur", data = df)
os_risk <- subos(var = "risk_urethra", data = df)

# patients with urethral recurrence
df_urethra <- d %>% 
  filter(urethra_rec == "yes") %>% 
  mutate(time_urethra_trt = time_length(date_surgery %--% date_urethra, "months") %>% round(digits = 1),
         fu = fu %>% round(digits = 1)) %>%
  select(sex, time_urethra_trt, urethra_trt, nutrfs, nutr, fu, status) %>% 
  arrange(time_urethra_trt)
