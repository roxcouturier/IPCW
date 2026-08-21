library(dplyr)
library(purrr)
library(tidyr)
library(mice)
library(survival)
library(broom)
library(ipw)
library(cobalt)
# 
# load("/Users/roxanecouturier/Desktop/Doctorat/Article 3 Version 2/graaph2023.RData")
# 
# #### 1. SÉLECTION ET TYPAGE DES VARIABLES ####
# 
# data <- graaph.2023 %>%
#   select(
#     greffe, delallo, Type.bmt, BRAS, randodt,
#     dc, survie, delrec, rec, suivi, RANDOMIZATION_R1,
#     age, SEXE, bmi, POIDSDG, TAILLEDG, CNSPRECOCE, CNSTARDIF, cNS, cNSTARDIF, PS, GBJ1,
#     CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
#     SdTUMORAL, MEDIASTINDG, TLP, TLPPRECOCE, TLPTARDIF,
#     GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
#     ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
#     TRANSCRIT, STRATE, IKZF1, NO_IKZF1, RES_IKZF1,
#     BLASTESJ1, CSR, IMMUNOJ1, IDBLASTPREPHASE,
#     MRD1FAIT, MRD1LAB, MRD2FAIT, mrd1dt, mrd2dt, mrd3dt, mrd4dt,
#     mrd1, mrd2, mrd3, mrd4
#   )
# 
# 
# data$PS         <- factor(data$PS)
# data$SdTUMORAL  <- as.factor(data$SdTUMORAL)
# data$PLQDIAG    <- as.numeric(data$PLQDIAG)
# data$GBJ1       <- as.numeric(data$GBJ1)
# data$BLASTESJ1  <- as.numeric(data$BLASTESJ1)
# data$GBDIAG     <- as.numeric(data$GBDIAG)
# data$BLASTESDIAG<- as.numeric(data$BLASTESDIAG)
# 
# data <- data %>%
#   mutate(
#     PS_regroup = case_when(
#       PS %in% c(2, 3) ~ "2+",
#       TRUE ~ as.character(PS)
#     ) %>% factor()
#   )
# 
# #### 2. VARIABLES DÉRIVÉES ####
# 
# data <- data %>%
#   mutate(
#     survie_month = survie / 30.4167,
#     allo         = ifelse(is.na(delallo), 0, 1),
#     trt          = ifelse(BRAS == 1, 1, 0),
#     pfs          = case_when(rec == 1 ~ 1, dc == 1 ~ 1, TRUE ~ 0),
#     pfs_time     = case_when(rec == 1 ~ delrec, dc == 1 ~ survie, TRUE ~ suivi) / 30.4167
#   ) %>%
#   mutate(
#     randodt = as.Date(randodt),
#     mrd1dt  = as.Date(mrd1dt),
#     mrd2dt  = as.Date(mrd2dt),
#     mrd3dt  = as.Date(mrd3dt),
#     mrd4dt  = as.Date(mrd4dt)
#   ) %>%
#   mutate(
#     mrd1_month = as.numeric(mrd1dt - randodt) / 30.4167,
#     mrd2_month = as.numeric(mrd2dt - randodt) / 30.4167,
#     mrd3_month = as.numeric(mrd3dt - randodt) / 30.4167,
#     mrd4_month = as.numeric(mrd4dt - randodt) / 30.4167
#   ) %>%
#   mutate(
#     age_cat40 = factor(ifelse(age < 40, "< 40", ">= 40"), levels = c("< 40", ">= 40")),
#     age_cat50 = factor(ifelse(age < 50, "< 50", ">= 50"), levels = c("< 50", ">= 50"))
#   ) %>%
#   mutate(patient_id = row_number())
# 
# 
# 
# #### 3. SUPPRESSION DES DATES BRUTES AVANT EXPORT (anonymisation) ####
# data_export <- data %>%
#   select(-randodt, -mrd1dt, -mrd2dt, -mrd3dt, -mrd4dt)
# 
# # Vérifier qu'il ne reste plus de colonne de type Date
# sapply(data_export, class)[sapply(data_export, function(x) inherits(x, "Date"))]
# # doit renvoyer un vecteur vide (named character(0))
# 
# #### 4. SAUVEGARDE (version sans dates, prête à mettre en ligne) ####
# saveRDS(data_export, file = "/Users/roxanecouturier/Desktop/Doctorat/Article 3 Version 2/data_graaph_anonymise.rds")

data <- readRDS("/Users/roxanecouturier/Desktop/Doctorat/Article 3 Version 2/data_graaph_anonymise.rds")
file.exists("/Users/roxanecouturier/Desktop/Doctorat/Article 3 Version 2/data_graaph_anonymise.rds")  # TRUE attendu

#### 3. IMPUTATION (mice) ####

imp_data <- data %>%
  select(
    patient_id,
    cNS, BLASTESJ1, GBJ1, PS, PS_regroup, SdTUMORAL,
    age, age_cat40, age_cat50, SEXE, bmi, POIDSDG, TAILLEDG,
    CNSPRECOCE, CNSTARDIF, cNSTARDIF,
    CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
    MEDIASTINDG, TLP, TLPPRECOCE, TLPTARDIF,
    GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
    ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
    TRANSCRIT, STRATE, IKZF1, RES_IKZF1,NO_IKZF1,
    CSR, IMMUNOJ1, IDBLASTPREPHASE,
    mrd1, dc, survie_month, allo, pfs, pfs_time, mrd1_month
  )

pred_matrix <- make.predictorMatrix(imp_data)
pred_matrix[, "patient_id"] <- 0
pred_matrix["patient_id", ] <- 0

meth <- make.method(imp_data)
meth["patient_id"] <- ""

imp <- mice(imp_data, m = 20, method = meth, predictorMatrix = pred_matrix, seed = 123)
# plot(imp)  # à décommenter pour vérifier la convergence

imputed_vars <- complete(imp, 1) %>%
  select(patient_id, cNS, BLASTESJ1, GBJ1, PS, PS_regroup, SdTUMORAL,
         MEDIASTINDG, GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
         TLPTARDIF, RES_IKZF1, TLPPRECOCE, TLP, NO_IKZF1)   # ajouté

data <- data %>%
  select(-cNS, -BLASTESJ1, -GBJ1, -PS, -PS_regroup, -SdTUMORAL,
         -MEDIASTINDG, -GBDIAG, -PNNDIAG, -PLQDIAG, -BLASTESDIAG,
         -TLPTARDIF, -RES_IKZF1, -TLPPRECOCE, -TLP, -NO_IKZF1) %>%   # ajouté
  left_join(imputed_vars, by = "patient_id")

#### 4. FORMAT LONG : mrd change de valeur à mrd1_month ####

make_long_mrd1 <- function(patient_id, survie_month, dc, mrd1_month, mrd1) {
  if (is.na(mrd1_month)) {
    tibble(
      patient_id = patient_id, tstart = 0, tstop = survie_month,
      mrd = 0, event = ifelse(dc == 1, 1, 0)
    )
  } else {
    bind_rows(
      tibble(patient_id = patient_id, tstart = 0, tstop = mrd1_month, mrd = 0, event = 0),
      tibble(patient_id = patient_id, tstart = mrd1_month, tstop = survie_month,
             mrd = mrd1, event = ifelse(dc == 1, 1, 0))
    )
  }
}

data_long <- pmap_dfr(
  list(patient_id = data$patient_id, survie_month = data$survie_month,
       dc = data$dc, mrd1_month = data$mrd1_month, mrd = data$mrd1),
  make_long_mrd1
)

data_long <- data_long %>%
  left_join(
    data %>%
      select(
        patient_id, age, SEXE, bmi, POIDSDG, TAILLEDG, delallo, BRAS,
        PS, PS_regroup, cNS, age_cat40, age_cat50, GBJ1, mrd1_month,
        CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
        SdTUMORAL, MEDIASTINDG, CNSPRECOCE, TLP, TLPPRECOCE, TLPTARDIF,
        GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
        ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
        TRANSCRIT, STRATE, IKZF1, NO_IKZF1, RES_IKZF1,
        BLASTESJ1, CSR, IMMUNOJ1, IDBLASTPREPHASE
      ),
    by = "patient_id"
  ) %>%
  mutate(
    allo_event = ifelse(!is.na(delallo) & delallo <= tstop, 1, 0),
    event1     = ifelse(!is.na(delallo) & delallo <= tstop, 0, event),
    tstop      = ifelse(!is.na(delallo) & delallo <= tstop, delallo, tstop),
    trt        = ifelse(BRAS == 1, 1, 0)
  )

#### 5. SCREENING UNIVARIÉ : effet des covariables sur OS et allogreffe ####

covariates <- c(
  "age", "SEXE", "bmi", "POIDSDG", "TAILLEDG", "trt",
  "CENTRE", "CENTREBMT", "FRATRIE", "FRATRIEBMT",
  "SdTUMORAL", "MEDIASTINDG", "CNSPRECOCE",
  "TLP", "TLPPRECOCE", "TLPTARDIF",
  "GBDIAG", "PNNDIAG", "PLQDIAG", "BLASTESDIAG",
  "ALBUMINEDIAG", "LDHDIAG", "ASATDIAG", "ALATDIAG", "PS_regroup",
  "CREATDIAG", "BILITDIAG",
  "TRANSCRIT", "STRATE", "IKZF1", "NO_IKZF1", "RES_IKZF1",
  "BLASTESJ1", "CSR", "IMMUNOJ1", "IDBLASTPREPHASE", "GBJ1", "cNS",
  "age_cat40", "age_cat50", "PS", "mrd"
)

run_univariate <- function(outcome, covariates, data) {
  map_dfr(covariates, function(cov) {
    f <- as.formula(paste0("Surv(tstart, tstop, ", outcome, ") ~ ", cov, " + cluster(patient_id)"))
    fit <- tryCatch(coxph(f, data = data), error = function(e) NULL)
    if (is.null(fit)) {
      return(tibble(covariate = cov, term = NA, HR = NA, lower = NA, upper = NA,
                    p.value = NA, note = "erreur/non convergent"))
    }
    tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      mutate(covariate = cov, .before = 1) %>%
      select(covariate, term, HR = estimate, lower = conf.low, upper = conf.high, p.value)
  })
}

res_os   <- run_univariate("event", covariates, data_long)
res_allo <- run_univariate("allo_event", covariates, data_long)

sig_os   <- res_os   %>% filter(!is.na(p.value), p.value < 0.05) %>% arrange(p.value)
sig_allo <- res_allo %>% filter(!is.na(p.value), p.value < 0.05) %>% arrange(p.value)

print(sig_os,   n = 40)
print(sig_allo, n = 40)

#### 6. ÉCLATEMENT DU JEU DE DONNÉES À TOUS LES TEMPS (pour IPCW) ####

data_long_df <- as.data.frame(data_long) %>% arrange(patient_id, tstart)

# marqueur de censure administrative (avant éclatement)
data_long_df <- data_long_df %>%
  arrange(patient_id, tstart) %>%
  group_by(patient_id) %>%
  mutate(is_last = row_number() == n(),
         censored = ifelse(is_last & event1 == 0, 1, 0)) %>%
  ungroup() %>%
  as.data.frame()

times <- sort(unique(data_long_df$tstop))

data.long <- survSplit(
  data_long_df, cut = times, start = "tstart", end = "tstop", event = "event1"
)

data.long.c <- survSplit(
  data_long_df, cut = times, start = "tstart", end = "tstop", event = "allo_event"
)
data.long$censored <- data.long.c$allo_event

#### 7. POIDS IPCW (exposure = allogreffe, traitée comme censure informative) ####

w_full <- ipwtm(
  exposure    = censored,
  family      = "survival",
  numerator   = ~ 1,
# denominator = ~  SdTUMORAL +GBJ1 +BLASTESDIAG +trt +mrd +MEDIASTINDG +MEDIASTINDG:mrd+ MEDIASTINDG:trt+ SdTUMORAL:mrd + BLASTESDIAG:mrd + GBJ1:mrd+
#   SdTUMORAL:trt + GBJ1:trt + BLASTESDIAG:trt+mrd:trt+ GBJ1:BLASTESDIAG + GBJ1:MEDIASTINDG,
denominator = ~  SdTUMORAL  +BLASTESJ1 +trt +  NO_IKZF1 +mrd, #MEDIASTINDG:mrd+ MEDIASTINDG:trt+ SdTUMORAL:mrd + BLASTESJ1:mrd + 
 #SdTUMORAL:trt + BLASTESJ1:trt+mrd:trt + BLASTESJ1:SdTUMORAL1,
  id          = patient_id,
  tstart      = tstart,
  timevar     = tstop,
  type        = "cens",
  data        = data.long
)$ipw.weights

summary(w_full)
sd(w_full)


library(ggplot2)

data.long$ipw_full <- w_full

ggplot(data.long, aes(x = ipw_full)) +
  geom_histogram(
    bins = 50,
    color = "black",
    fill = "#0072B2"
  ) +
  labs(
    title = "Distribution of IPCW weights",
    x = "IPCW weights",
    y = "Number of intervals"
  ) +
  theme_minimal()


data.long$ipw_full  <- w_full
data.long$ipw_trunc <- pmin(w_full, quantile(w_full, 0.99, na.rm = TRUE))

# # #### 8. VÉRIFICATION DE L'ÉQUILIBRE : bal.tab / love.plot ####
# # 
#  data_last <- data.long %>%
#    arrange(patient_id, tstart) %>%
#    group_by(patient_id) %>%
#    slice_tail(n = 1) %>%
#    ungroup()
# # #
#  table(data_last$allo_event, useNA = "always")
#  summary(data_last$ipw_full)
# # 
#  covariates_model <- c( "mrd","SdTUMORAL", "BLASTESJ1","NO_IKZF1")
# # #
#  bt <- bal.tab(
#   data_last[, covariates_model],
#   treat   = data_last$allo_event,
#   weights = data_last$ipw_full,
#    method  = "weighting",
#    un      = TRUE
#  )
#  bt
# # #
#  love.plot(
#   bt,
#    stat         = "mean.diffs",
#    threshold    = 0.1,
#    abs          = TRUE,
#    var.order    = "unadjusted",
#   colors       = c("#D55E00", "#0072B2"),
#   shapes       = c(16, 17),
#    sample.names = c("Before weighting", "After IPCW (split at all times)"),
#    title        = "Covariate balance across alloSCT status",
#   line         = FALSE
# )

#### 9. MODÈLES DE COX PONDÉRÉS PAR IPW (effet du bras sur l'OS) ####

cox_unweighted <- coxph(Surv(tstart, tstop, event1) ~ trt + cluster(patient_id),
                        data = data.long)
summary(cox_unweighted)

cox_ipw <- coxph(Surv(tstart, tstop, event1) ~ trt + cluster(patient_id),
                 data = data.long, weights = ipw_full)
summary(cox_ipw)


table(data.long$event1)


# 
# library(survival)
# library(purrr)
# library(dplyr)
# library(tibble)
# 
# vars <- c("SdTUMORAL", "BLASTESDIAG", "mrd","NO_IKZF1","trt")
# # 
# # 
#  fit_main <- coxph(
#    Surv(tstart, tstop, allo_event) ~ SdTUMORAL  +BLASTESDIAG +trt +mrd +NO_IKZF1,
#    data = data.long
#  )
#  
#  pairs <- combn(vars, 2, simplify = FALSE)
# # 
#  interaction_tests <- map_dfr(pairs, function(p) {
#    f_int <- as.formula(
#      paste0("Surv(tstart, tstop, allo_event) ~ SdTUMORAL  +BLASTESDIAG +trt +mrd +NO_IKZF1+",
#             p[1], ":", p[2])
#    )
#    fit_int <- tryCatch(coxph(f_int, data = data.long), error = function(e) NULL)
#      if (is.null(fit_int)) {
#      return(tibble(var1 = p[1], var2 = p[2], p_value = NA_real_, note = "erreur"))
#    }
# #   
#    test <- anova(fit_main, fit_int, test = "LRT")
#    # Récupère la dernière colonne (toujours la p-value, quel que soit son nom exact)
#    pval <- test[2, ncol(test)]
# #   
#    tibble(var1 = p[1], var2 = p[2], p_value = as.numeric(pval))
#  })
# # 
#  interaction_tests %>% arrange(p_value)
# 
# 
# 
# data.long$GBJ1_cat <- cut(
#   data.long$GBJ1,
#   breaks = quantile(data.long$GBJ1, 
#                     probs = c(0, .25, .5, .75, 1),
#                     na.rm = TRUE),
#   include.lowest = TRUE
# )
# 
# fit_cat <- coxph(
#   Surv(tstart, tstop, allo_event) ~ 
#     GBJ1_cat + trt + mrd + SdTUMORAL,
#   data = data.long
# )
# 
# summary(fit_cat)
# 
# 
# library(cobalt)
# library(dplyr)
# 
# covariates_model <- c("mrd", "SdTUMORAL", "BLASTESJ1", "NO_IKZF1")
# 
# temps <- sort(unique(data.long$tstop))
# 
# smd_list <- list()
# effectif_list <- list()
# 
# for(t in temps){
#   
#   data_t <- data.long %>%
#     filter(tstop == t)
#   
#   # Effectifs dans chaque groupe
#   effectifs <- table(data_t$allo_event)
#   
#   effectif_list[[as.character(t)]] <- data.frame(
#     time = t,
#     n_non_allo = ifelse("0" %in% names(effectifs), effectifs["0"], 0),
#     n_allo = ifelse("1" %in% names(effectifs), effectifs["1"], 0)
#   )
#   
#   # Vérifier qu'il reste les deux groupes
#   if(length(unique(data_t$allo_event)) == 2){
#     
#     bt <- bal.tab(
#       data_t[, covariates_model],
#       treat   = data_t$allo_event,
#       weights = data_t$ipw_full,
#       method  = "weighting",
#       un      = TRUE
#     )
#     
#     smd_list[[as.character(t)]] <- as.data.frame(bt$Balance) %>%
#       mutate(time = t)
#   }
# }
# 
# # Tableau des effectifs par temps
# effectifs_par_temps <- do.call(rbind, effectif_list)
# 
# effectifs_par_temps
# 
# 
# 
# covariates_model <- c("mrd", "SdTUMORAL", "BLASTESJ1", "NO_IKZF1")
# 
# # Calcul des SMD pour chaque temps
# smd_par_temps <- map2_dfr(
#   data_par_temps_filtre,
#   names(data_par_temps_filtre),
#   ~ {
#     
#     bt <- bal.tab(
#       .x[, covariates_model],
#       treat   = .x$allo_event,
#       weights = .x$ipw_full,
#       method  = "weighting",
#       un      = TRUE
#     )
#     
#     as.data.frame(bt$Balance) %>%
#       rownames_to_column("covariate") %>%
#       filter(covariate %in% covariates_model) %>%
#       mutate(time = .y)
#   }
# )

 
 
 ############## 
 
 
 library(dplyr)
 library(purrr)
 library(tidyr)
 library(mice)
 library(survival)
 library(broom)
 library(ipw)
 library(cobalt)
 
 load("/Users/roxanecouturier/Desktop/Doctorat/Article 3 Version 2/graaph2023.RData")
 
 #### 1. SÉLECTION ET TYPAGE DES VARIABLES ####
 
 data <- graaph.2023 %>%
   select(
     greffe, delallo, Type.bmt, BRAS, randodt,
     dc, survie, delrec, rec, suivi, RANDOMIZATION_R1,
     age, SEXE, bmi, POIDSDG, TAILLEDG, CNSPRECOCE, CNSTARDIF, cNS, cNSTARDIF, PS, GBJ1,
     CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
     SdTUMORAL, MEDIASTINDG, TLP, TLPPRECOCE, TLPTARDIF,
     GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
     ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
     TRANSCRIT, STRATE, IKZF1, NO_IKZF1, RES_IKZF1,
     BLASTESJ1, CSR, IMMUNOJ1, IDBLASTPREPHASE,
     MRD1FAIT, MRD1LAB, MRD2FAIT, mrd1dt, mrd2dt, mrd3dt, mrd4dt,
     mrd1, mrd2, mrd3, mrd4
   )
 
 data$PS         <- factor(data$PS)
 data$SdTUMORAL  <- as.factor(data$SdTUMORAL)
 data$PLQDIAG    <- as.numeric(data$PLQDIAG)
 data$GBJ1       <- as.numeric(data$GBJ1)
 data$BLASTESJ1  <- as.numeric(data$BLASTESJ1)
 data$GBDIAG     <- as.numeric(data$GBDIAG)
 data$BLASTESDIAG<- as.numeric(data$BLASTESDIAG)
 
 data <- data %>%
   mutate(
     PS_regroup = case_when(
       PS %in% c(2, 3) ~ "2+",
       TRUE ~ as.character(PS)
     ) %>% factor()
   )
 
 #### 2. VARIABLES DÉRIVÉES ####
 
 data <- data %>%
   mutate(
     survie_month = survie / 30.4167,
     allo         = ifelse(is.na(delallo), 0, 1),
     trt          = ifelse(BRAS == 1, 1, 0),
     pfs          = case_when(rec == 1 ~ 1, dc == 1 ~ 1, TRUE ~ 0),
     pfs_time     = case_when(rec == 1 ~ delrec, dc == 1 ~ survie, TRUE ~ suivi) / 30.4167
   ) %>%
   mutate(
     randodt = as.Date(randodt),
     mrd1dt  = as.Date(mrd1dt),
     mrd2dt  = as.Date(mrd2dt),
     mrd3dt  = as.Date(mrd3dt),
     mrd4dt  = as.Date(mrd4dt)
   ) %>%
   mutate(
     mrd1_month = as.numeric(mrd1dt - randodt) / 30.4167,
     mrd2_month = as.numeric(mrd2dt - randodt) / 30.4167,
     mrd3_month = as.numeric(mrd3dt - randodt) / 30.4167,
     mrd4_month = as.numeric(mrd4dt - randodt) / 30.4167
   ) %>%
   mutate(
     age_cat40 = factor(ifelse(age < 40, "< 40", ">= 40"), levels = c("< 40", ">= 40")),
     age_cat50 = factor(ifelse(age < 50, "< 50", ">= 50"), levels = c("< 50", ">= 50"))
   ) %>%
   mutate(patient_id = row_number())
 
 #### 3. IMPUTATION (mice) ####
 
 imp_data <- data %>%
   select(
     patient_id,
     cNS, BLASTESJ1, GBJ1, PS, PS_regroup, SdTUMORAL,
     age, age_cat40, age_cat50, SEXE, bmi, POIDSDG, TAILLEDG,
     CNSPRECOCE, CNSTARDIF, cNSTARDIF,
     CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
     MEDIASTINDG, TLP, TLPPRECOCE, TLPTARDIF,
     GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
     ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
     TRANSCRIT, STRATE, IKZF1, RES_IKZF1, NO_IKZF1,
     CSR, IMMUNOJ1, IDBLASTPREPHASE,
     mrd1, pfs, pfs_time, allo, mrd1_month
   )
 
 pred_matrix <- make.predictorMatrix(imp_data)
 pred_matrix[, "patient_id"] <- 0
 pred_matrix["patient_id", ] <- 0
 
 meth <- make.method(imp_data)
 meth["patient_id"] <- ""
 
 imp <- mice(imp_data, m = 20, method = meth, predictorMatrix = pred_matrix, seed = 123)
 # plot(imp)  # à décommenter pour vérifier la convergence
 
 imputed_vars <- complete(imp, 1) %>%
   select(patient_id, cNS, BLASTESJ1, GBJ1, PS, PS_regroup, SdTUMORAL,
          MEDIASTINDG, GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
          TLPTARDIF, RES_IKZF1, TLPPRECOCE, TLP, NO_IKZF1)
 
 data <- data %>%
   select(-cNS, -BLASTESJ1, -GBJ1, -PS, -PS_regroup, -SdTUMORAL,
          -MEDIASTINDG, -GBDIAG, -PNNDIAG, -PLQDIAG, -BLASTESDIAG,
          -TLPTARDIF, -RES_IKZF1, -TLPPRECOCE, -TLP, -NO_IKZF1) %>%
   left_join(imputed_vars, by = "patient_id")
 
 #### 4. FORMAT LONG POUR LA PFS : mrd change de valeur à mrd1_month ####
 
 make_long_mrd1_pfs <- function(patient_id, pfs_time, pfs, mrd1_month, mrd1) {
   if (is.na(mrd1_month)) {
     tibble(
       patient_id = patient_id, tstart = 0, tstop = pfs_time,
       mrd = 0, event = ifelse(pfs == 1, 1, 0)
     )
   } else {
     bind_rows(
       tibble(patient_id = patient_id, tstart = 0, tstop = mrd1_month, mrd = 0, event = 0),
       tibble(patient_id = patient_id, tstart = mrd1_month, tstop = pfs_time,
              mrd = mrd1, event = ifelse(pfs == 1, 1, 0))
     )
   }
 }
 
 data_long <- pmap_dfr(
   list(patient_id = data$patient_id, pfs_time = data$pfs_time,
        pfs = data$pfs, mrd1_month = data$mrd1_month, mrd1 = data$mrd1),
   make_long_mrd1_pfs
 )
 
 data_long <- data_long %>%
   left_join(
     data %>%
       select(
         patient_id, age, SEXE, bmi, POIDSDG, TAILLEDG, delallo, BRAS,
         PS, PS_regroup, cNS, age_cat40, age_cat50, GBJ1, mrd1_month,
         CENTRE, CENTREBMT, FRATRIE, FRATRIEBMT,
         SdTUMORAL, MEDIASTINDG, CNSPRECOCE, TLP, TLPPRECOCE, TLPTARDIF,
         GBDIAG, PNNDIAG, PLQDIAG, BLASTESDIAG,
         ALBUMINEDIAG, LDHDIAG, ASATDIAG, ALATDIAG, CREATDIAG, BILITDIAG,
         TRANSCRIT, STRATE, IKZF1, NO_IKZF1, RES_IKZF1,
         BLASTESJ1, CSR, IMMUNOJ1, IDBLASTPREPHASE
       ),
     by = "patient_id"
   ) %>%
   mutate(
     allo_event = ifelse(!is.na(delallo) & delallo <= tstop, 1, 0),
     event1     = ifelse(!is.na(delallo) & delallo <= tstop, 0, event),
     tstop      = ifelse(!is.na(delallo) & delallo <= tstop, delallo, tstop),
     trt        = ifelse(BRAS == 1, 1, 0)
   )
 
 #### 5. SCREENING UNIVARIÉ : effet des covariables sur la PFS et l'allogreffe ####
 
 covariates <- c(
   "age", "SEXE", "bmi", "POIDSDG", "TAILLEDG", "trt",
   "CENTRE", "CENTREBMT", "FRATRIE", "FRATRIEBMT",
   "SdTUMORAL", "MEDIASTINDG", "CNSPRECOCE",
   "TLP", "TLPPRECOCE", "TLPTARDIF",
   "GBDIAG", "PNNDIAG", "PLQDIAG", "BLASTESDIAG",
   "ALBUMINEDIAG", "LDHDIAG", "ASATDIAG", "ALATDIAG", "PS_regroup",
   "CREATDIAG", "BILITDIAG",
   "TRANSCRIT", "STRATE", "IKZF1", "NO_IKZF1", "RES_IKZF1",
   "BLASTESJ1", "CSR", "IMMUNOJ1", "IDBLASTPREPHASE", "GBJ1", "cNS",
   "age_cat40", "age_cat50", "PS", "mrd"
 )
 
 run_univariate <- function(outcome, covariates, data) {
   map_dfr(covariates, function(cov) {
     f <- as.formula(paste0("Surv(tstart, tstop, ", outcome, ") ~ ", cov, " + cluster(patient_id)"))
     fit <- tryCatch(coxph(f, data = data), error = function(e) NULL)
     if (is.null(fit)) {
       return(tibble(covariate = cov, term = NA, HR = NA, lower = NA, upper = NA,
                     p.value = NA, note = "erreur/non convergent"))
     }
     tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
       mutate(covariate = cov, .before = 1) %>%
       select(covariate, term, HR = estimate, lower = conf.low, upper = conf.high, p.value)
   })
 }
 
 # Ici "event" = "event1" (événement de PFS censuré à l'allogreffe)
 res_pfs  <- run_univariate("event1", covariates, data_long)
 res_allo <- run_univariate("allo_event", covariates, data_long)
 
 sig_pfs  <- res_pfs  %>% filter(!is.na(p.value), p.value < 0.05) %>% arrange(p.value)
 sig_allo <- res_allo %>% filter(!is.na(p.value), p.value < 0.05) %>% arrange(p.value)
 
 print(sig_pfs,  n = 40)
 print(sig_allo, n = 40)
 
 #### 6. ÉCLATEMENT DU JEU DE DONNÉES À TOUS LES TEMPS (pour IPCW) ####
 
 data_long_df <- as.data.frame(data_long) %>% arrange(patient_id, tstart)
 
 # marqueur de censure administrative (avant éclatement)
 data_long_df <- data_long_df %>%
   arrange(patient_id, tstart) %>%
   group_by(patient_id) %>%
   mutate(is_last = row_number() == n(),
          censored = ifelse(is_last & event1 == 0, 1, 0)) %>%
   ungroup() %>%
   as.data.frame()
 
 times <- sort(unique(data_long_df$tstop))
 
 data.long <- survSplit(
   data_long_df, cut = times, start = "tstart", end = "tstop", event = "event1"
 )
 
 data.long.c <- survSplit(
   data_long_df, cut = times, start = "tstart", end = "tstop", event = "allo_event"
 )
 data.long$censored <- data.long.c$allo_event
 
 #### 7. POIDS IPCW (exposure = allogreffe, traitée comme censure informative) ####
 
 w_full <- ipwtm(
   exposure    = censored,
   family      = "survival",
   numerator   = ~ 1,
   # denominator = ~  SdTUMORAL +GBJ1 +BLASTESDIAG +trt +mrd +MEDIASTINDG +MEDIASTINDG:mrd+ MEDIASTINDG:trt+ SdTUMORAL:mrd + BLASTESDIAG:mrd + GBJ1:mrd+
   #   SdTUMORAL:trt + GBJ1:trt + BLASTESDIAG:trt+mrd:trt+ GBJ1:BLASTESDIAG + GBJ1:MEDIASTINDG,
   denominator = ~  SdTUMORAL  +BLASTESJ1 +trt +  NO_IKZF1 + mrd + mrd:trt, #MEDIASTINDG:mrd+ MEDIASTINDG:trt+ SdTUMORAL:mrd + BLASTESJ1:mrd + 
   #SdTUMORAL:trt + BLASTESJ1:trt+mrd:trt + BLASTESJ1:SdTUMORAL1,
   id          = patient_id,
   tstart      = tstart,
   timevar     = tstop,
   type        = "cens",
   data        = data.long
 )$ipw.weights
 
 summary(w_full)
 sd(w_full)
 
 quantile(w_full, c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)
 
 data.long$ipw_full  <- w_full
 data.long$ipw_trunc <- pmin(w_full, quantile(w_full, 0.99, na.rm = TRUE))
 
 #### 9. VÉRIFICATION DE L'ÉQUILIBRE : bal.tab / love.plot ####
 # 
 # data_last <- data.long %>%
 #   arrange(patient_id, tstart) %>%
 #   group_by(patient_id) %>%
 #   slice_tail(n = 1) %>%
 #   ungroup()
 # 
 # table(data_last$allo_event, useNA = "always")
 # summary(data_last$ipw_full)
 # 
 # covariates_model <- c("SdTUMORAL", "BLASTESJ1", "trt", "mrd", "MEDIASTINDG","NO_IKZF1")
 # 
 # bt <- bal.tab(
 #   data_last[, covariates_model],
 #   treat   = data_last$allo_event,
 #   weights = data_last$ipw_full,
 #   method  = "weighting",
 #   un      = TRUE
 # )
 # bt
 # 
 # love.plot(
 #   bt,
 #   stat         = "mean.diffs",
 #   threshold    = 0.1,
 #   abs          = TRUE,
 #   var.order    = "unadjusted",
 #   colors       = c("#D55E00", "#0072B2"),
 #   shapes       = c(16, 17),
 #   sample.names = c("Before weighting", "After IPCW (PFS, split at all times)"),
 #   title        = "Covariate balance across alloSCT status (PFS)",
 #   line         = FALSE
 # )
 
 #### 10. MODÈLES DE COX PONDÉRÉS PAR IPW (effet du bras sur la PFS) ####
 
 cox_unweighted <- coxph(Surv(tstart, tstop, event1) ~ trt + cluster(patient_id),
                         data = data.long)
 summary(cox_unweighted)
 
 cox_ipw <- coxph(Surv(tstart, tstop, event1) ~ trt + cluster(patient_id),
                  data = data.long, weights = ipw_full)
 summary(cox_ipw)
 
 cox_ipw_trunc <- coxph(Surv(tstart, tstop, event1) ~ trt + cluster(patient_id),
                        data = data.long, weights = ipw_trunc)
 summary(cox_ipw_trunc)
 