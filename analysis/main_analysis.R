# Analysis on synthetic data -------------------------------------------
library(tidyverse)
library(BaBooN)
library(mice)
library(miceadds)
library(optmatch)
library(MatchIt)
library(WeightIt)
library(rms)
library(misaem)
library(twang)
library(CBPS)
library(cobalt)
library(naniar)
library(tictoc)

# The aim of this script is to run the main analysis on the synthetic
# dataset

# 1) Define inputs for the analysis ------------------------------------
seed <- 87963214
n_iter <- 20L
n_imp <- 20L
source(here::here("analysis", "misc_functions.R"))

# 2) Prepare the data --------------------------------------------------
df <- read_rds(here::here("data", "synthetic_data.rds"))

treat <- "treat"
covs <- df %>%
   dplyr::select(-treat) %>%
   names()

trt_frm <- as.formula(
   paste(treat, paste0(covs, collapse = " + "), sep = " ~ ")
)

# Evaluate the balance of the raw dataset ------------------------------
bal_pre <- bal.tab(
   trt_frm,
   data = df,
   continous = "std",
   binary = "std",
   s.d.denom = "pooled",
   abs = TRUE
)

love.plot(bal_pre)

# 3) Complete case analysis --------------------------------------------
db_cc <- df %>%
   na.omit() %>%
   mutate_if(is.factor, droplevels)

# Run the design stage analysis ----------------------------------------
cc_res <- design_stage(treatment = treat, covs = covs, data = db_cc)

cc_res$balance <- cc_res$balance %>%
   add_column(approach = "cc", .before = 1)

# 4) Missing indicator -------------------------------------------------
mind_res <- design_stage_mind(treatment = treat, covs = covs, data = df)

mind_res$balance <- mind_res$balance %>%
   add_column(approach = "mind", .before = 1)

# 5) SAEM --------------------------------------------------------------
saem_res <- design_stage_saem(
   treatment = treat, covs = covs, data = df, seed = seed
)

saem_res$balance <- saem_res$balance %>%
   add_column(approach = "saem", .before = 1)

# 6) GBM ---------------------------------------------------------------
gbm_res <- design_stage_gbm(
   treatment = treat, covs = covs, data = df, seed = seed
)

gbm_res$balance <- gbm_res$balance %>%
   add_column(approach = "gbm", .before = 1)

save(
   cc_res, mind_res, saem_res, gbm_res,
   file = here::here("analysis", "res_no_imp.rda")
)

rm(cc_res, mind_res, saem_res, gbm_res)

# 7) FCS with MICE -----------------------------------------------------
# 7A) Imputation stage -------------------------------------------------
imp_fcs <- mi_mice(
   data = df,
   seed = seed,
   n_imp = n_imp,
   n_iter = n_iter,
   cont_method = "pmm",
   binary_method = "logreg",
   cat_method = "polyreg"
)

# 7B) Design stage -----------------------------------------------------
list_imp_df <- imp_fcs$imp_list %>%
   set_names(x = ., nm = glue::glue("imp_{1:length(.)}"))

# plan(multisession)


tic()

mice_balance <- purrr::imap(
   .x = list_imp_df,
   ~ {

      dd <- design_stage(treatment = treat, covs = covs, data = .x)

      dd$balance <- dd$balance %>%
         tibble::add_column(imp = .y, .before = 1)

      message(glue::glue("{.y} is finished"))

      dd

   }
)

toc()

mice_res <- list(
   "imputation" = imp_fcs,
   "design_stage" = mice_balance
)

save(
   mice_res,
   file = here::here("analysis", "res_mice.rda")
)

rm(mice_res, imp_fcs, list_imp_df, mice_balance)

# 8) Bayesian Bootstrap ------------------------------------------------
# 8A) Imputation stage -------------------------------------------------
imp_bb <- mi_baboon(
   data = df,
   seed = seed,
   n_imp = n_imp,
   n_iter = n_iter
)

# 8B) Design stage -----------------------------------------------------
list_imp_df <- imp_bb$imp_list %>%
   set_names(x = ., nm = glue::glue("imp_{1:length(.)}"))

tic()

bb_balance <- purrr::imap(
   .x = list_imp_df,
   ~ {

      dd <- design_stage(treatment = treat, covs = covs, data = .x)

      dd$balance <- dd$balance %>%
         tibble::add_column(imp = .y, .before = 1)

      dd

   }
)

toc()

bb_res <- list(
   "imputation" = imp_bb,
   "design_stage" = bb_balance
)

save(
   bb_res,
   file = here::here("analysis", "res_bb.rda")
)

rm(bb_res, imp_bb, list_imp_df, bb_balance)

# 9) aregImpute --------------------------------------------------------
# 9A) Imputation stage -------------------------------------------------
imp_aregimp <- mi_aregimpute(
   data = df,
   seed = seed,
   n_imp = n_imp,
   n_iter = n_iter,
   n_knots = 3L,
   pmmtype = 2
)

# 9B) Design stage -----------------------------------------------------
list_imp_df <- imp_aregimp$imp_list %>%
   set_names(x = ., nm = glue::glue("imp_{1:length(.)}")) %>%
   purrr::map(
      .x = .,
      ~ .x %>%
         dplyr::mutate(
            dplyr::across(c(V2:V8), ~ as.double(.))
         ) %>%
         dplyr::mutate(
            dplyr::across(c(V9:V15), ~ factor(.))
         )
   )

aregimp_balance <- vector(mode = "list", length = length(list_imp_df))

for(i in 1:length(list_imp_df)) {

   aregimp_balance[[i]] <- design_stage(
      treatment = treat, covs = covs, data = list_imp_df[[i]]
   )

}

aregimp_res <- list(
   "imputation" = imp_aregimp,
   "design_stage" = aregimp_balance
)

save(
   aregimp_res,
   file = here::here("analysis", "res_aregimp.rda")
)

rm(aregimp_res, imp_aregimp, list_imp_df, aregimp_balance)
