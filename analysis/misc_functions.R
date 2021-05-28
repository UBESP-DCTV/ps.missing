# Misc functions -------------------------------------------------------
# PS with logistic -----------------------------------------------------
ps_logit <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   # Fit the treatment assignment model
   ps_fit <- stats::glm(trt_frm, data = data, family = binomial("logit"))

   # Return the estimated PS, both on logit and probability scales
   list(
      "prob" = ps_fit$fitted.values,
      "logit" = ps_fit$linear.predictors
   )

}

# PS with GBM ----------------------------------------------------------
ps_gbm <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   # Fit the treatment assignment model
   ps_fit <- twang::ps(
      trt_frm,
      data = as.data.frame(data),
      verbose = TRUE,
      stop.method = "es.mean"
   )

   # Return the estimated PS, both on logit and probability scales
   list(
      "prob" = ps_fit$ps[, 1],
      "logit" = qlogis(ps_fit$ps[, 1])
   )

}

# PS with CBPS ---------------------------------------------------------
ps_cbps <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   # Fit the treatment assignment model
   ps_fit <- CBPS::CBPS(trt_frm, data = as.data.frame(data), ATT = 0)

   # Return the estimated PS, both on logit and probability scales
   list(
      "prob" = ps_fit$fitted.values,
      "logit" = ps_fit$linear.predictor[, 1]
   )

}

# Estimate PS ----------------------------------------------------------
ps_est <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   # Create the formula of the treatment model
   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   # Estimate PS
   ps_logistic <- ps_logit(treatment, covs, data)
   ps_mlt <- ps_gbm(treatment, covs, data)
   ps_cs <- ps_cbps(treatment, covs, data)

   # Retrieve PS names
   ps_nm <- names(ps_logistic$prob)

   list(
      "logit" = ps_logistic,
      "gbm" = purrr::map(.x = ps_mlt, ~ purrr::set_names(.x, ps_nm)),
      "cbps" = purrr::map(.x = ps_cs, ~ purrr::set_names(.x, ps_nm))
   )

}

# NN matching ----------------------------------------------------------
match_nn <- function(treatment, covs, data, ps_list) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)
   assertive::assert_is_list(ps_list)

   # Create the formula of the treatment model
   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   ps_logit_scale <- list(
      "logit" = ps_list$logit$logit,
      "gbm" = ps_list$gbm$logit,
      "cbps" = ps_list$cbps$logit
   )

   # 2) Perform matching -----------------------
   purrr::map(
      .x = ps_logit_scale,
      ~ {

         # Perform matching
         mtch_obj <- MatchIt::matchit(
            trt_frm,
            data = data,
            method = "nearest",
            distance = .x,
            ratio = 1,
            caliper = 0.2 * sd(.x)
         )

         m_data <- MatchIt::match.data(mtch_obj)

         list("mtch_obj" = mtch_obj, "data" = data, "m_data" = m_data)

      }
   )

}

# Balance with NN matching ---------------------------------------------
bal_match_nn <- function(mtch_obj) {

   if(class(mtch_obj$mtch_obj) != "matchit") {
      usethis::ui_stop(
         "All the objects contained in `mtch_objs` must belong to the
         `matchit` class"
      )

   }

   # Create the balance table
   bal_tab <- cobalt::bal.tab(
      mtch_obj$mtch_obj,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs", "ovl.coefficients"),
      abs = TRUE,
      quick = FALSE
   )

   # Single SMDs
   single_smd <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::filter(variable != "distance")

   # Store the balancing statistics
   bal_df_single <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::filter(variable == "distance")

   # Average SMDs and Overlapping coefficients (pre and post)
   asmd_pre <- bal_df_single$Diff.Un
   ovl_pre <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::filter(variable == "distance") %>%
      .[["OVL.Un"]]

   asmd_post <- bal_df_single$Diff.Adj
   ovl_post <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::filter(variable == "distance") %>%
      .[["OVL.Adj"]]

   # Pre-C statistics
   tr <- mtch_obj$data$treat
   ps <- stats::plogis(mtch_obj$mtch_obj$distance)

   c_stat_pre <- pROC::auc(pROC::roc(tr, ps))

   # Post-C statistics
   m_data <- mtch_obj$m_data
   tr <- m_data$treat
   ps <- stats::plogis(m_data$distance)

   c_stat_post <- pROC::auc(pROC::roc(tr, ps))

   # Prop retained subjects
   prop <- (mtch_obj$mtch_obj$nn[2, 1] + mtch_obj$mtch_obj$nn[2, 2])/
      (mtch_obj$mtch_obj$nn[1, 1] + mtch_obj$mtch_obj$nn[1, 2])

   # Store everything into a tibble
   tibble::tibble(
      ps_method = "match_nn",
      single_smd = list(single_smd),
      smd_pre = asmd_pre,
      smd_post = asmd_post,
      ovl_pre = ovl_pre,
      ovl_post = ovl_post,
      c_stat_pre = as.double(c_stat_pre),
      c_stat_post = as.double(c_stat_post),
      prop_retained = prop
   )

}

# Full-matching --------------------------------------------------------
full_match <- function(treatment, covs, data, ps_list) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)
   assertive::assert_is_list(ps_list)

   # Create the formula of the treatment model
   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   ps_logit_scale <- list(
      "logit" = ps_list$logit$logit,
      "gbm" = ps_list$gbm$logit,
      "cbps" = ps_list$cbps$logit
   )

   # 2) Perform matching -----------------------
   purrr::map(
      .x = ps_logit_scale,
      ~ {

         # Perform matching
         mtch_obj <- MatchIt::matchit(
            trt_frm,
            data = data,
            method = "full",
            distance = .x,
            caliper = 0.2 * sd(.x)
         )

         m_data <- MatchIt::match.data(mtch_obj)

         list("mtch_obj" = mtch_obj, "data" = data, "m_data" = m_data)

      }
   )

}

# Balance with Full-matching -------------------------------------------
bal_full_match <- function(fm_obj) {


   if(class(fm_obj$mtch_obj)[1] != "matchit") {
      usethis::ui_stop(
         "`fm_obj` must belong be of `matchit` class"
      )
   }

   # Balance statistics with cobalt
   # Create the balance table
   bal_tab <- cobalt::bal.tab(
      fm_obj$mtch_obj,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs", "ovl.coefficients"),
      abs = TRUE,
      quick = FALSE
   )

   # Single SMDs
   single_smd <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::filter(variable != "distance")

   # Store the balancing statistics
   bal_df_single <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::filter(variable == "distance")

   # Average SMDs and Overlapping coefficients (pre and post)
   asmd_pre <- bal_df_single$Diff.Un
   ovl_pre <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::filter(variable == "distance") %>%
      .[["OVL.Un"]]

   asmd_post <- bal_df_single$Diff.Adj
   ovl_post <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::filter(variable == "distance") %>%
      .[["OVL.Adj"]]

   # Pre-C statistics
   tr <- fm_obj$data$treat
   ps <- stats::plogis(fm_obj$mtch_obj$distance)

   c_stat_pre <- pROC::auc(pROC::roc(tr, ps))

   # Post-C statistics
   m_data <- fm_obj$m_data

   c_stats_post <- purrr::map_dbl(
      .x = sort(unique(m_data$subclass)),
      ~ {

         dd <- m_data %>%
            dplyr::filter(subclass == .x)

         ps <- stats::plogis(dd$distance)
         tr <- dd$treat
         wts <- dd$weights

         wt_roc <- WeightedROC::WeightedROC(
            stats::plogis(dd$distance),
            dd$treat,
            dd$weights
         )

         WeightedROC::WeightedAUC(wt_roc)
      }
   )

   wts <- table(m_data$subclass)

   c_stat_post <- weighted.mean(x = c_stats_post, w = wts)

   # Prop retained subjects
   prop <- (fm_obj$mtch_obj$nn[2, 1] + fm_obj$mtch_obj$nn[2, 2])/
      (fm_obj$mtch_obj$nn[1, 1] + fm_obj$mtch_obj$nn[1, 2])

   # Store everything into a tibble
   tibble::tibble(
      ps_method = "full_match",
      single_smd = list(single_smd),
      smd_pre = asmd_pre,
      smd_post = asmd_post,
      ovl_pre = ovl_pre,
      ovl_post = ovl_post,
      c_stat_pre = as.double(c_stat_pre),
      c_stat_post = as.double(c_stat_post),
      prop_retained = prop
   )

}

# PS-IPTW --------------------------------------------------------------
ps_iptw <- function(treatment, covs, data, ps_list) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)
   assertive::assert_is_list(ps_list)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   ps_vector <- list(
      "logit" = ps_list$logit$prob,
      "gbm" = ps_list$gbm$prob,
      "cbps" = ps_list$cbps$prob
   )

   # 2) Perform weighting
   purrr::map(
      .x = ps_vector,
      ~ {

         # Perform weighting
         WeightIt::weightit(
            trt_frm,
            data = data,
            estimand = "ATE",
            ps = .x
         )
      }
   )

}

# Balance with PS-IPTW -------------------------------------------------
bal_ps_iptw <- function(ps_iptw_obj) {

   if(class(ps_iptw_obj) != "weightit") {
      usethis::ui_stop("`ps_iptw_objs` must be of `weightit` class")
   }

   # Balance statistics with cobalt
   # Create the balance table
   bal_tab <- cobalt::bal.tab(
      ps_iptw_obj,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs", "ovl.coefficients"),
      abs = TRUE,
      quick = FALSE
   )

   # Single SMDs
   single_smd <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::slice(-1)

   # Store the balancing statistics
   bal_df_single <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::slice(1)

   # SMD of PS and Overlapping coefficients (pre and post)
   asmd_pre <- bal_df_single$Diff.Un
   ovl_pre <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::slice(1) %>%
      .[["OVL.Un"]]

   asmd_post <- bal_df_single$Diff.Adj
   ovl_post <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      dplyr::slice(1) %>%
      .[["OVL.Adj"]]

   # Pre-C statistics
   tr <- ps_iptw_obj$treat
   ps_est <- ps_iptw_obj$ps

   c_stat_pre <- pROC::auc(pROC::roc(tr, ps_est))

   # Post-C statistics
   wts <- ps_iptw_obj$weights

   roc_post <- WeightedROC::WeightedROC(ps_est, tr, wts)
   c_stat_post <- WeightedROC::WeightedAUC(roc_post)

   # Prop retained subjects
   ss <- summary(ps_iptw_obj)$effective.sample.size
   prop <- (ss[2, 1] + ss[2, 2])/(ss[1, 1] + ss[1, 2])

   # Store everything into a tibble
   tibble::tibble(
      ps_method = "ps_iptw",
      single_smd = list(single_smd),
      smd_pre = asmd_pre,
      smd_post = asmd_post,
      ovl_pre = ovl_pre,
      ovl_post = ovl_post,
      c_stat_pre = as.double(c_stat_pre),
      c_stat_post = as.double(c_stat_post),
      prop_retained = prop
   )
}

# Design stage analysis ------------------------------------------------
design_stage <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   # 1) PS estimation --------------------------------------------------
   ps_list <- ps_est(treatment, covs, data)

   # 2) NN matching ----------------------------------------------------
   nn_match <- match_nn(treatment, covs, data, ps_list)

   nn_match_bal <- purrr::map2(
      .x = nn_match, .y = names(nn_match),
      ~ bal_match_nn(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   ) %>%
      dplyr::bind_rows()

   # 3) Full-matching --------------------------------------------------
   fm <- full_match(treatment, covs, data, ps_list)

   fm_bal <- purrr::map2_dfr(
      .x = fm, .y = names(fm),
      ~ bal_full_match(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   )

   # 4) PS-IPTW --------------------------------------------------------
   ipw <- ps_iptw(treatment, covs, data, ps_list)

   ipw_bal <- purrr::map2_dfr(
      .x = ipw, .y = names(ipw),
      ~ bal_ps_iptw(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   )

   # 5) Store balance results ------------------------------------------
   bal_res <- dplyr::bind_rows(nn_match_bal, fm_bal, ipw_bal)

   ps_res <- purrr::map(
      .x = ps_list,
      ~ tibble::tibble(
         ps = .x$prob,
         treat = data$treat
      )
   )

   ll <- list("balance" = bal_res, "ps" = ps_res)

}

# Preprocess data with missing indicator method ------------------------
preproc_mind <- function(data) {

   assertive::assert_is_data.frame(data)

   d_na <- data %>%
      dplyr::select(-treat) %>%
      dplyr::mutate(
         dplyr::across(.cols = dplyr::everything(), ~ is.na(.))
      ) %>%
      dplyr::rename_all(~ glue::glue("{.}_NA"))

   d_imp <- data %>%
      dplyr::mutate(
         dplyr::across(
            where(is.double),
            ~ dplyr::if_else(is.na(.), mean(., na.rm = TRUE), .)
         )
      ) %>%
      dplyr::mutate(
         dplyr::across(
            where(is.factor),
            ~ dplyr::if_else(
               is.na(.),
               mean(as.double(.) - 1, na.rm = TRUE),
               as.double(.) - 1
            )
         )
      )

   dd <- dplyr::bind_cols(d_imp, d_na) %>%
      # Missing indicator variables as 0 and 1
      dplyr::mutate_at(
         dplyr::vars(dplyr::ends_with("_NA")),
         ~ dplyr::if_else(., 1, 0)
      )

   # Retrieve the covariates for PS model and treatment indicator
   tr <- "treat"
   covs <- dd %>%
      dplyr::select(-treat) %>%
      names()

   list(
      "data" = dd,
      "treatment" = tr,
      "covs" = covs
   )
}

# Design stage missing indicator ---------------------------------------
design_stage_mind <- function(treatment, covs, data) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)

   # 1) Preproc data ---------------------------------------------------
   df_preproc <- preproc_mind(data)

   # 2) Estimate PS ----------------------------------------------------
   df_mind <- df_preproc$data

   ps_list <- ps_est(df_preproc$treatment, df_preproc$covs, df_mind)

   # 3) NN matching ----------------------------------------------------
   nn_match <- match_nn(
      df_preproc$treatment, df_preproc$covs, df_mind, ps_list
   )

   nn_match_bal <-

   nn_match_bal <- purrr::imap_dfr(
      .x = nn_match,
      ~ bal_match_nn(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   )

   # 4) Full-matching --------------------------------------------------
   fm <- full_match(
      df_preproc$treatment, df_preproc$covs, df_mind, ps_list
   )

   fm_bal <- purrr::imap_dfr(
      .x = fm,
      ~ bal_full_match(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   )

   # 5) PS-IPTW --------------------------------------------------------
   ipw <- ps_iptw(
      df_preproc$treatment, df_preproc$covs, df_mind, ps_list
   )

   ipw_bal <- purrr::imap_dfr(
      .x = ipw,
      ~ bal_ps_iptw(.x) %>%
         tibble::add_column(ps_est = .y, .after = 1)
   )

   # 6) Store balance results ------------------------------------------
   bal_res <- dplyr::bind_rows(nn_match_bal, fm_bal, ipw_bal)

   ps_res <- purrr::map(
      .x = ps_list,
      ~ tibble::tibble(
         ps = .x$prob,
         treat = data$treat
      )
   )

   list("balance" = bal_res, "ps" = ps_res)

}

# Design stage saem ----------------------------------------------------
design_stage_saem <- function(treatment, covs, data, seed) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)
   assertive::assert_is_a_number(seed)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   df_misaem <- data %>%
      dplyr::mutate(dplyr::across(where(is.factor), ~ as.double(.) - 1))

   ps_fit <- misaem::miss.glm(trt_frm, data = df_misaem)
   ps <- misaem:::predict.miss.glm(
      ps_fit, newdata = df_misaem %>% dplyr::select(-.data[[treatment]])
   )

   df_ps <- data %>%
      dplyr::mutate(
         ps = as.double(ps) %>%
            purrr::set_names(., nm = seq_len(length(ps))),
         ps_logit = stats::qlogis(ps)
      ) %>%
      dplyr::select(treat, ps, ps_logit)

   # 2) NN matching
   # 2A) Matching
   mtch_obj <- MatchIt::matchit(
      treat ~ ps,
      data = df_ps,
      method = "nearest",
      distance = df_ps$ps_logit,
      ratio = 1,
      caliper = 0.2 * sd(df_ps$ps_logit)
   )

   m_data <- MatchIt::match.data(mtch_obj)

   mtch_obj_nn <- list(
      "mtch_obj" = mtch_obj, "data" = df_ps, "m_data" = m_data
   )

   # 2B) Balance
   nn_match_bal <- bal_match_nn(mtch_obj_nn)

   # 3) Full-matching
   # 3A) Matching
   mtch_obj <- MatchIt::matchit(
      treat ~ ps,
      data = df_ps,
      method = "full",
      distance = df_ps$ps_logit,
      caliper = 0.2 * sd(df_ps$ps_logit)
   )

   m_data <- MatchIt::match.data(mtch_obj)

   mtch_obj_fm <- list(
      "mtch_obj" = mtch_obj, "data" = df_ps, "m_data" = m_data
   )

   # 3B) Balance
   bal_match_fm <- bal_full_match(mtch_obj_fm)

   # 4) PS-IPTW
   wt_obj <- WeightIt::weightit(
      treat ~ ps,
      data = df_ps,
      estimand = "ATE",
      ps = df_ps$ps
   )

   # 4B) Balance
   bal_tab_single_pre <- cobalt::bal.tab(
      trt_frm,
      data = data,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs"),
      abs = TRUE
   )

   bal_tab_single <- cobalt::bal.tab(
      trt_frm,
      data = data,
      weights = wt_obj$weights,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs"),
      abs = TRUE
   )

   bal_tab <- cobalt::bal.tab(
      treat ~ ps,
      data = df_ps,
      weights = wt_obj$weights,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs", "ovl.coefficients"),
      abs = TRUE,
      quick = FALSE
   )

   # Single SMDs
   single_smd_pre <- bal_tab_single_pre$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un)

   single_smd <- bal_tab_single$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::mutate(Diff.Un = single_smd_pre$Diff.Un)

   # Store the balancing statistics
   bal_df_single <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      filter(variable != "prop.score")

   # Average SMDs and Overlapping coefficients (pre and post)
   asmd_pre <- bal_df_single$Diff.Un
   ovl_pre <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      filter(variable == "ps") %>%
      .[["OVL.Un"]]

   asmd_post <- bal_df_single$Diff.Adj
   ovl_post <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      filter(variable == "ps") %>%
      .[["OVL.Adj"]]

   # Pre-C statistics
   tr <- data$treat
   ps_est <- wt_obj$ps

   c_stat_pre <- pROC::auc(pROC::roc(tr, ps_est))

   # Post-C statistics
   wts <- wt_obj$weights

   roc_post <- WeightedROC::WeightedROC(ps_est, tr, wts)
   c_stat_post <- WeightedROC::WeightedAUC(roc_post)

   # Prop retained subjects
   ss <- summary(wt_obj)$effective.sample.size
   prop <- (ss[2, 1] + ss[2, 2])/(ss[1, 1] + ss[1, 2])

   # Store everything into a tibble
   bal_ps_ipw <- tibble::tibble(
      ps_method = "ps_iptw",
      single_smd = list(single_smd),
      smd_pre = asmd_pre,
      smd_post = asmd_post,
      ovl_pre = ovl_pre,
      ovl_post = ovl_post,
      c_stat_pre = as.double(c_stat_pre),
      c_stat_post = as.double(c_stat_post),
      prop_retained = prop
   )

   # 5) Final balance table
   bal_res <- dplyr::bind_rows(nn_match_bal, bal_match_fm, bal_ps_ipw) %>%
      tibble::add_column(ps_est = "saem", .after = 1)

   list("balance" = bal_res, "ps" = df_ps)

}

# Design stage GBM -----------------------------------------------------
design_stage_gbm <- function(treatment, covs, data, seed) {

   assertive::assert_is_character(treatment)
   assertive::assert_is_character(covs)
   assertive::assert_is_data.frame(data)
   assertive::assert_is_a_number(seed)

   trt_frm <- stats::as.formula(
      paste(treatment, paste0(covs, collapse = " + "), sep = " ~ ")
   )

   # 1) Estimate PS with GBM, including missing data
   set.seed(seed)
   wt_obj <- WeightIt::weightit(
      trt_frm,
      data = data,
      method = "gbm",
      estimand = "ATE",
      missing = "surr"
   )

   df_ps <- data %>%
      dplyr::mutate(
         ps = wt_obj$ps,
         ps_logit = stats::qlogis(wt_obj$ps)
      ) %>%
      dplyr::select(treat, ps, ps_logit) %>%
      dplyr::mutate(
         dplyr::across(
            c(ps, ps_logit),
            ~ purrr::set_names(., seq_len(length(.)))
         )
      )

   # 2) NN matching
   # 2A) Matching
   mtch_obj <- MatchIt::matchit(
      treat ~ ps,
      data = df_ps,
      method = "nearest",
      distance = df_ps$ps_logit,
      ratio = 1,
      caliper = 0.2 * sd(df_ps$ps_logit)
   )

   m_data <- MatchIt::match.data(mtch_obj)

   mtch_obj_nn <- list(
      "mtch_obj" = mtch_obj, "data" = df_ps, "m_data" = m_data
   )

   # 2B) Balance
   nn_match_bal <- bal_match_nn(mtch_obj_nn)

   # 3) Full-matching
   # 3A) Matching
   mtch_obj <- MatchIt::matchit(
      treat ~ ps,
      data = df_ps,
      method = "full",
      distance = df_ps$ps_logit,
      caliper = 0.2 * sd(df_ps$ps_logit)
   )

   m_data <- MatchIt::match.data(mtch_obj)

   mtch_obj_fm <- list(
      "mtch_obj" = mtch_obj, "data" = df_ps, "m_data" = m_data
   )

   # 3B) Balance
   bal_match_fm <- bal_full_match(mtch_obj_fm)

   # 4) PS-IPTW
   # 4B) Balance
   bal_tab_single_pre <- cobalt::bal.tab(
      trt_frm,
      data = data,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs"),
      abs = TRUE
   )

   bal_tab_single <- cobalt::bal.tab(
      trt_frm,
      data = data,
      weights = wt_obj$weights,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs"),
      abs = TRUE
   )

   bal_tab <- cobalt::bal.tab(
      treat ~ ps,
      data = df_ps,
      weights = wt_obj$weights,
      continuous = "std",
      binary = "std",
      s.d.denom = "pooled",
      stats = c("mean.diffs", "ovl.coefficients"),
      abs = TRUE,
      quick = FALSE
   )

   # Single SMDs
   single_smd_pre <- bal_tab_single_pre$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un)

   single_smd <- bal_tab_single$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      dplyr::mutate(Diff.Un = single_smd_pre$Diff.Un)

   # Store the balancing statistics
   bal_df_single <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, Diff.Un, Diff.Adj) %>%
      filter(variable != "prop.score")

   # Average SMDs and Overlapping coefficients (pre and post)
   asmd_pre <- bal_df_single$Diff.Un
   ovl_pre <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      filter(variable == "ps") %>%
      .[["OVL.Un"]]

   asmd_post <- bal_df_single$Diff.Adj
   ovl_post <- bal_tab$Balance %>%
      tibble::as_tibble(rownames = "variable") %>%
      dplyr::select(variable, OVL.Un, OVL.Adj) %>%
      filter(variable == "ps") %>%
      .[["OVL.Adj"]]

   # Pre-C statistics
   tr <- data$treat
   ps_est <- wt_obj$ps

   c_stat_pre <- pROC::auc(pROC::roc(tr, ps_est))

   # Post-C statistics
   wts <- wt_obj$weights

   roc_post <- WeightedROC::WeightedROC(ps_est, tr, wts)
   c_stat_post <- WeightedROC::WeightedAUC(roc_post)

   # Prop retained subjects
   ss <- summary(wt_obj)$effective.sample.size
   prop <- (ss[2, 1] + ss[2, 2])/(ss[1, 1] + ss[1, 2])

   # Store everything into a tibble
   bal_ps_ipw <- tibble::tibble(
      ps_method = "ps_iptw",
      single_smd = list(single_smd),
      smd_pre = asmd_pre,
      smd_post = asmd_post,
      ovl_pre = ovl_pre,
      ovl_post = ovl_post,
      c_stat_pre = as.double(c_stat_pre),
      c_stat_post = as.double(c_stat_post),
      prop_retained = prop
   )

   # 5) Final balance table
   bal_res <- dplyr::bind_rows(nn_match_bal, bal_match_fm, bal_ps_ipw) %>%
      tibble::add_column(ps_est = "gbm", .after = 1)

   list("balance" = bal_res, "ps" = df_ps)

}

# MI with mice ---------------------------------------------------------
mi_mice <- function(
   data,
   seed,
   n_imp = 5L,
   n_iter = 20L,
   cont_method = "pmm",
   binary_method = "logreg",
   cat_method = "polyreg"
) {

   assertive::assert_is_data.frame(data)
   assertive::assert_is_a_number(seed)
   assertive::assert_is_a_number(n_imp)
   assertive::assert_is_a_number(n_iter)
   assertive::assert_is_character(cont_method)
   assertive::assert_is_character(binary_method)
   assertive::assert_is_character(cat_method)

   # Remove the treatment variable
   df_to_imp <- data

   # Quick initialization of mice algorithm
   init_mice <- mice::mice(df_to_imp, m = 1, maxit = 0)

   # Define the imputation methods: pmm for continuous variables,
   # logistic regression for binary and multinomial for categorical.
   # They are the defaul methods for mice

   # Run the mice algorithm
   set.seed(seed)
   mice_imp <- mice::mice(
      df_to_imp,
      m = n_imp,
      maxit = n_iter,
      seed = seed
   )

   # Save the Rhat for convergence assessment
   rhat_conv <- miceadds::Rhat.mice(mice_imp)

   # Imputed dataframes into a list
   list(
      "original_data" = data,
      "mice_obj" = mice_imp,
      "convergence" = rhat_conv,
      "all" = mice::complete(mice_imp, action = "all", include = TRUE),
      "imp" = mice::complete(mice_imp, action = "all", include = FALSE),
      "imp_list" = as.list(
         mice::complete(mice_imp, action = "all", include = FALSE)
      )
   )
}

# Reconvert categorical variables to factors ---------------------------
reconvert_factor <- function(data, factor_info) {

   assertive::assert_is_data.frame(data)
   assertive::assert_is_list(factor_info)

   numerical_data <- data %>%
      dplyr::select(!dplyr::matches(names(factor_info)))

   categorical_data <- purrr::imap_dfc(
      .x = factor_info,
      ~ {

         zz <- factor(data[, .y])
         levels(zz) <- .x

         zz

      }

   )

   dplyr::bind_cols(numerical_data, categorical_data) %>%
      dplyr::select(dplyr::matches(names(data)))

}

# List of imputed dfs with BaBooN as mids object -----------------------
list_to_mids_baboon <- function(original_data, imp_data) {

   assertive::assert_is_list(imp_data)
   assertive::assert_is_data.frame(original_data)

   imp_dfs <- purrr::map(
      .x = imp_data,
      ~ .x %>%
         tibble::add_column(.id = 1:nrow(.), .before = 1)
   ) %>%
      purrr::map2(
         .x = ., .y = 1:length(.),
         ~ .x %>% tibble::add_column(.imp = .y, .before = 1)
      ) %>%
      unname()

   all_dfs <- dplyr::bind_rows(
      original_data %>%
         dplyr::mutate(.imp = 0, .id = 1:nrow(.)),
      imp_dfs
   ) %>%
      mice::as.mids()

   list(
      "all" = all_dfs,
      "imp_list" = imp_dfs
   )
}

# MI with baboon -------------------------------------------------------
mi_baboon <- function(
   data,
   seed,
   n_imp = 5L,
   n_iter = 20L
) {

   assertive::assert_is_data.frame(data)
   assertive::assert_is_a_number(seed)
   assertive::assert_is_a_number(n_imp)
   assertive::assert_is_a_number(n_iter)

   df_to_imp <- data

   # Store into a separate object the information on the factors,
   # i.e. variables names and levels
   factor_info <- purrr::map(
      .x = names(dplyr::select_if(df_to_imp, is.factor)),
      ~ levels(df_to_imp[[.x]])
   ) %>%
      purrr::set_names(names(dplyr::select_if(df_to_imp, is.factor)))

   # Preprocess the dataframe. Transform factors into doubles
   df_pre_proc <- df_to_imp %>%
      dplyr::mutate_if(is.factor, ~ as.double(.))

   # Run the imputation algorithm
   baboon_imp <- BBPMM(
      Data = df_pre_proc,
      M = n_imp,
      nIter = n_iter,
      verbose = TRUE,
      setSeed = seed
   )

   # Assess convergence
   diag_conv <- impdiagnosticconversion(baboon_imp, "mcmc")

   # Assign the original values of the factors to each imputed
   imp_data <- purrr::map(
      .x = baboon_imp$impdata,
      ~ .x %>%
         reconvert_factor(factor_info = factor_info)
   )

   # Transform the imputed data into appropriate objects
   mids_obj <- list_to_mids_baboon(df_to_imp, imp_data)

   list(
      "original_data" = data,
      "baboon_obj" = baboon_imp,
      "convergence" = diag_conv,
      "all" = mids_obj$all,
      "imp_list" = mids_obj$imp_list
   )

}

# List of imputed dfs with aregImpute as mids object -------------------
list_to_mids_aregimp <- function(original_data, imp_data) {

   if(!class(imp_data) == "aregImpute") {
      usethis::ui_stop("`imp_data` must of class `aregImpute`")
   }
   assertive::assert_is_data.frame(original_data)

   imp_dfs <- purrr::map(
      .x = 1:imp_data$n.impute,
      ~ Hmisc::impute.transcan(
         imp_data, data = original_data,
         list.out = TRUE, imputation = .x
      ) %>%
         cbind.data.frame() %>%
         tibble::add_column(.id = 1:nrow(.), .before = 1)
   ) %>%
      purrr::imap(
         .x = .,
         ~ .x %>% tibble::add_column(.imp = .y, .before = 1)
      )

   all_dfs <- dplyr::bind_rows(
      original_data %>%
         dplyr::mutate(.imp = 0, .id = 1:nrow(.)) %>%
         dplyr::mutate_if(is.factor, ~ as.character(.)),
      purrr::map(
         .x = imp_dfs,
         ~ .x %>%
            dplyr::mutate_if(is.factor, ~ as.character(.))
      )
   ) %>%
      dplyr::mutate_if(is.character, ~ factor(.)) %>%
      mice::as.mids(.)

   list(
      "all" = all_dfs,
      "imp_list" = imp_dfs
   )
}


# MI with aregImpute ---------------------------------------------------
mi_aregimpute <- function(
   data,
   seed,
   n_imp = 5L,
   n_iter = 20L,
   n_knots = 3L,
   pmmtype = 2
) {

   assertive::assert_is_data.frame(data)
   assertive::assert_is_a_number(seed)
   assertive::assert_is_a_number(n_imp)
   assertive::assert_is_a_number(n_knots)
   assertive::assert_is_a_number(pmmtype)

   df_to_imp <- data

   # Prepare formula for aregImpute
   aregimp_frm <- as.formula(
      paste("~", paste0(names(df_to_imp), collapse = " + "), sep = " ")
   )

   # Run the aregImpute algorithm
   set.seed(seed)
   areg_imp <- Hmisc::aregImpute(
      aregimp_frm,
      data = df_to_imp,
      n.impute = n_imp,
      burnin = n_iter,
      nk = n_knots,
      type = "pmm",
      pmmtype = pmmtype
   )

   # Transform the imputed data into appropriate objects
   mids_obj <- list_to_mids_aregimp(df_to_imp, areg_imp)

   list(
      "original_data" = data,
      "aregimp_obj" = areg_imp,
      "all" = mids_obj$all,
      "imp_list" = mids_obj$imp_list
   )

}
