##########################################
##########################################
# Fa outs
##########################################
##########################################

##' @author Saulo F. S. Chaves 
##' Modified  by Gabriel Blasques 

fa.outs <- function(model, name.env, name.gen) {
  
  require(asreml)
  asreml.options(workspace = '6gb', pworkspace = '8gb')
  
  data <- as.data.frame(model$mf)
  
  temp   <- attr(attr(model$formulae$random, "factors"), "dimnames")[[2]]
  facall <- temp[grepl("\\bfa\\b", temp) & grepl(name.gen, temp)]
  
  genocall <- strsplit(facall, split = ":")[[1]][1]
  envcall  <- strsplit(facall, split = ":")[[1]][2]
  
  lev.env <- model$G.param[[facall]][[envcall]]$levels[
    -grep("Comp", model$G.param[[facall]][[envcall]]$levels)
  ]
  lev.gen <- model$G.param[[facall]][[genocall]]$levels
  
  num.env <- length(lev.env)
  num.gen <- length(lev.gen)
  
  summa <- summary(model)
  vc    <- summa$varcomp
  
  load <- vc[grep("fa\\d+", rownames(vc)), ]
  var  <- vc[grep("var",   rownames(vc)), 1]
  
  load$fa <- regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings <- do.call(
    cbind,
    lapply(split(load[, "component"], load$fa), c)
  )
  rownames(mat.loadings) <- lev.env
  
  svd.lambda <- svd(mat.loadings)
  D <- diag(svd.lambda$d^2, length(svd.lambda$d))
  
  if (sum(svd.lambda$u[, 1] < 0) / nrow(svd.lambda$u) > 0.5) {
    svd.lambda$u <- -svd.lambda$u
    svd.lambda$v <- -svd.lambda$v
  }
  
  mat.loadings.star <- svd.lambda$u
  colnames(mat.loadings.star) <- colnames(mat.loadings)
  rownames(mat.loadings.star) <- rownames(mat.loadings)
  
  coef_rand <- data.frame(summary(model, coef = TRUE)$coef.random)
  idx.scores <- grep(
    paste0("(?=.*", genocall, ")(?=.*Comp)"),
    rownames(coef_rand),
    perl = TRUE
  )
  scores <- coef_rand[idx.scores, , drop = FALSE]
  scores$fa <- regmatches(
    rownames(scores),
    regexpr("Comp\\d+", rownames(scores))
  )
  
  scor.vec <- unlist(
    lapply(split(scores[, "solution"], scores$fa), c)
  )
  
  scor.vec.star <- kronecker(
    sqrt(D) %*% t(svd.lambda$v),
    diag(num.gen)
  ) %*% scor.vec
  
  scor.mat.star <- matrix(
    scor.vec.star,
    nrow = num.gen,
    ncol = length(unique(scores$fa)),
    dimnames = list(lev.gen, unique(load$fa))
  )
  
  Gvcov <- mat.loadings.star %*% tcrossprod(D, mat.loadings.star) + diag(var)
  rownames(Gvcov) <- lev.env
  colnames(Gvcov) <- lev.env
  Gcor  <- cov2cor(Gvcov)
  
  expvar <- sum(diag(mat.loadings.star %*% tcrossprod(D, mat.loadings.star))) /
    sum(diag(Gvcov))
  
  expvar.j <- matrix(
    NA_real_,
    nrow = ncol(mat.loadings),
    ncol = num.env,
    dimnames = list(colnames(mat.loadings), lev.env)
  )
  
  for (i in 1:nrow(expvar.j)) {
    for (j in 1:ncol(expvar.j)) {
      expvar.j[i, j] <-
        100 * mat.loadings.star[j, i]^2 * diag(D)[i] /
        (sum(mat.loadings.star[j, ]^2 * diag(D)) + var[j])
    }
  }
  
  lambdacross <- mat.loadings.star %*% tcrossprod(D, mat.loadings.star)
  rownames(lambdacross) <- lev.env
  colnames(lambdacross) <- lev.env
  
  fullsv <- list()
  for (i in lev.env) {
    semivar <- matrix(
      NA_real_, nrow = num.env, ncol = 3,
      dimnames = list(lev.env, c("i", "j", "semivar"))
    )
    for (j in lev.env) {
      semivar[j, ] <- c(
        i, j,
        0.5 * (lambdacross[i, i] + lambdacross[j, j]) - lambdacross[i, j]
      )
    }
    fullsv[[i]] <- semivar
  }
  
  fullsv <- data.frame(do.call(rbind, fullsv))
  fullsv$semivar <- as.numeric(fullsv$semivar)
  fullsv <- reshape(
    data = fullsv,
    timevar = "i",
    idvar = "j",
    direction = "wide"
  )[ , -1]
  
  colnames(fullsv) <- sub("semivar\\.", "", colnames(fullsv))
  lambda_ASV <- 2 / (nrow(fullsv) * (nrow(fullsv) - 1)) *
    sum(fullsv[upper.tri(fullsv)])
  
  fullsv <- list()
  for (i in lev.env) {
    semivar <- matrix(
      NA_real_, nrow = num.env, ncol = 3,
      dimnames = list(lev.env, c("i", "j", "semivar"))
    )
    for (j in lev.env) {
      semivar[j, ] <- c(
        i, j,
        0.5 * (Gvcov[i, i] + Gvcov[j, j]) - Gvcov[i, j]
      )
    }
    fullsv[[i]] <- semivar
  }
  
  fullsv <- data.frame(do.call(rbind, fullsv))
  fullsv$semivar <- as.numeric(fullsv$semivar)
  fullsv <- reshape(
    data = fullsv,
    timevar = "i",
    idvar = "j",
    direction = "wide"
  )[ , -1]
  
  colnames(fullsv) <- sub("semivar\\.", "", colnames(fullsv))
  G_ASV <- 2 / (nrow(fullsv) * (nrow(fullsv) - 1)) *
    sum(fullsv[upper.tri(fullsv)])
  
  ASVR <- lambda_ASV / G_ASV
  
  modpred <- predict(
    model,
    classify   = paste(name.gen, name.env, sep = ":"),
    sed        = TRUE
  )
  
  blups <- modpred$pvals
  
  if ("std.error" %in% colnames(blups)) {
    blups <- blups[, !(colnames(blups) %in% "std.error"), drop = FALSE]
  }
  
  temp_blup <- data.frame(
    env = rep(lev.env, each = length(lev.gen)),
    gen = rep(lev.gen, times = length(lev.env)),
    marginal = kronecker(mat.loadings.star, diag(num.gen)) %*% scor.vec.star
  )
  
  colnames(temp_blup)[1:2] <- c(name.env, name.gen)
  
  blups <- merge(
    blups,
    temp_blup,
    by = c(name.env, name.gen)
  )
  colnames(blups)[colnames(blups) == "predicted.value"] <- "conditional"
  
  results <- list(
    rot.loads  = mat.loadings.star,
    Gvcov      = Gvcov,
    Gcor       = Gcor,
    diagnostics = c(
      expvar = round(expvar * 100, 3),
      ASVR   = round(ASVR * 100, 3),
      aic    = round(summa$aic, 3),
      bic    = round(summa$bic, 3),
      logl   = round(summa$loglik, 3)
    ),
    expvar_j   = expvar.j,
    rot.scores = scor.mat.star,
    blups      = blups
  )
  
  return(results)
}














##########################################
#########################################

##########################################
##########################################
library(tidyverse)
library(openxlsx)
library(cowplot)
library(forcats)
library(asreml)
library(dplyr)
library(purrr)
library(pls)


data <- readRDS("data1.RDS")



num.loc = nlevels(data$loc)
num.env = nlevels(data$env)
num.gen = nlevels(data$gen)
name.env = levels(data$env)
name.gen = levels(data$gen)




##########################################
##########################################
# Fa model Fit
##########################################
##########################################


mod4 = asreml(fixed = yield ~ repblk:env + env,
              random = ~ gen:fa(env, 4),
              residual = ~ dsum(~ id(units) | env),
              data = data,
              maxit = 60,
              workspace = '300Mb',
              na.action = na.method(x = "include", y = "include"))
mod4 = update(mod4)


fa.res <- fa.outs(mod4, name.env = "env", name.gen = "gen")
all <- fa.res$blups
fa4 = fa.res



##########################################
##########################################
# Env data agreg
##########################################
##########################################


data <- data

weather <- read.csv("Weather_Data.csv", sep = ";")

weather$Date <- as.Date(as.character(weather$Date), format = "%Y%m%d")
weather$Env  <- as.character(weather$Env)

data$env <- as.character(data$env)


dates_env <- data %>% filter(!is.na(planting_date), !is.na(harvest_date)) %>%
  count(env, planting_date, harvest_date, name = "n") %>% group_by(env) %>%
  arrange(desc(n), planting_date, harvest_date) %>% slice(1) %>% ungroup() %>%
  select(env, planting_date, harvest_date)


#Join weather data with crop cycle dates

weather_cycle <- weather %>% left_join(dates_env, by = c("Env" = "env")) %>%
  filter(!is.na(planting_date), !is.na(harvest_date)) %>%
  filter(Date >= planting_date, Date <= harvest_date) %>%
  mutate(dap = as.numeric(Date - planting_date) + 1,
         cycle_days = as.numeric(harvest_date - planting_date) + 1)


#Environmental variables

env_vars <- c("QV2M","T2MDEW","PS","RH2M","WS2M","GWETTOP","ALLSKY_SFC_SW_DWN",
              "ALLSKY_SFC_PAR_TOT","T2M_MAX","T2M_MIN","T2MWET","GWETROOT",
              "T2M","GWETPROF","ALLSKY_SFC_SW_DNI","PRECTOTCORR")


#Five-day window means

weather_5d <- weather_cycle %>%
  mutate(window = ceiling(dap / 5), 
         window = paste0("w", sprintf("%02d", window))) %>%
  group_by(Env, planting_date, harvest_date, cycle_days, window) %>%
  summarise(across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)), 
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Env, planting_date, harvest_date, cycle_days),
              names_from = window,
              values_from = all_of(env_vars),
              names_glue = "{.value}_{window}")

weather_5d <- as.data.frame(weather_5d)
rownames(weather_5d) <- weather_5d$Env

weather_5d <- weather_5d %>% 
  select(-Env, -planting_date, -harvest_date, -cycle_days)

weather_5d <- weather_5d[, colSums(is.na(weather_5d)) == 0, drop = FALSE]


#Fifteen-day window means

weather_15d <- weather_cycle %>%
  mutate(window = ceiling(dap / 15), 
         window = paste0("w", sprintf("%02d", window))) %>%
  group_by(Env, planting_date, harvest_date, cycle_days, window) %>%
  summarise(across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)), 
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Env, planting_date, harvest_date, cycle_days),
              names_from = window,
              values_from = all_of(env_vars),
              names_glue = "{.value}_{window}")

weather_15d <- as.data.frame(weather_15d)
rownames(weather_15d) <- weather_15d$Env

weather_15d <- weather_15d %>% 
  select(-Env, -planting_date, -harvest_date, -cycle_days)

weather_15d <- weather_15d[, colSums(is.na(weather_15d)) == 0, drop = FALSE]
 
#Thirty-day window means


weather_30d <- weather_cycle %>%
  mutate(window = ceiling(dap / 30),
         window = paste0("w", sprintf("%02d", window))) %>%
  group_by(Env, planting_date, harvest_date, cycle_days, window) %>%
  summarise(across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_wider(id_cols = c(Env, planting_date, harvest_date, cycle_days),
              names_from = window,
              values_from = all_of(env_vars),
              names_glue = "{.value}_{window}")

weather_30d <- as.data.frame(weather_30d)
rownames(weather_30d) <- weather_30d$Env

weather_30d <- weather_30d %>%
  select(-Env, -planting_date, -harvest_date, -cycle_days)

weather_30d <- weather_30d[, colSums(is.na(weather_30d)) == 0, drop = FALSE]


#Whole-cycle mean

weather_cycle_mean <- weather_cycle %>%
  group_by(Env, planting_date, harvest_date, cycle_days) %>%
  summarise(across(all_of(env_vars), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

weather_cycle_mean <- as.data.frame(weather_cycle_mean)
rownames(weather_cycle_mean) <- weather_cycle_mean$Env

weather_cycle_mean <- weather_cycle_mean %>%
  select(-Env, -planting_date, -harvest_date, -cycle_days)

weather_cycle_mean <- weather_cycle_mean[, colSums(is.na(weather_cycle_mean)) == 0, drop = FALSE]


#Export to Excel

wb <- createWorkbook()

addWorksheet(wb, "mean_5_days")
writeData(wb, "mean_5_days", weather_5d, rowNames = TRUE)

addWorksheet(wb, "mean_15_days")
writeData(wb, "mean_15_days", weather_15d, rowNames = TRUE)

addWorksheet(wb, "mean_30_days")
writeData(wb, "mean_30_days", weather_30d, rowNames = TRUE)

addWorksheet(wb, "mean_whole_cycle")
writeData(wb, "mean_whole_cycle", weather_cycle_mean, rowNames = TRUE)

setColWidths(wb, "mean_5_days", cols = 1:(ncol(weather_5d) + 1), widths = "auto")
setColWidths(wb, "mean_15_days", cols = 1:(ncol(weather_15d) + 1), widths = "auto")
setColWidths(wb, "mean_30_days", cols = 1:(ncol(weather_30d) + 1), widths = "auto")
setColWidths(wb, "mean_whole_cycle", cols = 1:(ncol(weather_cycle_mean) + 1), widths = "auto")

saveWorkbook(wb, "weather_aggregated_by_crop_cycle_wide.xlsx", overwrite = TRUE)



###############################
# CV LOO
###############################
EnvExp <- read.xlsx("features_G2F.xlsx", sheet = "mean_5_days", rowNames = TRUE)
EBLUPs_matrix <- matrix(NA, nrow = num.gen, ncol = num.env, dimnames = list(name.gen, name.env))


for (excluded_env in name.env) {
  
  train_data <- subset(data, env != excluded_env)
  
  EnvExp <- read.xlsx("features_G2F.xlsx", sheet = "mean_5_days", rowNames = TRUE)
  train_covamb <- EnvExp[rownames(EnvExp) != excluded_env, ]
  
  
  mod <- asreml(fixed = yield ~ repblk:env + env,
                random = ~ gen:fa(env, 4),
                residual = ~ dsum(~ id(units) | env),
                maxit = 60,
                workspace = '300Mb',
                data = train_data,
                na.action = na.method(x = "include", y = "include"))
  mod <- update(mod)
  
  fa.res <- fa.outs(mod, name.env = "env", name.gen = "gen")
  
  data.pls.lamb <- data.frame(
    lambda = I(fa.res$rot.loads),
    CovAmb = I(scale(as.matrix(train_covamb))))
  
  pls_model <- pls::plsr(lambda ~ CovAmb, validation = "none", data = data.pls.lamb, ncomp = 7)
  coef_pls <- coef(pls_model, intercept = TRUE)[,,1]
  
  EnvExp <- scale(EnvExp)
  matcovamb_pred <- EnvExp[rownames(EnvExp) == excluded_env, , drop = FALSE]
  
  pred_lambda <- matcovamb_pred %*% coef_pls[-1, ]
  pred_lambda <- cbind(
    fa1 = pred_lambda[,"fa1"] + coef_pls["(Intercept)", "fa1"],
    fa2 = pred_lambda[,"fa2"] + coef_pls["(Intercept)", "fa2"],
    fa3 = pred_lambda[,"fa3"] + coef_pls["(Intercept)", "fa3"],
    fa4 = pred_lambda[,"fa4"] + coef_pls["(Intercept)", "fa4"]
  )
  
  EBLUPs_pred <- as.matrix(pred_lambda[,"fa1"]) %*% fa.res$rot.scores[,"fa1"] +
    as.matrix(pred_lambda[,"fa2"]) %*% fa.res$rot.scores[,"fa2"] +
    as.matrix(pred_lambda[,"fa3"]) %*% fa.res$rot.scores[,"fa3"] +
    as.matrix(pred_lambda[,"fa4"]) %*% fa.res$rot.scores[,"fa4"] +
    mean(data$yield, na.rm = TRUE)
  
  colnames(EBLUPs_pred) <- levels(data$gen)
  EBLUPs_matrix[, excluded_env] <- EBLUPs_pred
}

EBLUPs_matrix_df <- as.data.frame(EBLUPs_matrix)


correlations <- list()
environments <- unique(all$env)

for (env in environments) {
  subset_all <- all %>% filter(env == !!env)
  if (env %in% colnames(EBLUPs_matrix_df)) {
    common_gen <- intersect(subset_all$gen, rownames(EBLUPs_matrix_df))
    if (length(common_gen) > 1) {
      true_values <- subset_all %>% filter(gen %in% common_gen) %>% arrange(gen) %>% pull(conditional)
      predicted_values <- EBLUPs_matrix_df[common_gen, env, drop = FALSE] %>% arrange(rownames(.)) %>% pull()
      cor_value <- cor(true_values, predicted_values, method = "spearman")
      correlations[[env]] <- cor_value
    }
  }
}

correlation_df <- data.frame(
  Environment = names(correlations),
  Spearman_Correlation = unlist(correlations)
)

meanCV4 <- mean(correlation_df$Spearman_Correlation, na.rm = TRUE)


rmse_list <- list()

for (env in environments) {
  subset_all <- all %>% filter(env == !!env)
  if (env %in% colnames(EBLUPs_matrix_df)) {
    common_gen <- intersect(subset_all$gen, rownames(EBLUPs_matrix_df))
    if (length(common_gen) > 1) {
      true_values <- subset_all %>% 
        filter(gen %in% common_gen) %>% 
        arrange(gen) %>% 
        pull(conditional)
      predicted_values <- EBLUPs_matrix_df[common_gen, env, drop = FALSE] %>% 
        arrange(rownames(.)) %>% 
        pull()
      rmse <- sqrt(mean((true_values - predicted_values)^2, na.rm = TRUE))
      rmse_list[[env]] <- rmse
    }
  }
}

rmse_df <- data.frame(
  Environment = names(rmse_list),
  RMSEP = unlist(rmse_list)
)

meanRMSEP <- mean(rmse_df$RMSEP, na.rm = TRUE)



###############################
# CV LOO MARS
###############################
EnvExp <- read.xlsx("MARS_features_G2F.xlsx", sheet = "mean_5_days_MARS", rowNames = TRUE)
EBLUPs_matrix <- matrix(NA, nrow = num.gen, ncol = num.env, dimnames = list(name.gen, name.env))


for (excluded_env in name.env) {
  
  train_data <- subset(data, env != excluded_env)
  
  EnvExp <- read.xlsx("MARS_features_G2F.xlsx", sheet = "mean_5_days_MARS", rowNames = TRUE)
  train_covamb <- EnvExp[rownames(EnvExp) != excluded_env, ]
  
  
  mod <- asreml(fixed = yield ~ repblk:env + env,
                random = ~ gen:fa(env, 4),
                residual = ~ dsum(~ id(units) | env),
                maxit = 60,
                workspace = '300Mb',
                data = train_data,
                na.action = na.method(x = "include", y = "include"))
  mod <- update(mod)
  
  fa.res <- fa.outs(mod, name.env = "env", name.gen = "gen")
  
  data.pls.lamb <- data.frame(
    lambda = I(fa.res$rot.loads),
    CovAmb = I(scale(as.matrix(train_covamb))))
  
  pls_model <- pls::plsr(lambda ~ CovAmb, validation = "none", data = data.pls.lamb, ncomp = 11)
  coef_pls <- coef(pls_model, intercept = TRUE)[,,1]
  
  EnvExp <- scale(EnvExp)
  matcovamb_pred <- EnvExp[rownames(EnvExp) == excluded_env, , drop = FALSE]
  
  pred_lambda <- matcovamb_pred %*% coef_pls[-1, ]
  pred_lambda <- cbind(
    fa1 = pred_lambda[,"fa1"] + coef_pls["(Intercept)", "fa1"],
    fa2 = pred_lambda[,"fa2"] + coef_pls["(Intercept)", "fa2"],
    fa3 = pred_lambda[,"fa3"] + coef_pls["(Intercept)", "fa3"],
    fa4 = pred_lambda[,"fa4"] + coef_pls["(Intercept)", "fa4"]
  )
  
  EBLUPs_pred <- as.matrix(pred_lambda[,"fa1"]) %*% fa.res$rot.scores[,"fa1"] +
    as.matrix(pred_lambda[,"fa2"]) %*% fa.res$rot.scores[,"fa2"] +
    as.matrix(pred_lambda[,"fa3"]) %*% fa.res$rot.scores[,"fa3"] +
    as.matrix(pred_lambda[,"fa4"]) %*% fa.res$rot.scores[,"fa4"] +
    mean(data$yield, na.rm = TRUE)
  
  colnames(EBLUPs_pred) <- levels(data$gen)
  EBLUPs_matrix[, excluded_env] <- EBLUPs_pred
}

EBLUPs_matrix_df <- as.data.frame(EBLUPs_matrix)


correlations <- list()
environments <- unique(all$env)

for (env in environments) {
  subset_all <- all %>% filter(env == !!env)
  if (env %in% colnames(EBLUPs_matrix_df)) {
    common_gen <- intersect(subset_all$gen, rownames(EBLUPs_matrix_df))
    if (length(common_gen) > 1) {
      true_values <- subset_all %>% filter(gen %in% common_gen) %>% arrange(gen) %>% pull(conditional)
      predicted_values <- EBLUPs_matrix_df[common_gen, env, drop = FALSE] %>% arrange(rownames(.)) %>% pull()
      cor_value <- cor(true_values, predicted_values, method = "spearman")
      correlations[[env]] <- cor_value
    }
  }
}

correlation_df <- data.frame(
  Environment = names(correlations),
  Spearman_Correlation = unlist(correlations)
)

meanCV4MARS <- mean(correlation_df$Spearman_Correlation, na.rm = TRUE)


rmse_list <- list()

for (env in environments) {
  subset_all <- all %>% filter(env == !!env)
  if (env %in% colnames(EBLUPs_matrix_df)) {
    common_gen <- intersect(subset_all$gen, rownames(EBLUPs_matrix_df))
    if (length(common_gen) > 1) {
      true_values <- subset_all %>% 
        filter(gen %in% common_gen) %>% 
        arrange(gen) %>% 
        pull(conditional)
      predicted_values <- EBLUPs_matrix_df[common_gen, env, drop = FALSE] %>% 
        arrange(rownames(.)) %>% 
        pull()
      rmse <- sqrt(mean((true_values - predicted_values)^2, na.rm = TRUE))
      rmse_list[[env]] <- rmse
    }
  }
}

rmse_df <- data.frame(
  Environment = names(rmse_list),
  RMSEP = unlist(rmse_list)
)

meanRMSEPMARS <- mean(rmse_df$RMSEP, na.rm = TRUE)




