pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom,
  rslurm
)

lib_vers           <- 4
sampling_frac      <- 0.1

a <- read_csv(
  here("data", paste0("seed_effects_data_",sampling_frac,".csv"))
)

seed_func <- function(seed_value, library_version, outcome =            a$pree_acog, 
                      exposure =           a$fv_totdens_2_5, 
                      covariates =         a[,-c(1,2)], 
                      covariates_augment = data.frame(a[,-c(1,2)], 
                                                      exposure = a$fv_totdens_2_5)){
  
  set.seed(seed_value)
  
  print(paste("Now running seed number ", seed_value))
  
  sl_lib_list <- list(
    sl_lib <- 
      c("SL.mean", "SL.glm"),
    c("SL.mean", "SL.glm", "SL.glmnet"),
    c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger"),
    c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.xgboost")
  )
  
  sl_lib <- unlist(sl_lib_list[library_version])
  
  print(paste("This is the SL Library: ", sl_lib))
  
  print(paste("Setting up CV folds: seed number ", seed_value))
  
  num.folds <- 10
  
  folds <- sort(seq(nrow(covariates)) %% num.folds) + 1
  fold_dat <- tibble(id = 1:nrow(covariates),folds)
  fold_index <- split(fold_dat$id,fold_dat$folds)
  
  print(paste("Outcome Model: seed number ", seed_value))
  
  fit_mu <- CV.SuperLearner(Y = outcome,
                            X = covariates_augment, 
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index, shuffle=F),
                            control = SuperLearner.control(saveCVFitLibrary = TRUE),
                            parallel = "seq",
                            verbose = T)
  
  print(paste("PS Model: seed number ", seed_value))
  
  fit_pi <- CV.SuperLearner(Y = exposure,
                            X = covariates,
                            method = "method.NNLS", 
                            family = binomial,
                            SL.library = sl_lib,
                            cvControl = list(V = num.folds, validRows = fold_index, shuffle=F),
                            control = SuperLearner.control(saveCVFitLibrary = FALSE),
                            parallel = "seq",
                            verbose = T)
  
  print(paste("Constructing Cross-Fit Predictions: seed number ", seed_value))
  
  pscore <- as.matrix(fit_pi$SL.predict)
  
  mu_hat <- as.matrix(fit_mu$SL.predict)
  
  mu_hat1 <- NULL
  for(i in 1:num.folds){
    mu_hat1 <- rbind(mu_hat1, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_augment[fold_index[[i]],], exposure = 1), 
                             onlySL=T)$pred)
  }
  
  mu_hat0 <- NULL
  for(i in 1:num.folds){
    mu_hat0 <- rbind(mu_hat0, 
                     predict(fit_mu$AllSL[[i]],
                             newdata = base::transform(
                               covariates_augment[fold_index[[i]],], exposure = 0), 
                             onlySL=T)$pred)
  }
  
  ## aipw
  aipw_func <- function(exposure, outcome, pscore, mu_hat, mu_hat0, mu_hat1){
    aipw_score <- ((2*exposure - 1)*(outcome - mu_hat))/((2*exposure - 1)*pscore + (1 - exposure)) + (mu_hat1 - mu_hat0)
    return(aipw_score)
  }
  
  print(paste("Constructing AIPW Scores: seed number ", seed_value))
  
  aipw_score <- aipw_func(exposure, 
                          outcome, 
                          pscore, 
                          mu_hat, 
                          mu_hat0, 
                          mu_hat1)
  colnames(aipw_score) <- NULL
  
  aipw_psi <- mean(aipw_score)
  
  aipw_se <- sd(aipw_score)/sqrt(nrow(covariates))
  
  aipw_ate <- c(aipw_psi, aipw_se)
  
  print(paste("Constructing Data Output Objects: seed number ", seed_value))
  
  aipw_sl_Q.res <- apply(
    data.frame(
      do.call(rbind,lapply(1:num.folds, function(x) fit_mu$AllSL[[x]]$coef)),
      seed = seed_value
    ), 2, mean
  )
  aipw_sl_Q.res <- data.frame(t(aipw_sl_Q.res))
  aipw_sl_Q.res$model <- "Q.mod"
  aipw_sl_g.res <- apply(
    data.frame(
      do.call(rbind,lapply(1:num.folds, function(x) fit_pi$AllSL[[x]]$coef)),
      seed = seed_value
    ), 2, mean
  )
  aipw_sl_g.res <- data.frame(t(aipw_sl_g.res))
  aipw_sl_g.res$model <- "g.mod"
  
  aipw_res <- data.frame(algorithm = "AIPW",
                         aipw_psi, aipw_se, 
                         LCL = aipw_psi - 1.96*aipw_se,
                         UCL = aipw_psi + 1.96*aipw_se,
                         p.value_z = 2*pnorm(-abs(aipw_psi/aipw_se)),
                         p.value_t = 2*pt(-abs(aipw_psi/aipw_se), df = nrow(covariates) - 1),
                         seed = seed_value)
  
  names(aipw_res) <- c("Algorithm",
                       "Estimate",
                       "SE", "LCL", "UCL",
                       "p.val_z", "p.val_t",
                       "Seed")  
  
  print(paste("Intermediate Write to File: seed number ", seed_value))
  
  write_csv(aipw_res, 
            here("data",paste0("aipw_res_",lib_vers,"_",array_num,".csv")),
            append = T)
  
  write_csv(aipw_sl_Q.res, 
            here("data",paste0("aipw_sl_Q.res_",lib_vers,"_",array_num,".csv")),
            append = T)
  
  write_csv(aipw_sl_g.res, 
            here("data",paste0("aipw_sl_g.res_",lib_vers,"_",array_num,".csv")),
            append = T)
  
  print(paste("Function Output: seed number ", seed_value))
  
  res <- list(
    aipw_res,
    aipw_sl_g.res,
    aipw_sl_Q.res
  )
  
  return(res)
  
}

#seed_value, library_version


seed_data <- read_csv(
  here("data", "random_seed_values.csv")
)

seed_data <- seed_data$random_seeds[1:5]


pars <- expand.grid(seed_data,
                    library_version = 1:4)
head(pars, 3)
tail(pars, 3)

names(pars) <- c("seed_value","library_version")

sjob <- slurm_apply(seed_func, pars, jobname = 'seed_effects',
                    nodes = 1, cpus_per_node = 10, submit = FALSE)