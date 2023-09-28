#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","xgboost",
              "here","foreach","doParallel")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

remotes::install_github("yqzhong7/AIPW")
library(AIPW)
library(SuperLearner)
library(tmle3)
library(sl3)

a <- rio::import(here("data","numom_processed_2023_03_28.csv"))

miss_func <- function(x){
  mean(is.na(x))
}

apply(a, 2, miss_func)

val <- unique(a$multivit[!is.na(a$multivit)])
mode_multivit <- val[which.max(tabulate(match(a$multivit, val)))]
a$multivit[is.na(a$multivit)] <- mode_multivit

apply(a, 2, miss_func)

outcome <- as.matrix(a$pree_acog)

exposure <- as.matrix(a$fv_totdens_2_5)

covariate_list <-  c("momage",              "bmiprepreg",          "smokecigspre",       
                     "gravidity",           "pa_totmetwk_v1",     
                     "puqe_tot",            "realm_tot",            
                     "momaccult",           "epds_tot_v1",         "stress_tot_v1",
                     "anx_tot",             "walk_nat",            "adi_nat",
                     "povperc",             "momracehisp2",        "momracehisp3",
                     "momracehisp4",        "momracehisp5",        "momracehisp6",
                     "momracehisp7",        "momracehisp8",
                     "smokerpre1",          "multivit",
                     "marital2",            "marital4",            "marital5",
                     "insurance2",          "insurance3",          "momeduc2",
                     "momeduc3",            "momeduc4",            "momeduc5",
                     "momeduc6",            "artcat2",             "artcat3",
                     "alc_bingeprepreg1",   "prediab1",            "prehtn1",
                     "pregplanned1",        "sleepsat1",
                     "hei2015_total_nofv")

covariates <- a %>% select(all_of(covariate_list))

number_seeds <- as.numeric(args[1])

submission_list <- as.numeric(args[2])

start.time <- Sys.time()

parallel::detectCores()

n.cores <- parallel::detectCores()

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK",
  outfile=""
)

parallel::clusterEvalQ(my.cluster, 
                       {library(AIPW)
                         library(SuperLearner)
                         library(tmle3)
                         library(sl3)
                         library(ranger)
                         library(xgboost)})

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()

par_res <- foreach(i = 1:number_seeds, .packages = c("AIPW",
                                                     "SuperLearner",
                                                     "xgboost")) %dopar% {
  
  set.seed(i)
  
  SL.xgboost_grid <- create.Learner("SL.xgboost",
                                    tune=list(minobspernode = c(10,50)),
                                    params = list(verbose=1))
  
  sl_lib <- SL.xgboost_grid$names
  
  AIPW_SL <- AIPW$new(Y = outcome,
                      A = exposure,
                      W = covariates, 
                      Q.SL.library = sl_lib,
                      g.SL.library = sl_lib,
                      k_split = 5,
                      verbose=T,
                      save.sl.fit = T
  )$
    fit()$
    summary(g.bound = 0.025)
  
  aipw_sl_Q.res <- apply(
    data.frame(
      do.call(rbind,lapply(1:5, function(x) AIPW_SL$libs$Q.fit[[x]]$coef)), 
      seed = i
    ), 2, mean
  )
  aipw_sl_Q.res <- data.frame(t(aipw_sl_Q.res))
  aipw_sl_Q.res$model <- "Q.mod"
  aipw_sl_g.res <- apply(
    data.frame(
      do.call(rbind,lapply(1:5, function(x) AIPW_SL$libs$g.fit[[x]]$coef)), 
      seed = i
    ), 2, mean
  )
  aipw_sl_g.res <- data.frame(t(aipw_sl_g.res))
  aipw_sl_g.res$model <- "g.mod"
  
  aipw_res <- data.frame(algorithm = "AIPW",
                         t(AIPW_SL$result[3,1:4]), 
                         seed = i)
  names(aipw_res) <- c("Algorithm", 
                       "Estimate",
                       "SE", "LCL", "UCL",
                       "Seed")
  res <- list(
              aipw_res,
    rbind(aipw_sl_Q.res,
          aipw_sl_g.res)
  )
  
  return(res)
  
} # end of foreach loop

parallel::stopCluster(cl = my.cluster)

seed_res <- do.call(rbind, lapply(1:number_seeds, function(x) par_res[[x]][[1]]))
seed_res$submission <- submission_list

write_csv(seed_res, file = here("data", "seed_results.csv"))

aipw_res <- do.call(rbind, lapply(1:number_seeds, function(x) par_res[[x]][[2]]))
aipw_res$submission <- submission_list

write_csv(aipw_res, file = here("data", "aipw_sl_coef_results.csv"))

end.time <- Sys.time()

duration_time <- end.time - start.time

dur_dat <- data.frame(run_time = duration_time,
                      number_seeds = number_seeds,
                      submission = submission_list)

dur_dat

write_csv(dur_dat, 
          file = here("data","run_time.csv"), 
          append = T)
