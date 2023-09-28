packages <- c("data.table","tidyverse","skimr",
              "here","foreach","doParallel",
              "DoubleML","mlr3","mlr3learners","mlr3pipelines")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

remotes::install_github("yqzhong7/AIPW")
remotes::install_github("ehkennedy/npcausal")
library(AIPW)
library(npcausal)
library(SuperLearner)
library(ggplot2)
library(tmle)
library(tmle3)
library(sl3)

a <- rio::import(here("data","numom_processed_2023_03_28.csv"))

VIM::aggr(a)

miss_func <- function(x){
  mean(is.na(x))
}

apply(a, 2, miss_func)

val <- unique(a$multivit[!is.na(a$multivit)])
mode_multivit <- val[which.max(tabulate(match(a$multivit, val)))]
a$multivit[is.na(a$multivit)] <- mode_multivit

apply(a, 2, miss_func)

index <- sample(1:nrow(a), 2000)

a <- a[index, ]

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

number_seeds <- 12

start.time <- Sys.time()

parallel::detectCores()

n.cores <- 6

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

foreach::getDoParRegistered()

par_res <- foreach(i = 1:number_seeds) %dopar% {
  
  create.SL.glmnet <- function(alpha = c(0.25, 0.50, 0.75)) {
    for(mm in seq(length(alpha))){
      eval(parse(text = paste('SL.glmnet.', alpha[mm], '<- function(..., alpha = ', alpha[mm], ') SL.glmnet(..., alpha = alpha)', sep = '')), envir = .GlobalEnv)
    }
    invisible(TRUE)
  }
  
  set.seed(i)
  
  create.SL.glmnet(alpha = 0)
  
  AIPW_SL <- AIPW$new(Y = outcome,
                      A = exposure,
                      W = covariates, 
                      Q.SL.library = c("SL.glmnet.0", "SL.glmnet","SL.glm", "SL.ranger", "SL.xgboost"),
                      g.SL.library = c("SL.glmnet.0", "SL.glmnet","SL.glm", "SL.ranger", "SL.xgboost"),
                      k_split = 5,
                      verbose=FALSE)$
    fit()$
    summary(g.bound = 0.025)
  
  aipw_res <- data.frame(algorithm = "AIPW",
                         t(AIPW_SL$result[3,1:4]), 
                         seed = i)
  names(aipw_res) <- c("Algorithm", 
                       "Estimate",
                       "SE", "LCL", "UCL",
                       "Seed")
  
  ate.res <- ate(y=outcome,
                 a=exposure,
                 x=covariates,
                 nsplits = 5,
                 sl.lib = c("SL.glmnet.0", "SL.glmnet","SL.glm", "SL.ranger", "SL.xgboost")
                 )
  
  npcausal_res <- data.frame(algorithm = "npcausal",
                             ate.res$res[3,2:5],
                             seed = i)
  names(npcausal_res) <- c("Algorithm", 
                           "Estimate",
                           "SE", "LCL", "UCL",
                           "Seed")
  
  
  sl_ <- make_learner(Stack, unlist(list(make_learner(Lrnr_glmnet, alpha = 1),
                                         make_learner(Lrnr_glmnet, alpha = 0),
                                         make_learner(Lrnr_xgboost),
                                         make_learner(Lrnr_glm),
                                         make_learner(Lrnr_ranger)),
                                    recursive = TRUE))
    
  learner <- Lrnr_sl$new(learners = sl_,
                           metalearner = Lrnr_nnls$new(convex=T))
  learner_list <- list(Y = learner,
                       A = learner)
  
  ######################################################################
  
  # PREPARE THE THINGS WE WANT TO FEED IN TO TMLE3
  ate_spec <- tmle_ATE(treatment_level = 1, 
                       control_level = 0)
  
  nodes_ <- list(W = names(covariates), 
                 A = "fv_totdens_2_5", 
                 Y = "pree_acog")
  
  # RUN TMLE3 
  tmle_res <- tmle3(ate_spec, 
                    a, 
                    nodes_, 
                    learner_list)
  tmle3_res <- data.frame(algorithm = "tmle3",
                          tmle_res$summary[,4:7],
                          seed = i)
  names(tmle3_res) <- c("Algorithm", 
                        "Estimate",
                        "SE", "LCL", "UCL",
                        "Seed")
  
  # Gruber TMLE
  sl.lib <- c("SL.glmnet.0", "SL.glmnet","SL.glm", "SL.ranger", "SL.xgboost")
  tmleG <- tmle(
    Y = outcome,
    A = exposure,
    W = covariates,
    Q.SL.library = sl.lib,
    g.SL.library = sl.lib,
    v = 5
  )
  tmleG_res <- data.frame(algorithm = "tmle",
                          t(c(tmleG$estimates$ATE$psi, tmleG$estimates$ATE$var.psi, tmleG$estimates$ATE$CI)),
                          seed = i)
  names(tmleG_res) <- c("Algorithm", 
                        "Estimate",
                        "SE", "LCL", "UCL",
                        "Seed")
  
  # DML
  # dml_data = DoubleMLData$new(a,
  #                             y_col = "aqoc",
  #                             d_cols = "pg",
  #                             x_cols = names(covariates))
  # 
  # lrn_1 = lrn("regr.cv_glmnet", alpha = 1, s = "lambda.min")
  # lrn_10 = po("learner_cv", lrn_1$clone())
  # lrn_10$id = "glmnet1_cv"
  # 
  # lrn_2 = lrn("regr.cv_glmnet", alpha = 0, s = "lambda.min")
  # lrn_20 = po("learner_cv", lrn_2$clone())
  # lrn_20$id = "glmnet0_cv"
  # 
  # lrn_3 = lrn("regr.lm")
  # lrn_30 = po("learner_cv", lrn_3$clone())
  # lrn_30$id = "lm_cv"
  # 
  # lrn_4 = lrn("regr.ranger")
  # lrn_40 = po("learner_cv", lrn_4$clone())
  # lrn_40$id = "ranger_cv"
  # 
  # lrn_5 = lrn("regr.xgboost")
  # lrn_50 = po("learner_cv", lrn_5$clone())
  # lrn_50$id = "xgboost_cv"
  # 
  # level_0 = gunion(list(lrn_10, lrn_20, lrn_30, lrn_40, lrn_50, po("nop")))
  # combined = level_0 %>>% po("featureunion", 6)
  # stack = combined %>>% po("learner", lrn_1$clone())
  # stack$plot(html = FALSE)
  # 
  # stacklrn = as_learner(stack)
  # 
  # dml_res = DoubleMLPLR$new(dml_data,
  #                           ml_l=stacklrn,
  #                           ml_m=stacklrn)
  # dml_res$fit()
  # 
  # dml2_res <- data.frame(algorithm = "dml2",
  #                         t(c(dml_res$all_coef,dml_res$se,
  #                             dml_res$all_coef-1.96*dml_res$se,
  #                             dml_res$all_coef+1.96*dml_res$se)),
  #                         seed = i)
  # names(dml2_res) <- c("Algorithm", 
  #                       "Estimate",
  #                       "SE", "LCL", "UCL",
  #                       "Seed")
  
  
  res <- rbind(aipw_res, 
               npcausal_res,
               tmle3_res,
               tmleG_res)#,
               #dml2_res)
  
  return(res)
  
}

parallel::stopCluster(cl = my.cluster)

seed_res <- do.call(rbind, par_res)

seed_res

write_csv(seed_res, file = here("data", "seed_results.csv"))

hist(seed_res$Seed)

seed_res$Seed

end.time <- Sys.time()

duration_time <- end.time - start.time

dur_dat <- data.frame(run_time = duration_time,
                      number_seeds = number_seeds)

dur_dat

write_csv(dur_dat, 
          file = here("data","run_time.csv"), 
          append = T)


# time-analysis
# 
dur <- read_csv(here("data","run_time.csv"))
plot(dur)
mod <- lm(run_time ~ number_seeds, data = dur)
predict(mod, newdata = data.frame(
  num_sims = 10000
))
# 
# 
# 


seed_res %>% 
  group_by(Algorithm) %>% 
  summarise(tibble(min = min(Estimate), 
                   mean = mean(Estimate),
                   median = median(Estimate),
                   max = max(Estimate),
                   sd = sd(Estimate)))

ggplot(seed_res) + 
  geom_histogram(aes(Estimate)) +
  facet_wrap(~Algorithm)
