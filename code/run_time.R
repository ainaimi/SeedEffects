#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

userLib <-  "~/R/R_LIBS_USER"
.libPaths(userLib)

packages <- c("data.table","tidyverse","skimr","here", "rio")

for (package in packages) {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package, repos='http://lib.stat.cmu.edu/R/CRAN')
  }
}

for (package in packages) {
  library(package, character.only=T)
}

a <- rio::import(here("data","run_time.csv"))

names(a) <- c("time", "seeds", "run")

pred_mod <- lm(time ~ seeds, data = a)

pred_dat <- as.matrix((predict(pred_mod, newdata = data.frame(seeds = c(5000, 10000, 20000, 50000)))/60)/24)

row.names(pred_dat) <- c("5000 seeds",
                            "10,000 seeds",
                            "20,000 seeds",
                            "50,000 seeds")

colnames(pred_dat) <- "days"

print(pred_dat)