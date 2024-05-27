pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom,
  parallel
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

set.seed(123)

n <- 100

simulation_function <- function(index){
  
  sim_res <- data.frame(index, 
                        mean1 = mean(rnorm(n)),
                        mean2 = mean(rbinom(n, size = 1, prob = .25)),
                        mean3 = mean(runif(n)))
  
  if(index == 1){
    write_csv(sim_res, here("data", "write_csv_parallel_test.csv"))
  } else{
    write_csv(sim_res, here("data", "write_csv_parallel_test.csv"), append = T)
  }

}

# simulation_function(index = 1)

simulation_res <- mclapply(1:1e6, 
                           function(x) simulation_function(index = x), 
                           mc.cores = 10)

obj1 <- read_csv(here("data", "write_csv_parallel_test_v1.csv"))
obj2 <- read_csv(here("data", "write_csv_parallel_test.csv"))

head(obj1)
head(obj2)

sum(is.na(obj1))
sum(is.na(obj2))