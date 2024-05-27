pacman::p_load(
  rio,          
  here,         
  skimr,        
  tidyverse,     
  lmtest,
  sandwich,
  broom,
  xgboost,
  rmutil,
  SuperLearner,
  Matrix,
  mvtnorm
)

thm <- theme_classic() +
  theme(
    legend.position = "top",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
theme_set(thm)

# adapted from https://bit.ly/3V96pff
# set the seed for reproducibility
set.seed(123)

# generate the observed data
n=10000
sigma <- abs(matrix(runif(25,0,1), ncol=5))
sigma <- forceSymmetric(sigma)
sigma <- as.matrix(nearPD(sigma)$mat)
x <- rmvnorm(n, mean=c(0,.25,.15,0,.1), sigma=sigma)
modelMat<-model.matrix(as.formula(~ (x[,1]+x[,2]+x[,3]+x[,4]+x[,5])^3))
beta<-runif(ncol(modelMat)-1,0,1)
beta<-c(2,beta) # setting intercept
mu <- 1-plogis(modelMat%*%beta) # true underlying risk of the outcome
y<-rbinom(n,1,mu)

hist(mu);mean(y)

x<-data.frame(x)
D<-data.frame(x,y)

head(D)

# Specify the number of folds for V-fold cross-validation
folds = 5
## split data into 5 groups for 5-fold cross-validation 
## we do this here so that the exact same folds will be used in 
## both the SL fit with the R package, and the hand coded SL
index <- split(1:n,1:folds)
splt  <- lapply(1:folds,function(ind) D[index[[ind]],])
# view the first 6 observations in the first [[1]] and second [[2]] folds
head(splt[[1]])
head(splt[[2]])

#-------------------------------------------------------------------------------
# Fit using the SuperLearner Package
#-------------------------------------------------------------------------------

## reset the seed
set.seed(NULL)
set.seed(321)

# Fit using the SuperLearner package, specify 
#		outcome-for-prediction (y), the predictors (x), the loss function (L2),
#		the library (sl.lib), and number of folds 
fitY1 <- SuperLearner(Y = y, verbose = T,
                      X = data.frame(x), 
                      method = "method.NNLS", 
                      family = "binomial",
                      SL.library = "SL.glmnet",
                      cvControl=list(V = folds,
                                     shuffle = F, 
                                     validRows = index))

# View the output: 'Risk' column returns the CV-MSE estimates
#		'Coef' column gives the weights for the final SuperLearner (meta-learner)
fitY1
# Now predict the outcome for all possible x 'doses'
yS<-predict(fitY1,
            newdata=data.frame(x),
            onlySL=T)$pred

# Create a dataframe of all x 'doses' and predicted SL responses
Dl1<-data.frame(x,yS)

Dl1$seed = 123

## reset the seed
set.seed(NULL)
set.seed(876543)

# Fit using the SuperLearner package, specify 
#		outcome-for-prediction (y), the predictors (x), the loss function (L2),
#		the library (sl.lib), and number of folds 
fitY2 <- SuperLearner(Y = y, verbose = T,
                      X = data.frame(x), 
                      method = "method.NNLS", 
                      family = "binomial",
                      SL.library = "SL.glmnet",
                      cvControl=list(V = folds, 
                                     validRows = index))

# View the output: 'Risk' column returns the CV-MSE estimates
#		'Coef' column gives the weights for the final SuperLearner (meta-learner)
fitY2
# Now predict the outcome for all possible x 'doses'
yS<-predict(fitY2,
            newdata=data.frame(x),
            onlySL=T)$pred

# Create a dataframe of all x 'doses' and predicted SL responses
Dl2<-data.frame(x,yS)

Dl2$seed = 876543

# examine fits
fitY1
fitY2

## combine results
plot_dat <- rbind(Dl1, Dl2)

## plot
ggplot(plot_dat) +
  geom_point(aes(x = X1, y = yS, 
                 group = factor(seed), 
                 color = factor(seed),
                 shape = factor(seed)), 
             alpha = .5) +
  scale_color_manual(values = c("red","blue"))

plot(Dl1$yS, Dl2$yS)


X <- model.matrix(~-1 + ., data.frame(x))
newX <- model.matrix(~-1 + ., data.frame(x))
Y = y


set.seed(NULL)
set.seed(876543)
mod1 <- glmnet::cv.glmnet(x = X, y = Y, 
                  lambda = NULL, type.measure = "deviance", nfolds = 10, 
                  family = "binomial", alpha = 1, nlambda = 100)
pred1 <- predict(mod1, newx = newX, 
                type = "response", 
                s = "lambda.min")
set.seed(NULL)
set.seed(123)
mod2 <- glmnet::cv.glmnet(x = X, y = Y, 
                  lambda = NULL, type.measure = "deviance", nfolds = 10, 
                  family = "binomial", alpha = 1, nlambda = 100)
pred2 <- predict(mod2, newx = newX, 
                 type = "response", 
                 s = "lambda.min")


mod1
mod2

plot(pred1, pred2)














