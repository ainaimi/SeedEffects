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
  SuperLearner
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
n=1000
x = runif(n,0,8)
y = 5 + 4*sqrt(9 * x)*as.numeric(x<2) + as.numeric(x>=2)*(abs(x-6)^(2)) + rlaplace(n)

# to plot the true dose-response curve, generate sequence of 'doses' from 0 to 8 at every 0.1,
#	then generate the true outcome
xl<-seq(0,8,.1)
yl<-5 + 4 * sqrt(9 * xl)*as.numeric(xl<2) + as.numeric(xl>=2)*(abs(xl-6)^(2))

D<-data.frame(x,y)  	# observed data
Dl<-data.frame(xl,yl)   # for plotting the true dose-response curve

# Specify the number of folds for V-fold cross-validation
folds= 5
## split data into 5 groups for 5-fold cross-validation 
## we do this here so that the exact same folds will be used in 
## both the SL fit with the R package, and the hand coded SL
index<-split(1:1000,1:folds)
splt<-lapply(1:folds,function(ind) D[index[[ind]],])
# view the first 6 observations in the first [[1]] and second [[2]] folds

#-------------------------------------------------------------------------------
# Fit using the SuperLearner Package
#-------------------------------------------------------------------------------

## reset the seed
set.seed(NULL)
set.seed(123)

# Fit using the SuperLearner package, specify 
#		outcome-for-prediction (y), the predictors (x), the loss function (L2),
#		the library (sl.lib), and number of folds 
fitY1<-SuperLearner(Y=y,
                   X=data.frame(x), 
                   method="method.NNLS", 
                   SL.library = "SL.nnet",
                   cvControl=list(V= folds,
                                  validRows=index))

# View the output: 'Risk' column returns the CV-MSE estimates
#		'Coef' column gives the weights for the final SuperLearner (meta-learner)
fitY1
# Now predict the outcome for all possible x 'doses'
yS<-predict(fitY1,
            newdata=data.frame(x=xl),
            onlySL=T)$pred

# Create a dataframe of all x 'doses' and predicted SL responses
Dl1<-data.frame(xl,yS)

Dl1$seed = 123

## reset the seed
set.seed(NULL)
set.seed(876543)

# Fit using the SuperLearner package, specify 
#		outcome-for-prediction (y), the predictors (x), the loss function (L2),
#		the library (sl.lib), and number of folds 
fitY2<-SuperLearner(Y=y,
                   X=data.frame(x), 
                   method="method.NNLS", 
                   SL.library = "SL.nnet",
                   cvControl=list(V= folds,
                                  validRows=index))

# View the output: 'Risk' column returns the CV-MSE estimates
#		'Coef' column gives the weights for the final SuperLearner (meta-learner)
fitY2
# Now predict the outcome for all possible x 'doses'
yS<-predict(fitY2,
            newdata=data.frame(x=xl),
            onlySL=T)$pred

# Create a dataframe of all x 'doses' and predicted SL responses
Dl2<-data.frame(xl,yS)

Dl2$seed = 876543

# examine fits
fitY1
fitY2

## combine results
plot_dat <- rbind(Dl1, Dl2)

## plot
ggplot(plot_dat) +
  geom_point(aes(x = xl, y = yS, 
                 group = factor(seed), 
                 color = factor(seed),
                 shape = factor(seed)), 
             alpha = .5) +
  scale_color_manual(values = c("red","blue"))