##################################################################################

   ##########################################################################
   #                   BIOVEG - Statistical analyses in R                   #
   #                                                                        #
   #     Écio Souza Diniz - Universidade Federal de Viçosa, Brazil          #
   #     Jan Thiele - Thünen Institute of Biodiversity, Germany             #
   #                                                                        #
   #         contact: eciodiniz@gmail.com and jan.thiele@thuenen.de         #    
   #                                                                        #  
   ##########################################################################    


############################### RECOMMENDATIONS ####################################

# 1º Complete an introductory R course  to understand the basic operations used 
# in the software: e.g. how to import files, attach data, install packages and 
# create objects and graphs.

# 2º It is important to have a basic understanding of the principles and functions 
# of the statistical analysis you wish to perform. 

# 3º If your analysis does not execute successfully using these syntaxes, it is 
# useful to consult R forums (e.g. Stack Overflow and Stack Exchange – Cross 
# Validated) check if there is something wrong, missing or if this is simply not 
# the appropriate analysis for your data.  

# 4º Be sure about your choice of analysis.  Refer to existing literature to help 
# decide which analysis is most appropriate for your data. 
 
# 5º If possible, it is highly recommended to read the references at the end of  
# this script, especially those which are highlighted.

# 6º Try not to be overwhelmed by the various methodologies of statistical analysis 
# which might be used to investigate the data.  Again, consult with existing 
# literature to help decide which analyses are most suitable. 
 
# 7º Remember that this script was not created by statisticians.  The adaptable 
# nature of the R means that there are often several ways to complete the same 
# analysis.  This script describes one way to perform typical statistical analysis 
# of the kind of data commonly generated in biological research. 



################################ ANALYSES: Content #################################

### NORMALITY TEST
  #  Shapiro Wilk
### DATA TRANSFORMATION
  #  Log, square, inverse and Arcsine square root
### ONE SAMPLE TESTS
  #  T-test and Wilcoxon Signed Rank Test 
### TWO SAMPLE TESTS
  #  Independent samples: Wilcoxon Rank-Sum Test (Mann-Whitney U-test) and T-test
  #  Dependent samples: Wilcoxon Rank-Sum Test and Paired T-test
### ANALYSES OF VARIANCE
  #  Anova and Tukey post-hoc / Kruskal-Wallis and Dunn's test post-hoc
### INDEPENDENCE TEST 
  #  Chi-squared
### CORRELATION ANALYSES
  #  Pearson, Spearman and VIF (Variance Inflation Factor)
### DIAGNOSTIC FOR LINEARITY
  #  Component+Residual (Partial Residual) Plots
### REGRESSION MODELS
  #  LM (Linear Model)
  #  GLM (Generalized Linear Models
  #  GLMM (Generalized linear mixed models)
  #  GLMM-PQL (Fit GLMM using Penalized Quasi-Likelihood - PQL)
  #  Logistic Regression
  #  SPATIAL AUTOCORRELATION TEST: Moran's I
  #  LME (Linear mixed models)
  #  GLS (Generalized Least squares)
  #  QUADRATIC MODEL
  #  GNM (Generalized Non-Linear Models)
  #  GAM (Generalized Aditive Models)
  #  GAMM (Generalized Aditive Mixed Models)
### SUPERVISED LEARNING
  #  Random Forest
### MODELS SELECTION
  #  Akaike Criterion: AIC and QAIC
### MODEL AVERAGING
  #  Average of models selected by AIC or QAIC
### POST-HOC TEST TO REGRESSION MODELS 
  #  Tukey's HSD (honest significant difference) 
  #  LmerTest (mixed models)
### MULTIVARIATE ANALYSIS
  #  NMDS, ANOSIM, CCA, PCA, PERMANOVA and INDVAL (Indicator species analysis)



################################### Read test datasets #############################

# BIOVEG = your source (file) of data. 
# Exemplo: 

setwd("C:/<WorkDirectory>")

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

################################### NORMALITY TEST #################################

# One of the most commonly used tests in Botany is Shapiro-Wilk. This test verifies 
# whether one specific variable from your data is normally (gaussian) distributed. 
# If the p-value is below the critical level (e.g. below 0.05) the variable can be 
# described as deviating significantly from the normal (gaussian) distribution.


# usage:

shapiro.test(BIOVEG$Basal.area)



############################### DATA TRANSFORMATION ################################

# In some cases it is necessary to transform your data, either to achieve normality 
# or to get a better fit for the residuals of your model. Here the  mention some of 
# the most commonly applied transformations, but there are many others (e.g. Box cox) 
# and it's always good to read about other possibilities to maximize the fit and 
# robustness of your analysis.


# plus: 

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

Basal.plus <- BIOVEG$Basal.area + 5
cbind(BIOVEG$Basal.area, Basal.plus) # just to have a look at the result

# log: 

# mainly indicated to try to make right-skewed data normal and by 
# default included in Poisson GLM (i.e. the default link function is log). But,
# many authors have suggested not to transform count (poisson) data. Read about
# it if possible. 
# O.B.S.: in case your data contain zeros or negative values, you have to add a
# constant to all values so that all of them are positive before you can take 
# the logarithm.
 
Basal.log <- log(BIOVEG$Basal.area) 
cbind(BIOVEG$Basal.area, Basal.log)

x <- rnorm(10, mean = 2, sd = 2) # toy example of a variable with values =< 0
log(x) # does not work
c <- abs(min(x)) + 0.1 # the constant you have to add to x to make all values > 0
x.log <- log(x + c) # works
cbind(x, x.log)

# square:
 
Basal.sqrt <- sqrt(BIOVEG$Basal.area)
par(mfrow=c(2,1)) # let's check the transformation using histograms
hist(BIOVEG$Basal.area)
hist(Basal.sqrt)

# Inverse:
 
Basal.inv <- 1 / BIOVEG$Basal.area

# Arcsine square root:
 
# Can be used for data that range between zero and 1. Often indicated for 
# binomial data (proportions/ percentages), but it is preferable to use a model 
# with binomial distribution rather than an arcsine square root transformation 
# in this case.

Hogweed <- read.table("Hogweed.csv", header=T, sep=";", dec=".")
head(Hogweed)

Hogweed.asr <- asin(sqrt(Hogweed$Hogweed.cover/100))*2/pi 
plot(Hogweed$Hogweed.cover/100, Hogweed.asr)
par(mfrow=c(2,1))
hist(Hogweed$Hogweed.cover/100)
hist(Hogweed.asr)



################################ ONE SAMPLE TEST ####################################

###### Parametric ##########

### T-TEST - One sample ######

# Indication: to test one sample mean against zero to verify whether there is 
# a significant difference.

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# Usage:

t.test(BIOVEG$Basal.area, mu=0) # test against zero



###### Nonparametric ##########

### Wilcoxon Signed Rank Test ######

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# Usage:

wilcox.test(BIOVEG$Basal.area, mu=0)



################################ TWO SAMPLE TEST ######################################


###### Independent Two samples: Parametric ##########

### T-TEST - Two sample ###

# Indication: to compare two independent samples. 

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# usage: 

t.test(BIOVEG$Basal.area[BIOVEG$type=="Aluvial"], 
       BIOVEG$Basal.area[BIOVEG$type=="Ombrofila"],alternative="two.sided")

t.test(BIOVEG$Basal.area[BIOVEG$type=="Aluvial"], 
       BIOVEG$Basal.area[BIOVEG$type=="Ombrofila"], alternative="less", 
       var.equal=TRUE)


### Two Dependent samples: Parametric ######

### PAIRED T-TEST - Two sample ####

# Indication: to compare two dependent samples. 


t.test(BIOVEG$Basal.area[BIOVEG$type=="Aluvial"], 
       BIOVEG$Basal.area[BIOVEG$type=="Ombrofila"],paired = TRUE)




###### Independent samples: Nonparametric #########


### Wilcoxon Rank-Sum Test (or Mann-Whitney U-test) ######


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# Usage:

## independent 2-group Mann-Whitney U Test

# where y is numeric and A is a binary factor (presence or absence, 0 or 
# 1, or similar) we can use the formula notation in wilcox.test

BIOVEG$Aluvial <- ifelse(BIOVEG$type=="Aluvial", 1, 0) # create binary factor
wilcox.test(BIOVEG$Basal.area ~ BIOVEG$Aluvial) 

# independent 2-group Mann-Whitney U Test; here we specify the two groups to be  
# compared using indexing (square brackets) 

wilcox.test(BIOVEG$Basal.area[BIOVEG$type=="Aluvial"], 
            BIOVEG$Basal.area[BIOVEG$type=="Ombrofila"], paired=FALSE)




###### Dependent Two samples:  Nonparametric ##########

######## Wilcoxon Rank-Sum Test ######

# Indication: the same function as the Wilcoxon test and two sample T-test, 
# but addressing paired samples. This applies to data which has been collected on 
# the same plots, for example, either on the same date or at different dates (e.g. 
# comparison between basal area in 2010 and 2014. The basal area of 2014 in general 
# is dependent on how much was accumulated in 2010).

# Another test dataset describes the soil nutrient content of a number of different 
# plots.  The quantities of each nutrient (N, P, K) were derived from a single soil 
# sample per plot.  


Hogweed <- read.table("Hogweed.csv", header=T, sep=";", dec=".")
head(Hogweed)


# Usage: 

wilcox.test(Hogweed$P_mg_100g, Hogweed$K_mg_100g, paired=TRUE) 

# When you have dependent samples, in "paired" insert TRUE.  For independent samples 
# insert FALSE.



################################ ANALYSIS OF VARIANCE ###############################



# Note: All Analyses below relate to variance, one-sample and two-sample. 
# The majority of this analysis can be completed using the basic pre-loaded 
# built-in "stats" package.  To perform Dunn's test (Post hoc), which follows 
# Kruskal-Wallis,the "dunn.test" package must be installed. 



###### Parametric tests: One-Way Anova ########

# Indication: tests for differences between the means of more than two groups 
# e.g. Diameter in treatment1, Diameter in treatment2, Diameter in treatment3)


#usage:

# Note: transformations of the dependent variable, e.g. log, can be specified 
# in the syntax of ANOVA and other functions directly. 

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

aov(log(Diameter) ~ type, data=BIOVEG)
results <- aov(log(Diameter) ~ type, data=BIOVEG)
summary(results)



## Tukey test (post hoc) ###

# Indication: to identify when there is a significant difference between categories
# previously shown as overall result by parametric ANOVA.

# Usage:

a1 <- aov(log(BIOVEG$Diameter) ~ BIOVEG$type)
posthoc <- TukeyHSD(x=a1, 'BIOVEG$type', conf.level=0.95)
print(posthoc)
plot(posthoc)



###### Nonparametric tests: Kruskal-Wallis ########

# Indication: use following a significant result in the Kruskal-Wallace test 
# (a non-parametric analysis, the underlying distribution is not assumed in 
# advance). 
 


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

#Usage:

r4<- kruskal.test(BIOVEG$Basal.area ~ BIOVEG$type)
r4


## Dunn's Test (post hoc) ###

# It performs the same function as a TukeyHSD test following ANOVA, but here is 
# following Kruskal-Wallis. 
 
# Usage:

install.packages("dunn.test")
library(dunn.test)

dunn.test(BIOVEG$Basal.area, g=BIOVEG$type, kw=TRUE, label=TRUE, wrap=TRUE, 
          table=TRUE, list=FALSE, rmc=FALSE, alpha=0.05)


# when the command line extends beyond more than one row in the text, be sure 
# to select all rows to run the complete command.




########################## Independence test - Chi-squared ##########################

# The Chi-squared test is one of the most commonly used tests in biology and
# biomedicine because it can test if one variable is dependent or not on a certain 
# category.  In the example below we use the phenology stages of a species and the 
# abundance of individuals counted at each stage over 4 months (Feb, Mar, Apr, May). 
# Our aim is to check if the abundances at each stage are dependent on the months. 
# In case Chi-squared is significant and, thus, shows dependence, you need to split 
# the table into two in a 2 x 2 format and run the test again. This post-hoc test is 
# called: Partition of Chi-squared. Further Information about this test is widely 
# available online.

# install and load the MASS package 



install.packages("MASS")
library(MASS)       

# to build a contingency from columns of a data table you can use:
tbl <- table(Data$flowering, Data$fruits, Data$Month)

# Here, we will just create a contingency table "manually":

tbl <- rbind(c(7,1,3,9), c(8,5,2,0))
colnames(tbl) <- c("Feb", "Mar", "Apr", "May")                 
rownames(tbl) <- c("Flowering", "Fruits")
tbl

# Usage: 

chisq.test(tbl) # normal p-value of Chi-squared test

chisq.test(tbl, simulate.p.value=TRUE) # simulated p-value (Monte Carlo test)

# "Posthoc" test comparing only two month

tbl.1 <- tbl[,c("Mar", "Apr")]
tbl.1
chisq.test(tbl.1)



############################### CORRELATION ANALYSIS ################################



### Pearson's correlation test ###


## use PerformanceAnalytics package

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# Usage:

library(PerformanceAnalytics)
chart.Correlation(BIOVEG[,c(5, 6, 7, 8)]) 
correlacao<- cor(BIOVEG[,c(5, 6, 7, 8)]) 
round(correlacao, digits=2)

### Correlation test between just two variables ###

## use the default R function "cor.test"

cor.test(BIOVEG$Richness, BIOVEG$Basal.area)


# The number between parentheses are the column numbers of your data file 
# (e.g. Excel file). This is how you select the variables whose correlation 
# you want to check.

# REMEMBER: correlation just shows you if one data value increases/ decreases 
# with the other, but it does not say that they are dependent on each other. 
# To check dependence between variables you need to test their relation and for 
# that there are regression models. However, even regression models do not prove 
# that there is a causal relationships between the variables. Correlation and 
# regression cannot discern between real causal effects and spurious correlations 
# (that are caused by a third variable/ other variables).


### Spearman's correlation test ###

# There is a perpetual discussion among statisticians about whetherPearson 
# analysis assumes strictly normality, but a large majority of high impact 
# literature continues to use it, even in non-normal cases. Some researchers 
# recommend Spearman, but, for a number of reasons, its application is more 
# limited than Pearson. Please read about this subject further if it applies 
# to your analysis.

# If you want to use Spearman rank correlation just use the 'method' argument:


chart.Correlation(BIOVEG[,c(5, 6, 7, 8)], method="spearman")
correlacao<- cor(BIOVEG[,c(5, 6, 7, 8)], method="spearman") 
round(correlacao, digits=2)

### Correlation test between just two variables ###

## use the default R function "cor.test"

cor.test(BIOVEG$Richness, BIOVEG$Basal.area, method="spearman")



### VIF - Variance Inflation Factor ###

# Among the most robust approaches used nowadays, the variance inflation factor 
# (VIF) is the quotient of the variance in a model with multiple terms (predictors) 
# by the variance of a model with one term alone. This is a great method to quantify
# the severity of multicollinearity in an ordinary least squares regression analysis
# (e.g. LM, LME, GLS, etc). VIF provides an index that returns an outcome measuring
# how much the variance (the square of the estimate's standard deviation) of an 
# estimated model coefficient is increased because of collinearity.

# In Ecology, an usual threshold to consider collinearity among predictors acceptable
# is VIF < 0.5 

# For further details and information on the statistical backbone and potential of
# VIF, please read the great revision and simulations of Dormann et al. 2013 and
# Zurr et al.2010.

# Dormann et al. Collinearity: A review of methods to deal with it and a simulation 
# study evaluating their performance. Ecography, 36(1): 27-46, 2013.

# zuur AF, Ieno EN & Elphick CS. A protocol for data exploration  to  avoid  common  
# statistical  problems.  Methods in Ecology and Evolution 1: 3–14, 2010. 


BIODIV.vif<- lm(Species.richness ~ P_mg_100g + Hogweed.cover + N_perc +
                  Habitat.type, data=Hogweed)

vif(BIODIV.vif)

# The outcome shows that all predictors have VIF < 5, thus they can be kept together
# in the same model.



########################### Diagnostic for linearity ################################

# To check for linear relationships among dependent (Y) and predictor variables (X) 
# we can use the function "crPlots" of the package "car". This function returns a 
# graph that bases the diagnostic on the residuals of the model. The more closely 
# the two lines (solid and dashed) fit each other, the more linear is the relationship
# among Y and X. 

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

library(car) # graphical diagnostic for linearity

DBH<- glm(Hogweed.cover ~ P_mg_100g + K_mg_100g + Veg.cover, data=Hogweed, 
                                                       family="gaussian")
crPlots(DBH, col.lines=c("red", "green"))


# In the diagnostic we can see that the predictor "Veg.cover" shows the best, and an 
# acceptable, fit for linearity. This is indicated by the close fit between the red 
# line (best fit) and green line. Conversely, the linearity of "P_mg_100g" and 
# "K_mg_100g" deviates considerably.


###################################### REGRESSION ###################################

# Indication: Regression analysis is a statistical approach for estimating the 
# relationships between a dependent variable (axis Y and commonly called the 
# 'outcome variable') and one or more independent variables (predictors).There
# are linear and nonlinear regression models, although the linear ones are the
# most common form of used regression analysis.


######### LM (Linear Model) ###############

# This is the most famous type of regression of analysis and relies in the 
# assumptions of the Ordinary least squares: 
# normal distribution,  homoscedasticity, linearity between, predictors and 
# dependent variable, weak multicollinearity among predictors, independence of 
# errors and residuals not correlated.

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)


install.packages("car")
library(car)

# Usage: 

BIODIV.pre<- lm(Species.richness ~ P_mg_100g + Habitat.type, data=Hogweed)
summary(BIODIV.pre)
Anova(BIODIV.pre)
par(mfrow=c(2,2))
plot(BIODIV.pre) # check diagnostic plots
BIODIV.lm.aov <- anova(BIODIV.pre, test="F")
BIODIV.lm.aov

# Residuals

x11(width=12, height=12)
par(mfrow=c(2,2))
plot(BIODIV.pre) 
hist(resid(BIODIV.pre),breaks=20) 

# Scatterplot graph: 

b.0.BIODIV<- BIODIV.pre$coefficients[1]
b.1.BIODIV<- BIODIV.pre$coefficients[2]
plot(Hogweed$Species.richness ~ Hogweed$P_mg_100g, las=1, 
     ylab="Species richness", xlab="P (mg/100 g)", cex.lab=1.1)
curve(b.0.BIODIV + b.1.BIODIV*x, add=TRUE, col="red", lwd=2)
mtext(ifelse(BIODIV.lm.aov[1,5]<0.001, "p < 0.001", 
      paste("p = ", round(BIODIV.lm.aov[1,5],3))))



# To construct the curve (or line) in this LM example, the dependent variable, 
# Species richness, is on the y axis and the metric predictor variable, phosphorous, 
# is on the x axis. The regression line is calculated using the intercept (b.0.BIODIV) 
# plus the estimate of phosphorous (b.1.BIODIV) multiplied by x. The insertion of the 
# P-value of ANOVA is placed according to the location of values in the ANOVA output. 
# That is, the "P_mg_100g" is being tested in relation to "Species.richness" and it's 
# located in the first row of the ANOVA output and its P-value on the 5th column in  
# this row. Thus: BIODIV.lm.aov[1,5]. If you have done any data transformation before 
# (log, square ...) you must include it in the syntax of the scatterplot, e.g. 
# log(Hogweed$Species.richness). If the variable include zeros, add a constant, 
# e.g. +1, before taking the log.


########################## GLM (Generalized Linear Models) #############################

# Note: in GLM you must specify the distribution family of the dependent variable. 
# Thus:


# Gaussian (normal) distribution: Continuous data (e.g. body weight, basal area, 
# height). If the residuals don't satisfactorily fit the Gaussian distribution, you 
# can try Gamma distribution. BUT, first of all read about both kinds of distribution 
# before you decide which to use . Gamma doesn't accept negative values.

# Poisson: count data (e.g. species richnees). First check if the Poisson model shows 
# overdispersion. If it does, try "quasipoisson" distribution. 

# Binomial: this distribution is often used to model the number of successes (1) or 
# failure (0), but also presence (1) or absence (0), in a sample. Furthermore, data
# regarding percentage or proportion in a sample.

# Negative BinomiaL: count data in case of overdispersion in "poisson" or "quasipoisson".

# A Poisson or Binomial distribution model might be likely to show under- or 
# overdispersion, indicating that your model fit is not satisfactory. After running the 
# "summary", check for evidence of overdispersion. You do this by checking residuals 
# and degrees of freedom: residual deviance/degrees of freedom << 1 (e.g. < 0.7) 
# indicate underdispersion, while residual deviance/degrees of freedom >> 1 (e.g. > 1.3) 
# indicate overdispersion. In such cases, you need to use "quasipoisson" or "quasibinomial" 
# distributions instead. The model will then estimate a dispersion parameter that adjusts 
# the variation in the data to the fixed variation of the Poisson or Binomial distribution.
# If you are analyzing count data (integers), you can also use the Negative Binomial 
# distribution instead of "quasipoisson". To use Negative Binomial distribution, you use 
# the function "glm.nb", which is the function for negative binomial glm from the package 
# MASS.

# To use Negative Binomial distribution, you just insert like this example: 
# "glm.nb" is the function for negative binomial glm from the package MASS.


install.packages(“MASS”)
library(MASS)


Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

riqueza.glm<- glm.nb(Species.richness ~ P_mg_100g + Land.use, data=Hogweed) 
summary(riqueza.glm)
Anova(riqueza.glm)


### GLM usage ###:

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

riqueza.glm<- glm(Species.richness ~ P_mg_100g + Land.use, data=Hogweed,
                  family="poisson")
summary(riqueza.glm)
library(car) 
riqueza.aov <- Anova(riqueza.glm)
riqueza.aov


## Checking the residuals and see if it does fit good to the model 

# residuals

par(mfrow=c(2,2))
plot(fitted(riqueza.glm), resid(riqueza.glm))
qqnorm(resid(riqueza.glm))
hist(resid(riqueza.glm))

# Boxplot

par(mfrow=c(1,1))
boxplot(Hogweed$Species.richness ~ Hogweed$Land.use, las=1, ylab="riqueza")
mtext(ifelse(riqueza.aov[2,3]<0.001, "p < 0.001", 
paste("p = ", round(riqueza.aov[2,3],3))))

# Scatterplot graph: 

b.0.riqueza<- riqueza.glm$coefficients[1]
b.1.riqueza<- riqueza.glm$coefficients[2]
plot(Hogweed$Species.richness ~ Hogweed$P_mg_100g, las=1, ylab="riqueza", 
     xlab="P (mg/ 100g)", cex.lab=1.1)
curve(exp(b.0.riqueza + b.1.riqueza*x), add=TRUE, col="red", lwd=2)
mtext(ifelse(riqueza.aov[1,3]<0.001, "p < 0.001", 
      paste("p = ", round(riqueza.aov[1,3],3))))


# Note 3.: the graph construction follows the principles described with LM previously. 
# If you’ve used a transformation or a link function (other than "=") in the model, 
# you can include it in the scatterplot syntax, e.g. log(Species.richness). But if you 
# want to plot the original values of the dependent variable, you can do that, too, you 
# just need to back-transform the model predictions to the original scale. In this case 
# you need to transform from the log-scale (used in the Poisson model) to the untransformed 
# species richness used in the scatterplot. That's why we need the "exp" function in the 
# example above.



##################### GLMM (Generalized Linear Mixed Models) ##########################

# Note 1: Like GLM, the GLMM also require to specify the distribution family (e.g. 
# Poisson), but they don't accept "quasipoisson" or "quasibinomial". If poisson 
# doesn't fit the data because of overdispersion, try negative binomial distribution  
# as explained in GLM section. The syntax is the same: just insert the argument "nb". 
# Example: basalarea1<- glmer.nb(..., Do not write family name). 


# GLMM indication: when you have pseudoreplicates. It's very common in ecological 
# field research to use subplots within each site and these subplots cannot be 
# considered  as real replicates. So, a mixed model is appropriate to treat this 
# situation. The definition of what variable is a fixed effect (e.g. in the example 
# below it is "Hogweed.cover" (cover of giant hogweed) and "Habitat.type") and what 
# is random effect (e.g. in the example below "Study.area") depends on your research 
# questions and study design. Usually, design factors, such as "block" or "study area" 
# will be random effects, and the variables that you are really interested in 
# (treatments, soil parameters etc.) will be fixed effects. BUT, read about random 
# and fixed effects before constructing your model. 

# When the command line doesn't end in the same row, but continues in the next one 
# (like the first command of the GLMM example below), select and press the command of 
# both two lines together. We broke some command lines here  to fit in the R script 
# style.


### GLMM usage:


Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

library(lme4) # glmer
library(car) # Anova

Richness <- glmer(Species.richness ~ scale(Hogweed.cover) + Habitat.type + 
                  (1|Study.area), data=Hogweed, family="poisson")
# Here, the metric predictor variable (Hogweed.cover) was scaled to avoid 
# problems in the calculation of the model

summary(Richness)
Richness.aov <- Anova(Richness, Type="II") 
Richness.aov

# Residuals

par(mfrow=c(2,2))
plot(fitted(Richness), resid(Richness))
qqnorm(resid(Richness))
qqline(resid(Richness))
hist(resid(Richness))



### Scatterplots ### 

# With the syntax below you can construct graphs of the relation of the depenent 
# variable with a metric predictor variable. This can be done for all levels of 
# the categorical predictor variable (here: Habitat.type). To construct the graph 
# you can just select and run all of the command lines below. 

b0<- coef(summary(Richness))[1,1] # intercept (habitat "agricultural grassland")
b1<- coef(summary(Richness))[2,1] # estimate of hogweed cover
b2<- coef(summary(Richness))[3,1] # estimate of habitat "ruderal grassland"
b3<- coef(summary(Richness))[4,1] # estimate of habitat "tall herbs"
b4<- coef(summary(Richness))[5,1] # estimate of habitat "wasteland"
b5<- coef(summary(Richness))[6,1] # estimate of habitat "woodland"

x11(width=14, height=12)

plot(Hogweed$Species.richness ~ Hogweed$Hogweed.cover, las=1, 
     ylab="Species richness", xlab="Cover of giant hogweed", cex.lab=1.3) 


# The model that we apply here used poisson distribution which, by default, 
# includes a log-link (i.e. log transformation of the dependent variable). 
# Thus, we need to use exp() to back-transform the values predicted by the model, 
# to the original values of species richness. 


mean <- mean(Hogweed$Hogweed.cover)
sd <- sd(Hogweed$Hogweed.cover)

curve(exp(b0 + b1*(x-mean)/sd), add=TRUE, col="green") # curve for agricultural grassland
curve(exp(b0 + b2 + + b1*(x-mean)/sd), add=TRUE, col="cyan") # curve for ruderal grassland
curve(exp(b0 + b3 + + b1*(x-mean)/sd), add=TRUE, col="blue") # curve for tall-herb stands
curve(exp(b0 + b4 + + b1*(x-mean)/sd), add=TRUE, col="orange") # curve for wasteland
curve(exp(b0 + b5 + + b1*(x-mean)/sd), add=TRUE, col="darkgreen") # curve for woodland
mtext(ifelse(Richness.aov[1,3]<0.001, "p < 0.001", 
      paste("p = ", round(Richness.aov[1,3],3))))


# Note 2: the graph construction follows the principles previously described for LM. 
# For further details type the following into your internet search engine: "how construct
# scatterplots for Mixed models - GLMM? If you have transformed your data, include it in 
# the scatterplot syntax, e.g. log(Data.h3$area basal+1), OR back-transform the predicted 
# values to the original scale (as in the previous example ).




############## GLMM-PQL (Fit GLMM using Penalized Quasi-Likelihood - PQL) ###############

# This type of GLMM is used when you need to compute a mixed model for a non-normal 
# dispersion (e.g. poisson,  binomial) , it is necessary to use a generalized mixed 
# model that allows the inclusion of correlation structures to account for spatial 
# or temporal effects in the random effects and residuals of the model. The usual 
# functions lmer and glmer of the package lme4 do not allow you to include such 
# correlation structures. Thus, you can adjust Partial Quasi-likelihoods for glmm with 
# the function glmmPQL of the package MASS. The most used spatial and temporal  
# correlation structures for a wide range of situations are, respectively, corExp 
# (exponential spatial correlation structure) and corAR1 (correlation structure of first 
# order for time  series) of the package nlme.

# The glmmPQL approximates a quasi-likelihood by iterative fitting of (re)weighted  
# linear mixed effects (lme) models based on the fit of Glm. Thus, it estimates the  
# fixed effects parameters by fitting a Glm with its incorporated correlation 
# (variance-covariance) structure. This constructs an Lme model and refits it to 
# re-estimate the variance-covariance structure, utilizing the variance structure from 
# the previous Glm. The iterations continue until the fit improvement is below a   
# threshold or a defined number of iterations has occurred. Therefore, this approach  
# accommodates heterogeneity and spatial/temporal autocorrelation in the model. However, 
# it worth being aware of the fit of PQL generalization for Poisson models, since it  
# commonly performs poorly. If your performance with glmmPQL is poor using the Poisson  
# model, take a look at glmmTMB and glmmADMB packages and functions.

# For more information on usage of glmmPQL, other models and correlation structures to 
# account for spatial and temporal autocorrelation of residuals, read Dormann et al. 
# (2007).



########## GLMM - PQL: Temporal autocorrelation ###########

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

# install.packages("MASS") 
# install.packages("nlme")
# install.packages("lme4")


library(MASS) # PQL Estimation Of Generalized Linear Mixed Models (glmmPQL)
library(nlme) # corAR1 (autocorrelation structure of order 1. 
              # For temporal autocorrelation correction)
library(lme4) # glmer

RichnessPQL1 <- glmmPQL(Species.richness ~ scale(Hogweed.cover) + 
                         Habitat.type, random=~1| Study.area, data=Hogweed,
                        family="poisson", correlation=corAR1())

summary(RichnessPQL1)
class(RichnessPQL1)="lme"
anova(object=RichnessPQL1,test="Chisq")


# Residuals

par(mfrow=c(1,1))
plot(fitted(RichnessPQL1), resid(RichnessPQL1))
qqnorm(resid(RichnessPQL1))
qqline(resid(RichnessPQL1))
hist(resid(RichnessPQL1))

# Note: Fitting a "lme" class for the glmmPQL is not recommended as the most appropriate 
# way to acquire p-values from anova of the model above. Here we can compute a lme model 
# for the same variables of the glmmPQL model and compare the outcomes of the anova of 
# both models. Despite lme requiring data to conform to assumptions of normality, we 
# can use this model for a dataset where Y is discrete (Poisson), once the residual fit 
# is as close as possible to normal, which were already complied using the residuals of 
# glmmPQL model.. Let's see.

RichnessPQL1.2 <- lme(Species.richness ~ scale(Hogweed.cover) + 
                         Habitat.type, random=~1| Study.area, 
                         data=Hogweed, correlation = corAR1())

summary(RichnessPQL1.2)
anova(RichnessPQL1.2)

# Residuals

par(mfrow=c(1,1))
plot(fitted(RichnessPQL1.2), resid(RichnessPQL1.2))
qqnorm(resid(RichnessPQL1.2))
qqline(resid(RichnessPQL1.2))
hist(resid(RichnessPQL1.2))


# Overall, we can see that both the lme and glmmPQL models in this example display 
# well-fitting residuals and the anova outcomes of both provide very similar results.



#### Wald Test ##### 

# The recommended way to test hypotheses (H0/H1) for glmmPQL, is to compute the Wald 
# test which assesses  the significance of groups of parameters. The main difference 
# between the Wald approach and a typical anova, is that Wald tests encompass marginal 
# tests of parameters, i.e. as they show up in the summary of the model in Wald outcome 
# versus nested model comparisons (either F or Chisq (likedlihood ratio)) in anova() 
# function.


RichnessPQL1 <- glmmPQL(Species.richness ~ scale(Hogweed.cover) + 
                          Habitat.type, random=~1| Study.area, data=Hogweed,
                        family="poisson", correlation=corAR1())
summary(RichnessPQL1)


# Since the main statistic of interest is the Wald approach (t for PQL) for the fixed 
# effects parameters, we can firstly fit a crude Wald. The summary of the model gives 
# us crude (Wald) estimates of the p-values for each individual parameter: 

coef(summary(RichnessPQL1))

# One could be interested in the investigation of the main effects as a whole. Thus, 
# we use the function "wald.test" of the package aod in order test the influence of 
# the factor as whole using the Wald Z and χ² (Wald Chi-squared)  test, which is an 
# approach analogous to an ANOVA. In function "wald.test"   the following parameters 
# are used: 

# fixed effects parameter estimates (b) the variance-covariance  matrix (Sigma) 
# and (c) a specification of which fixed factor terms (Terms) to combine for the 
# Wald statistic. For more information, please, read the R documentation of the 
# package aod in the wald.test section.

# You define the terms by designating the order in which the coefficients of 
# interest appear in the summary of the model. For instance, the predictor 
# scale(Hogweed.cover) is the second coefficient after the intercept and 
# Habitat.typewoodland is sixth and last. If we are interested in the main effects
# of the groups of parameters as whole, then we write 'Terms=2:6'. It is important 
# to notice that the Wald X² test is suitable only in the absence of overdispersion 
# in the model. When the model is overdispersed the suitable approach is to use Wald
# F and χ2 tests.

# install.packages("aod")
# library(aod)

wald<-wald.test(b=fixef(RichnessPQL1),
                  Sigma=vcov(RichnessPQL1),Terms=2:6)

wald

# The global p-value of Wald confirms the overall significant effect of the numeric 
# (Hogweed.cover) and categorical (Habitat.type) predictors, as found in the other 
# models (glmmPQL and lme), with Species.richness as response variable.

# To assess the confidence intervals for model estimates based on Wald  statistics.

library(gmodels)
ci(RichnessPQL1, method="Wald")


########## GLMM - PQL: Spatial autocorrelation ###########



BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)


library(MASS) # PQL Estimation Of Generalized Linear Mixed Models (glmmPQL)
library(nlme) # corExp (autocorrelation structure of order 1. 
              # For temporal autocorrelation correction)
library(lme4) # glmer


RichnessPQL2 <- glmmPQL(Richness ~ Mortality + Recruitment, random=~1| 
                          Site, data=BIOVEG, family=poisson(), 
                          correlation = corExp())

summary(RichnessPQL2)
class(RichnessPQL2)="lme"
anova(RichnessPQL2,test="Chisq")

# Residuals

par(mfrow=c(1,1))
plot(fitted(RichnessPQL2), resid(RichnessPQL2))
qqnorm(resid(RichnessPQL2))
qqline(resid(RichnessPQL2))
hist(resid(RichnessPQL2))


# The model RichnessPQL2 above did not show an acceptable residual fit in the 
# Q.Q Plot and Histogram of residuals. In this case we then need to calculate 
# a glmmPQL negative binomial. Firstly, we compute a glm model for the same 
# set of fixed effect variables to fit the Theta parameter of negative binomial
# dispersion.

RichnessPQL2.1 <- glm.nb(Richness ~ Mortality + Recruitment, data=BIOVEG)
theta<-RichnessPQL2.1$theta

RichnessPQL2.1.1 <- glmmPQL(Richness ~ Mortality + Recruitment, random=~1| 
                          Site, data=BIOVEG, family=quasipoisson(), 
                        correlation = corExp())

summary(RichnessPQL2.1.1)
class(RichnessPQL2.1.1)="lme"
anova(object=RichnessPQL2.1.1,test="Chisq")

# The outcomes of both the glmmPQL model with Poisson and Negative binomial 
# dispersion were the same.

# Similar to the first example of glmmPQL above, we run a lme model for the same
# set of variables and compare the residuals and anova outcomes of both models. 
# Despite the fact that lme assumes a normal distribution, we can still consider 
# this model if the residuals reasonably conform with a normal fit.

RichnessPQL2.2 <- lme(Richness ~ Mortality + Recruitment, random=~1| 
                        Site, data=BIOVEG, correlation 
                      = corExp())
summary(RichnessPQL2.2)
anova(RichnessPQL2.2)

# Residuals

par(mfrow=c(1,1))
plot(fitted(RichnessPQL2.2), resid(RichnessPQL2.2))
qqnorm(resid(RichnessPQL2.2))
qqline(resid(RichnessPQL2.2))
hist(resid(RichnessPQL2.2))

# The anova output of the lme model, similar to those of the glmmPQL with Richness
# as response variable, shows that none of the predictors are significantly related
# to this response variable. However, the residual fit of this lme model is also not
# as good as those  of the glmmPQL models for this response variable. Thus, we can 
# try to log-transform Richness, though log transformations of discrete data are 
# largely not advised and their use demands a cautious approach.

RichnessPQL2.2.2 <- lme(log(Richness) ~ Mortality + Recruitment, random=~1| 
                        Site, data=BIOVEG, correlation 
                      = corExp())
summary(RichnessPQL2.2.2)
anova(RichnessPQL2.2.2)

# Residuals

par(mfrow=c(1,1))
plot(fitted(RichnessPQL2.2.2), resid(RichnessPQL2.2.2))
qqnorm(resid(RichnessPQL2.2.2))
qqline(resid(RichnessPQL2.2.2))
hist(resid(RichnessPQL2.2.2))

# For this lme model, the log transformation of Richness worked well. The residuals
# are better fitted and the anova output maintains non-significance for the predictor 
# variables, similar to the previous glmmPQL and lme models. This suggests that there 
# is no great distortion of the empirical condition of the data in the model after 
# the log transformation.


#### Wald Test #####

# Now, let's compute the Wald test for the model "RichnessPQL2.1.1"


summary(RichnessPQL2.1.1)

wald2<-wald.test(b=fixef(RichnessPQL2.1.1),
                Sigma=vcov(RichnessPQL2.1.1),Terms=2:3)

wald2

# The global p-value of Wald confirms the non-significant effects of the predictors
# (Mortality and Recruitment) as found in the previous glmmPQL and lme models for 
# this set of variables.

# Confidence Interval under Wald fit.

ci(RichnessPQL2.1.1, method=Wald)


###### Useful links and examples: 


# https://rpubs.com/bbolker/glmmchapter

# http://www.flutterbys.com.au/stats/tut/tut11.2a.html

# https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

# https://bbolker.github.io/mixedmodels-misc/ecostats_chap.html

# https://www.rdocumentation.org/packages/MASS/versions/7.3-47/topics/glmmPQL


###### Recommended Reading: 

# Box, G.E.P., Jenkins, G.M., and Reinsel G.C. (1994) "Time Series Analysis: 
# Forecasting and Control", 3rd Edition, Holden-Day.


# Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", 
# Springer, esp. pp. 235, 397.


# Cressie, N.A.C. (1993), "Statistics for Spatial Data", J. Wiley & Sons.

# Venables, W.N. and Ripley, B.D. (2002) "Modern Applied Statistics with S", 
# 4th Edition, Springer-Verlag.



############################### LOGISTIC REGRESSION #################################


# Hypothesis test: in a period of 5 years do the continuous variables (mean diameter 
# and mean precipitation) affect the likelihood of tree recruitment, i.e., survival, 
# growth and achieve a minimum diameter criteria (e.g. > 15 cm). That is, we want to 
# know if in a five years-period of forest monitoring, do two continuous variables 
# (Mean_diameter + Precipitation_mean) and one categorical (type of vegetation) 
# determine the success of a  trees individual likelihood of recruitment (0 
# indicates the absence of recruitment in a subplot and 1 the  presence).

Lreg <- read.table("Logit2.csv", header=T, sep=";", dec=".")
head(Lreg)
 
Recrut.binar<- glm(Rec.binar ~ Mean_diameter + Precipitation_mean + 
                                type1, data=Lreg,family="binomial")
Recrut.binar

# According to the output, the model is logit(pi) = 
# 4.43 + -0.14*Mean_diameter + 0.001*Precipitation_mean -  -0.82*typeOmbrophilous
# + 0.12*typeSemideciduos


### LRT (likelihood ratio test) ####

# To test the overall model fit and hypothesis 

# LRT  compares the full model with a reduced model where the explanatory variables 
# of interest are omitted and provides p-values of the tests calculated using 
# Chi-squared distribution.

Rec.bin.reduced<- glm(Rec.binar ~ 1, data=Lreg, family="binomial")
anova(Rec.bin.reduced,Recrut.binar, test="Chisq")

# The LRT is 33.399 with a p-value 9.895e-07. Thus, we have strong evidence, based 
# on high significance, that there is influence of some of the predictors in the 
# success or failure of tree recruitment.

summary(Recrut.binar) # perform tests on the individual regression parameters

# We can see that mean diameter (z = -4.029 and p=5.6e-05 ***) appears to have a 
# significant impact on the likelihood of recruitment success or failure of the 
# trees. The same way the vegetation type Ombrophilous (z = -2.628  and p=0.00859**) 
# seems to have significant negative impact (Estimate = -0.825271) on the recruitment 
# success. Semideciduous forest (z = -2.628 and p = 0.00859 **)does not seem to have a 
# significant effect on recruitment, once mean diameter and ombrophilous are included 
# in the model. The same can be said about the effect of mean precipitation (z = -1.505 
# and p = 0.13224) on recruitment, since it's effect is controlled by the stronger 
# influence of mean diameter and ombrophilous type. 

# We can check for significance of influence from the predictors by executing a test 
# for the full model.

anova(Recrut.binar, test="Chisq")

# To compute the odds of successful tree recruitment as a function of mean diameter 
# and ombrophilous type

exp(coef(Recrut.binar))

## To create a 95% confidence interval for the estimate, type:

exp(confint.default(Recrut.binar))

# We see that the odds ratio corresponding to mean diameter is 0.87 (95% CI:  (0.80, 0.92)) 
# and Ombrophilous is 0.43 (95% CI: (0.23, 0.811))

# To use the fitted logistic regression curve to estimate the predictive power and 
# accuracy of the model. The most commonly used way to check the predictive ability of a 
# logistic model is to fit the ROC curve and its AUC (area under the curve), which are 
# typical performance measurements for a binary classifier. The ROC is a curve which plots 
# the true positive rate (TPR) against the false positive rate (FPR) at various threshold 
# settings, whereas the AUC is the area under the ROC curve. A commonly used rule to assess 
# the robustness of predictive power of the model, is to consider a good predictive ability 
# when AUC is closer to 1 (1 is ideal) than to 0.5.

# install.packages("ROCR")
library(ROCR)

p <- predict(Recrut.binar, type="response") # predict scorees
pr <- prediction(p, Lreg$Rec.binar)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)

auc <- performance(pr, measure = "auc")
auc <- auc@y.values[[1]]
auc

# or we can use the package plotROC for a better presentation of the AUC and true and 
# false positive rates.


# install.packages("plotROC")

library(plotROC)

plotROC(Lreg$Rec.binar, p)


# Here we got a reasonable predictive accuracy (0.64), but it could be better 
# (e.g. > 0.70).


#################### Moran's I: spatial autocorrelation test #########################

# Spatial autocorrelation refers to the correlation between the values of a variable 
# due their proximity in geographical space. This violates the assumption that they 
# are from independent observations. 

# The function "moransI" of the package "lctools" allows us to estimate the local 
# spatial correlation of the residuals from a linear model. 

# In order to execute the classical Moran's I test with this function, we first must 
# calculate and provide metric coordinates, such as UTM coordinates, for each row of 
# data. Firstly you must  convert the geographical coordinates belonging to your data 
# to UTM. After this, if you are using subplots (e.g. treatment with 50 subplots of 10 
# meters per 10 meters) you must add the distance of each subplot from the reference 
# point of the site to the UTM coordinates of the site. This is done separately for X 
# direction (East-West) and Y direction (North-South). If the reference point is the 
# South-West corner of the site, you add the distances of the subplot in W-E and S-N 
# direction to the site coordinates (e.g. 10 m and 30 m). But, if the reference point 
# is the North-East corner, you subtract the distances from the site coordinates.

# Each row of your data will contain real coordinates based on the sample units (here 
# represented by subplot) and you will be able to estimate spatial autocorrelation 
# in your linear model. 

# It is also important to note that Moran's I ranges from -1 to 1. Values close to -1 
# imply negative spatial autocorrelation (low values tend to have neighbours with high
# values and vice versa) and values close to 1 represent positive spatial autocorrelation
# (spatial clusters of similarly low or high values among neighbour subplots, as in the 
# present example). 

# The function "moransI" provides the Moran's I value (named as "Morans.I") and the two 
# tailed p-value (named "p.value.resampling" in the output) to test the significance of 
# Moran's I and, consequently, significance of spatial autocorrelation. 

# In order to graphically demonstrate the results for Moran's I, you can create a 
# correlogram with the function "correlog" of the package "ncf" as shown in the example 
# below. There, a univariate spatial correlogram based on the spatial distances between 
# the subplots is provided. Additionally, a plot with the corresponding p-values of the 
# Moran's I statistic is drawn, including a line of the threshold (p=0.05) to facilitate
# assessment of significance.



install.packages("lctools")
install.packages("ncf")

library(lctools) # Moran's I
library(ncf) # correlogram


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)


### LM to base the Moran's I computation

Mean.Diameter.pre <- lm(Mean_diameter ~ Mortality + Recruitment + Site, data=BIOVEG)

### Moran's I


Moran.M.diameter <- moransI(Coords=cbind(BIOVEG$local.utm.x, BIOVEG$local.utm.y), 
                   Bandwidth=4, resid(Mean.Diameter.pre), WType="Binary") 

Moran.M.diameter$Morans.I
Moran.M.diameter$p.value.resampling # significance of the Moran's I
Moran.M.diameter$p.value.randomization # another variant of p-value


## Graph: Correlogram

Cor.Mean_diameter <- correlog(x=BIOVEG$local.utm.x, y=BIOVEG$local.utm.y, 
                       as.vector(resid(Mean.Diameter.pre)), increment=15, resamp=1000)

x11(width=12, height=12) ####
par(mfrow=c(1,2))
plot(Cor.Mean_diameter$mean.of.class[1:6], Cor.Mean_diameter$correlation[1:6], 
  xlab="mean of distance class (m)", ylab="correlation (Moran's  I)", las=1)
abline(0, 0)
points(Cor.Mean_diameter$x.intercept, 0, pch="|")
legend("topright", "x intercept of autocorrelation", pch="|", bty="n", cex=0.6)
 plot(Cor.Mean_diameter$mean.of.class[1:6], Cor.Mean_diameter$p[1:6], ylim=c(0,1), 
  xlab="mean of distance class (m)", ylab="p-value of Moran's  I", las=1)
abline(0.05, 0)
dev.off() # turn off the graphics device



########################### LME (Linear Mixed Effects) ##############################


##### With spatial autocorrelation's correction included #####


# A LME model is a mixed model is implemented to deal with pseudoreplicates by 
# including random effects of, e.g., sites that meet the assumption of normal 
# distribution of residuals. In addition to random effects, it can handle spatial 
# autocorrelation. Therefore, if in a previous test (e.g. Moran's I) you find 
# significant spatial autocorrelation in a Linear Mixed Model (LMM; with random 
# and fixed effects and normal distribution), you must use LME with a spatial 
# correlation strucure in the model. If you have non-normal data with spatial 
# auto-correlation please check the function glmmPQL of the package MASS to 
# calculate a GLMM with a spatial correlation structure.

# There are different variants of spatial correlation structures. One good example 
# is Exponential correlation (function "corExp"). The correlation structure must be
# included in the syntax of the lme command as in the analyses below. If you don't 
# include a correlation function (e.g. corExp) in the LME model, it is just like a 
# common LMM. However, like GLMM, this would not be enough to resolve the 
# autocorrelation and provide reliable p-values. In order to test spatial 
# autocorrelation by Moran's I, you must first construct a LM with the same 
# predictors that will be used in the LME (see LM model used in the section 
# Linearity Tests and Moran's I). To better understand this approach of handling 
# spatial autocorrelation using the proper regression models, read Dormann et al. 
# (2007).

# Dormann, C.F., M. McPherson, J., B. Araújo, M., Bivand, R., Bolliger, J., Carl, 
# G.et al. (2007). Methods to account for spatial autocorrelation in the analysis 
# of species distributional data: A review. Ecography (Cop.)., 30, 609–628

# Note 1: LMM can be calculated using either Maximum Likelihood (ML) or Restricted 
# Maximum Likelihood (REML). If we are focusing on estimation and testing of fixed 
# effects, we choose "ML". In addition to correlation of the residuals, lme can also 
# deal with heterogeneous variance of residuals (heteroscedasticity), e.g. when 
# variances differ among groups, such as sites. For this purpose, we can use the 
# "weights" argument. 

install.packages("nlme")
install.packages("car")

library(nlme) # LME models
library(car) # Anova type II gives the predictors the same chance of 
             # being significant


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)


Mean_diameter.lme<- lme(Mean_diameter ~ Mortality + Recruitment, data=BIOVEG, 
              method="ML", random= ~1|Site, corr=corExp(form=~local.utm.x+ 
              local.utm.y|Site), weights=varIdent(form=~1|Site), na.action=na.omit)

summary(Mean_diameter.lme)

## Anova from car package

An<- Anova(Mean_diameter.lme, Type = II)
An

## Residuals

x11(width=12, height=12) ####
par(mfrow=c(2,2))
plot(fitted(Mean_diameter.lme), resid(Mean_diameter.lme))
lines(smooth.spline(fitted(Mean_diameter.lme), resid(Mean_diameter.lme), 
                                          spar=c(1)), col="red", lwd=2)
qqnorm(resid(Mean_diameter.lme))
qqline(resid(Mean_diameter.lme))
hist(resid(Mean_diameter.lme))
dev.off() ## turn off the graphics device

## Scatterplot 

b.0.Mean_diameter<- Mean_diameter.lme$coefficients$fixed[1]
b.1.Mean_diameter<- Mean_diameter.lme$coefficients$fixed[2]
b.2.Mean_diameter<- Mean_diameter.lme$coefficients$fixed[3]

x11(width=14, height=8)
par(mfrow=c(1,2))
plot(BIOVEG$Mean_diameter ~ BIOVEG$Mortality, las=1, ylab="Mean diameter", xlab="Mortality", 
     cex.lab=1.3)
curve(b.0.Mean_diameter + b.1.Mean_diameter*x, add=TRUE, col="red")
mtext(ifelse(An[1,3]<0.001, "p < 0.001", paste("p = ", round(An[1,3],3))))
plot(BIOVEG$Mean_diameter ~ BIOVEG$Recruitment, las=1, ylab="Mean diameter", xlab="Recruitment", 
     cex.lab=1.3)
curve(b.0.Mean_diameter + b.2.Mean_diameter*x, add=TRUE, col="red")
mtext(ifelse(An[2,3]<0.001, "p < 0.001", paste("p = ", round(An[2,3],3))))
dev.off()




######################### GLS (Generalized Least Squares) #############################


# GLS is another variant of linear regression. It is mostly used to perform regression
# analysis when there is significant spatial (or temporal) correlation between the 
# residuals (or if there is heteroscedasticity), but you do not have any random effect 
# in the model. In the absence of random effects, you cannot use lme, but you can use 
# gls. GLS models fit correlation structures of the variance-covariance matrix to model
# the (spatial) dependence of observations, just as LME models do. Thus, gls is an 
# efficient tool to deal with spatial autocorrelation when you would otherwise use lm 
# (read Dormann et al. 2007).

# Dormann et al. Methods to account for spatial autocorrelation in the analysis
# of species distributional data: a review. Ecography 30: 609628, 2007.

# N.B. GLS assumes a normal distribution.

# In cases such as in the example below, when you just want to know whether there 
# is variation of a variable among treatments and study sites, for  instance, but your 
# data comes from pseudoreplicates (e.g. subplots in our  example below), GLS is a good 
# option. Notice that in this case, you use site as a fixed effect. Otherwise, if you 
# want to include site as a random effect, you must use a mixed effect model 
# (GLMM or LME). Before calculating a GLS model, you should use Moran's I to check the 
# significance of spatial autocorrelation. If significant autocorrelation is found, 
# you must include a correlation structure in the model (e.g. correlation=corExp).

## Pre-test: Spatial autocorrelation (Moran'I)

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

## LM

Mean.Diameter.pre <- lm(Mean_diameter ~ Site, data=BIOVEG) 
summary(Mean.Diameter.pre)
shapiro.test(BIOVEG$Mean_diameter)

# Residuals

x11(width=12, height=12)
par(mfrow=c(2,2))
plot(Mean.Diameter.pre)
dev.off() ## turn off the graphics device
x11(width=12, height=12)
hist(resid(Mean.Diameter.pre),breaks=20) 
dev.off() ## turn off the graphics device

## Moran's I

library(lctools) # Moran's I

Moran.M.diameter <- moransI(Coords=cbind(BIOVEG$local.utm.x, BIOVEG$local.utm.y), 
                   Bandwidth=4, resid(Mean.Diameter.pre), WType="Binary") 

Moran.M.diameter$Morans.I
Moran.M.diameter$p.value.resampling # significance of the Moran's I
Moran.M.diameter$p.value.randomization



############# GLS ###################


# Since we found significant spatial autocorrelation in the pre-test with 
# Moran's I, a correlation (corExp) structure must be added in the syntax 
# of the model

library(nlme)
library(car)

Mean.diam.site<- gls(Mean_diameter ~ Site, correlation=corExp(form=~1|Site), 
                     data=BIOVEG)
 
summary(Mean.diam.site)
An<- Anova(Mean.diam.site)
An

# Residuals

x11(width=12, height=12)
par(mfrow=c(2,2))
plot(fitted(Mean.diam.site), resid(Mean.diam.site))
qqnorm(resid(Mean.diam.site))
hist(resid(Mean.diam.site))
dev.off()

# Boxplot

x11(width=12, height=12)
par(mfrow=c(1,1))
boxplot(BIOVEG$Mean_diameter ~ BIOVEG$Site, las=1, ylab="Mean Diameter")
mtext(ifelse(An[1,3]<0.001, "p < 0.001", paste("p = ", round(An[1,3],3))))
dev.off()


################################## QUADRATIC MODEL ####################################

#### Quadratic Hogweed model

# A quadratic regression is an extension of simple linear regression aiming to  
# find the best fit equation a set of data (dependent and predictor variables)  
# shaped like a parabola (e.g. U-shape curve). A quadratic fit is useful When 
# the coefficients and, consequently, the curve appear to fit the data better
# than the linear model does. 


Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

# Create linear model of species richness vs. Hogweed cover

M.lin <- glm(Species.richness ~ Hogweed.cover + Habitat.type, data=Hogweed, 
             family="quasipoisson")
summary(M.lin) 
anova(M.lin, test = "F")
# Hogweed cover has a significant negative linear effect on species richness

# Model diangostics

par(mfrow = c(2, 2))
plot(M.lin) 
# The diagnostic graphs look ok

# Plot the data and the regression line

par(mfrow = c(1, 1))
plot(log(Species.richness) ~ Hogweed.cover, data = Hogweed, 
     xlab = "Hogweed cover", ylab = "Species richness (log scale)")
abline(glm(Species.richness ~ Hogweed.cover + Habitat.type, data=Hogweed, 
           family="quasipoisson"), lwd = 2)
# The plot suggests that a quadartic model mignt fit the data even better

# Create variable of centered and squared Hogweed cover

Mean.cover <- mean(Hogweed$Hogweed.cover)
Hogweed.cover.sq <- (Hogweed$Hogweed.cover - Mean.cover)^2

# Calculate the model with original and, additionally, squared Hogweed cover

M.quad <- glm(Species.richness ~ Hogweed.cover + Hogweed.cover.sq + Habitat.type, 
              data=Hogweed, family="quasipoisson")
summary(M.quad)
anova(M.quad, test="F")

# Compare linear and quadratic model

anova(M.lin, M.quad, test = "F")
# Squared hogweed cover is (marginally) not significant.

# Anyway, let's create a scatter plot with regression curve

plot(Hogweed$Species.richness ~ Hogweed$Hogweed.cover)

b0 <- coef(M.quad)[1] # Intercept
b1 <- coef(M.quad)[2] # Hogweed cover linear
b2 <- coef(M.quad)[3] # Hogwee cover squared

# In the curve syntax, we need to center and square "x" before multiplying it with
# the coefficient of Hogweed.cover.sq.

curve(exp(b0 + b1*x + b2*(x-Mean.cover)^2), add=T, col="red", lwd=2)
	

######################## GNM (Generalized Non-Linear Models) #########################


# If you find significant non-linear relationships among variables, you can first  
# try transforming (e.g. log, square) the predictor variable. This will not always 
# be enough to resolve non-linearity issues and it can turn p-values of the linearity 
# tests non-significant. If non-linearity is not resolved using this procedure, models
#  such as LM, GLM and GLS are not the appropriate model for this purpose.  Indeed, 
# a better way to deal with non-linearity may be to apply a model that uses a  
# nonlinear function of the predictor variables, i.e. on the right-hand side of the 
# model formula. GNM models still require you to specify the distribution family of the 
# dependent variable (e.g. Gaussian, Poisson, Quasi-poisson). Furthermore, it is 
# possible to use an option to coerce a GLM model to a GNM model, as explained at the 
# end of this section.

# If the relationship between some of the predictor and the dependent variables  is 
# nonlinear, is necessary to use a proper model that accounts for the nonlinear effects.


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)


## Base LM model to apply linearity test

Abundance.pre <- lm(Abundance ~ Recruitment, data=BIOVEG) 
summary(Abundance.pre)


## Residuals

x11(width=12, height=12)
par(mfrow=c(2,2))
plot(Abundance.pre)
dev.off() ## turn off the graphics device
hist(resid(Abundance.pre),breaks=20) 
dev.off() ## turn off the graphics device
plot(BIOVEG$Site, BIOVEG$Abundance) # variation among sites


### Tests of linearity

library(car)

crPlots(Abundance.pre, col.lines=c("red", "green"))

# If one find some nonlinear relationship between predictor and
# dependent variable using the diagnostic for linearity showed
# previously here, thus a model that account for nonlinear terms
# is necessary to carry the statistical prediction.



############# GNM model #############

install.packages("gnm")
library(gnm)


Exp <- function(expression, inst = NULL){list(predictors = list(substitute
(expression)), term = function(predictors, ...) {paste("exp(", predictors, ")", 
                   sep = "")}, call = as.expression(match.call()), match = 1)} 
class(Exp) <- "nonlin"


N.gnm<- gnm(Abundance ~ Exp(Recruitment), data=BIOVEG, method = "gnmFit", 
start = NULL, checkLinear = TRUE, verbose = FALSE, family="quasipoisson")

summary(N.gnm)
An <- anova(N.gnm, test = "F")
An

## Residuals

x11(width=12, height=12)
par(mfrow = c(2,2), oma = c(0, 0, 3, 0))
title <- paste(deparse(N.gnm$formula, width.cutoff = 50), collapse = "\n")
plot(N.gnm, sub.caption = title) 
par(mfrow = c(1,1))
hist(resid(N.gnm))

# OR

par(mfrow=c(2,2))
plot(fitted(N.gnm), resid(N.gnm))
qqnorm(resid(N.gnm))
qqline(resid(N.gnm))
hist(resid(N.gnm))
dev.off()


## Deviance profile (optional) - just similar to a summary

prof<- profile(N.gnm, which = ofInterest(N.gnm), alpha = 0.05, trace = TRUE)
prof
x11(width=12, height=8)
plot(prof)
dev.off()


## scatterplot 


x11(width=12, height=12)
plot(BIOVEG$Abundance ~ BIOVEG$Recruitment, ylab="Abundance", 
xlab= "Annual Recruitment", las = 1, cex.lab=1.2)
points(BIOVEG$Recruitment, fitted(N.gnm), col = "red", lwd = 2)
Smooth <- smooth.spline(BIOVEG$Recruitment, fitted(N.gnm))
lines(Smooth, col = "red", lwd = 2)
mtext(ifelse(An[2,6]<0.001, "p < 0.001", paste("p = ", round(An[2,6],3))))
dev.off()


# Note 1: If you need to coerce objects of class "glm" to an object of class 
# "gnm" follow the procedure below


Nall.glm<- glm(Abundance ~ Recruitment, data=BIOVEG, family="quasipoisson")
N.gnm2<- update(asGnm(Nall.glm))
N.gnm2
anova(N.gnm2, test = "F")
summary(N.gnm2)


######################## GAM (Generalized Aditive Models) ##############################

# For a GAM model example, let's assume that we only have fixed effects data(e.g. 
# vegetation type). We will use the same variables as used for the model in a pre- 
# test for linearity ,using the function "crPlots" of the package "car". Since 
# "P_mg_100g"  and "K_mg_100g" demonstrated the greatest deviation from linearity, 
# we apply a smoothing term (s) for them in the model.

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

install.packages("mgcv")
library(mgcv) # carry out GAM model

spp.gam<- gam(log(Species.richness) ~ s(P_mg_100g) + s(K_mg_100g) + Veg.cover, 
                                              data=Hogweed, family="gaussian")

summary(spp.gam)

# In the summary, we can see two types of output: parametric coefficients, i.e., 
# the slope of the predictor that indicated a linear relationship and was not 
# subjected to smoothing and the non-linear predictors which were smoothed, showing 
# their F and p-values for significance. The ‘edf’ refers to the estimated degrees of 
# freedom and basically, the larger the number, the wigglier the fitted model.  Values 
# close to 1 indicate a linear term.

anova(spp.gam) # the output is quite similar to the summary

# Below we calculate a measure for correlation (analogous to collinearity) known as
# concurvity, between the predictors in the model, to consider non-linear correlations. 
# The range is from 0 to 1. The closer the value is to 1, the higher the correlation. 
# To know more about concurvity read:

# Marra, G. & Wood, S. N. Practical variable selection for generalized additive models. 
# Comput. Stat. Data Anal. 55, 2372–2387 (2011).

concurvity(spp.gam, full=FALSE)

# Residuals

par(mfrow = c(2,2))
gam.check(spp.gam) 

# Here the graphical residuals look generally acceptable, since the histogram is not
# that nice (i.e. symmetrical). 

# Together with the graphical residuals, the model returns a written output that shows 
# parameters for each predictor with its smoothing term e.g. "k", "edf", "K-index" and 
# its p-value. "k-index" is the ratio of neighbor differencing scale estimate to fitted 
# model scale estimate; "edf" is the previously mentioned estimated degrees of freedom 
# and "k'" is the maximum possible EDF for the term. These parameters test whether the 
# basis dimension for smoothing is adequate for the predictor. 

# For further details regarding this kind of residual diagnostic, please consult the 
# documentation of the function "gam.check" of the package "mgcv". One way to correctly 
# apply a smoothing term for the predictor is to choose the most appropriate and insert 
# it into the model using the argument "bs". Look for "smooth.terms" in the mgcv 
# documentation.

# In order to test for improvement of the fit of the model, we included the type of 
# smoothing with the argument "bs" and used "K" to choose the number of knots (k) 
# once cubic regression splines have a set number of knots. The default of knots 
# automatically defined by the model “spp.gam” was 9. Let’s test the model fit with 
# 6 and 9 knots. Here we apply  "cc" (cyclic cubic regression splines) smoothing to 
# "bs" because of the observed tendency of non-linear fitting in a pre-diangostic 
# graph with the function "crPLots" of package "car". Here we test of redefine the   
# knots for the predictor P_mg_100g because it was diagnosed as being the most  
# non-linear, and improving its fit also may improve the overall fit of the model.

par(mfrow = c(2, 2))
spp.gam2<- gam(log(Species.richness) ~ s(P_mg_100g, bs='cc', k=6) + 
               s(K_mg_100g) + Veg.cover, data=Hogweed, family="gaussian")
plot(spp.gam2)
spp.gam3<- gam(log(Species.richness) ~ s(P_mg_100g, bs='cc', k=9) + 
                 s(K_mg_100g) + Veg.cover, data=Hogweed, family="gaussian")
plot(spp.gam3)

# We can see that setting models with 6 and 9 knots provides similar model fitting.
# Thus, the first model (spp.gam) with 9 knots defined automatically by the model
# is already sufficient for a good model fit. 

# If we still call the residuals for each model test (6 and 9 knots), 
# the degrees of freedom (EDF) of the smoothers are close especially 
# for the predictor "K_mg_100g" and they k-index is nonsignificant.


par(mfrow = c(2,2))
gam.check(spp.gam2)

par(mfrow = c(2,2))
gam.check(spp.gam3) 


##################### GAMM (Generalized Aditive Mixed Models) ###########################

# To simulate a GAMM model we take the dependent and  almost all predictor variables as 
# we used in the previous GAM example. Our intention here is to also include random 
# effects in our modeling, rather than only fixed effects as we did in GAM analysis.

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

library(mgcv) # carry out GAMM model

spp.gamm<- gamm(log(Species.richness) ~ s(Hogweed.cover) + s(P_mg_100g), random = 
                          list(Study.area = ~ 1), data=Hogweed, family="gaussian")

summary(spp.gamm)

# Similar to the summary of the GAM model, we have two types of output: parametric 
# coefficients  i.e., the slope of the predictor that indicated a linear relationship 
# and was not smoothed and the non-linear predictors with smoothing and their F and 
# p-values for significance.

anova(spp.gamm$gam)

# Below we again conduct a measure for by concurvity between the predictors in the 
# model, to consider non-linear correlations.

concurvity(spp.gamm$gam, full=FALSE)

# Residuals

par(mfrow = c(2,2))
gam.check(spp.gamm$gam) 

# Here the graphical residuals looks fine in general. However, the predictor
# "Hogweed.cover" has significant low p-value for k. 
# In order to remember: "k-index" is the ratio of the neighbor differencing scale 
# estimate to fitted model scale estimate; "edf" is the previously mentioned 
# estimated degrees of freedom and "k" is the maximum possible EDF for the term. 
# These parameters test whether the basis dimension for a smoothing is adequate 
# for the predictor.Therefore, here you can try to improve the model by setting
# different number of knots (k) and maybe choosing another smooth term for this
# predictor (e.g. "gp" (Gaussian process), "tp" (Thin plate regression splines)
# and others. Please, check this out searching for "smooth.terms" documentation
# of the package mgcv.



################################## SUPERVISED LEARNING ##################################

############# Random Forest #############


# The Random Forest model is a supervised learning technique that generates multiple 
# models from a training dataset on the given data and then simply combines their output 
# rules (predictive variables), thus generating a robust high performance model that 
# corrects overfitting and balances variance inequalities. That is, this model improves 
# the predictive power from decision trees, reducing their variances by calculating their 
# averages. The decision tree is a type of modeling that operates with information gain 
# on each node, classifying data points with greater information increment on each node. 
# When all nodes are depleted for their information gain ratings, the model achieves its 
# optimal performance result. Thus, the Random Forest model is considered, in many cases, 
# to be the most robust decision tree approach.

# In its results, the Random Forest provides a percentage explanation of the Y variance 
# from the whole set of predictors considered.  It also provides the mean square residual 
# error in the predictions from an approach called “Out of Bag Error Estimation (OBB)”, 
# a robust and efficient method of fitting and error estimation.


Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
Hogweed$Plot <- as.integer(Hogweed$Plot) # Transform predictor Plot that is factor
                                        # in integer, since random Forest cannot 
                                        # handle more the 53 categories
head(Hogweed)


# install.packages("randomForest")

library(randomForest) 

# Below we split our data in train (70% of data) and test (30% of data) sets.

set.seed(100)
train <- sample(nrow(Hogweed), 0.7*nrow(Hogweed), replace = FALSE)
train <- sort(train)
TrainSet <- Hogweed[train,] # split the data in train set
TestSet <- Hogweed[-train,] # split the data in test set
head(TrainSet)# visualize if all variables are present in train set
head(TestSet)# visualize if all variables are present in test set
summary(TrainSet)# visualize if all data for train set is correct, i.e. 
# Min, Median, Mean, Max...
summary(TestSet) # visualize if all data for test set is correct


# Verify the best number of predictors (mtry) per node (split) in the branches 
# of the decision trees. The function "mtry" shows from which number of predictors 
# per node the prediction errors are lower, based on OBB error (Out of bag) estimation.

tuneRF(Hogweed, Hogweed$Species.richness) # check how many predictors per node reduce
                                          # error. The outcome shows that increasing 
                                          # the number of predictors per node of the trees
                                          # decreases the error. Thus, we include all
                                          # predictors (16) per node in the model. 
                                          

#### Random Forest Model #####

set.seed(100)
RF1 <- randomForest(Species.richness ~ ., data=Hogweed,
                    subset=train, mtry = 16, importance=TRUE, ntree=500)



print(RF1) # The model with train subset performed poorly regarding the percentage of 
           # variance explained in the prediction of species richness by the selected
           # set of predictors. 

# Let's see if the number of decision (500) chosen for the model were enough to reduce 
# the errors:

plot(RF1) # The output shows that 500 trees are sufficient to reduce error, since from
          #  approximately 100 trees, the error starts to decrease.

### Importance of predictors ###

# Here we check which predictors are the least and most important when predicting the 
# variation of the dependent variable. The variable importance in Random Forest is 
# computed based on Node Purity and Mean squared error. Node purity measures the total 
# decrease of impurities in the nodes (splits) by computing the average of each predictor 
# over the nodes of all trees,  being measured, in regression, by the sum of the squares 
# residuals. Of these two approaches for importance, Mean squared error (%IncMSE) is the 
# most robust and informative measure, where the higher number, the more important the 
# variable. Conversely, IncNodePurity is biased and should only be used if the extra 
# computation time of calculating %IncMSE is challenging. Since the %IncMSE is calculated 
# based on OBB estimation, the expectation is that the MSE will increase, especially if 
# the variable has some importance, thus the label "% IncMSE".

import<- importance(RF1)
import  


# The outcome states that the most important predictors based on %IncMSE. Remember that 
# the higher number, the more important the variable.


### save output of importance 

write.table(import, "Importance.csv", row.names=T, sep=",", dec=".") 


## Note: Prior to selecting the predictor variables for the Random Forest model, verify
# whether there is high correlation among them. The presence of highly correlated 
# predictors may bias the computation of importance of each one, for predicting the 
# dependent variable, leading to suboptimal predictor variables being artificially 
# and arbitrarily preferred. For more information about this topic, please read the 
# recommended papers below:

# Strobl C, Boulesteix AL, Zeileis A, Hothorn T.Bias in random forest variable 
# importance measures: illustrations, sources and a solution. BMC Bioinformatics, 
# 2007 Jan 25;8:25

# Strobl C et al. Conditional variable importance for random forests. 
# BMC Bioinformatics, BMC Bioinformatics. 2008 Jul 11;9:307

# Shows and save importance graphically

library(lattice)
tiff(filename="Fig_IMPORTANCE2019.jpg", width=180, height=80, units="mm", res=300)

varImpPlot(RF1, sort=TRUE, main=deparse(substitute(train)),  
           cex.lab=1.5, las=1, bty="n", cex=0.5,pch=19, col='black') 

dev.off() ## turn off the graphics device 


###### Comparison of the level of error of the model in relation to the number of 
# predictors per node for the train (70%) and test (30%) set ######

## The OBB (out of bag error estimation) is computed to the train set and compared 
# to the test set (test.err).

### Here we simulate whether the inclusion of all predictors per node of the trees 
# delivers a better or worse performance for the random forest prediction of both
# train and test set.
# Then, we reduce the number of trees from 500 to 400 to better fit  the whole
# number of predictors per node.


oob.err=double(16) # 16 regards to total number of predictors
test.err=double(16)

for(mtry in 1:16) 
{
  rf=randomForest(Species.richness ~ . , data = Hogweed , subset = train,mtry=mtry,ntree=400) 
  oob.err[mtry] = rf$mse[400] # Error of all trees adjusted
  
  pred<-predict(rf,Hogweed[-train,]) # Prediction in test set for each tree
  # Mean Squared Test Error:
  test.err[mtry]= with(Hogweed[-train,], mean( (Species.richness - pred)^2))
  
  cat(mtry," ") # show the outcome in R console
  
}

test.err
oob.err


matplot(1:mtry , cbind(oob.err,test.err), pch=19 , col=c("red","blue"),type="b",
        ylab="Mean Squared Error",xlab="Number of Predictors Considered at each Split")
legend("center",legend=c("Out of Bag Error","Test Error"),pch=19, col=c("red","blue"))


# The graph shows that the mean squared error increases slightly in both Out of bag 
# error and test error as the number of predictors increases at each split (node). 

# In addition, we can check the performance using cross-validation prediction
# of models with the number of predictors reducing sequentially. In such 
# reduction, the predictors are ranked by variable importance via a nested 
# cross-validation procedure. 


rfcv(Hogweed, Hogweed$Species.richness, cv.fold=5) 

# The section "$error.cv" in the output show us that reducing number of 
# predictors decreases the rates of error in the performance of the model. 
# This outcome make sense and is congruent with large amount of predictors
# with low importance for the model prediction observed above.


#### Prediction performance ####

# In order to verify the percentage of explanation for the Random Forest model we can 
# adjust a linear model to get a value for R-squared. We do this for the test  set created 
# above. If 30% of data from the test produces a high performance, your prediction is good.

predTest <- predict(RF1, TestSet)
r2_test <- lm(predTest ~ TestSet$Species.richness)
summary(r2_test) # The adjusted R-squared shows that the prediction performance of the 
                 # model is not good, since it predicts less than 50 % of the dependent 
                 # variable.


# Nevertheless, you do not need necessarily to use a linear model to check on prediction
# performance. You could also calculate the correlation coefficient and, then, square it:

cor(predTest, TestSet$Species.richness) ^ 2


##### RMSE #####

# We can adjust a RMSE (Root mean squared error) for the model. The RMSE is particularly 
# useful when it is compared with other models to see if it performs better or worse 
# depending if it has a higher or lower error of prediction.

rfvalpred = predict(RF1,newdata=TestSet)
rmse = sqrt(mean((TestSet$Species.richness-rfvalpred)^2))
rmse


############################### Model Selection (Akaike) ###############################

# For both Akaike and GLM / GLMM the predictor variables must be tested for the 
# degree of correlation, because if these variables are highly correlated, it 
# implies collinearity, and the real contributions of each predictor to the dependent 
# variable may be confounded. One conventional way to test correlation, and collinearity, 
# is the Pearson's correlation (described in this script). Pearson correlation 
# coefficients can be between -1 and 1, and the accepted limit for collinearity is |0.60|. 
# Correlations above 0.60 indicate too much collinearity between predictors. On the other 
# hand, many authors accept correlations up to  0.70 and use variables with this degree 
# of correlation in the model.


# Akaike's function: select the best regression model according to the relationship 
# between the dependent and predictor variables. The best model is normally assumed 
# to be the one with the lowest AICc value. However, the quantity of models considered 
# to be good does not necessarily follow a specific rule, since it depends on your 
# purpose and research approach. It also can depend on the statistical line you follow. 
# Some researchers suggest that a difference between AIC values smaller than 2 (i.e. 
# delta AIC =< 2) indicates that the models are almost equally good. Others use a 
# threshold level of delta AIC =< 4.


### GLM's Akaike selection ####

## Akaike: basic syntax

## install and load "MuMIn" package

# Usage

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

install.packages("MuMIn")
library(MuMIn)
options(na.action = "na.fail")
div <- glm(Veg.height ~ P_mg_100g + K_mg_100g + Land.use, family=gaussian, 
           data=Hogweed)
tested.div1 <- dredge(div)
tested.div1 


### LME's Akaike selection ###

# Firstly you must carry out the pre-tests with LM in the model to be executed as 
# LME model. If normality and linearity are accepted, then, you check for spatial 
# autocorrelation by Moran's I. If there is significant autocorrelation you must
# model include a spatial correlation structure (e.g. corExp) within the global  
# LME that will be submited to selection by Akaike such as the example below.


# Example of LM model to check for linearity and autocorrelation for the global
# LME model to execute Akaike's selection:

# Firstly you must conduct pre-tests with LM’s in the model to be executed as LME 
# model. If normality and linearity are acceptable,  you must then check for spatial
# autocorrelation by Moran's I. If there is significant autocorrelation your model 
# must include a spatial correlation structure (e.g. corExp) within the global LME, 
# which will be submitted to selection by Akaike criteria, as described in the 
# example below.

Mean.Diameter.pre <- lm(Mean_diameter ~ Mortality + Recruitment + Richness + Site, 
                                                                      data=BIOVEG) 

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)


## LME's Akaike selection

library(nlme) # LME models
library(MuMIn) # Dredge function to execute the selection

options(na.action = "na.fail")
div <- lme(Mean_diameter ~ Mortality + Recruitment + Richness, data=BIOVEG, method=  
              "ML", random= ~1|Site, corr=corExp(form=~local.utm.x+local.utm.y|Site), 
                                  weights=varIdent(form=~1|Site), na.action=na.fail)


tested.div1 <- dredge(div)
tested.div1



####### QAIC - QUASI Akaike ########

# When the generalized model shows oversdispersion, the typical AIC value cannot be 
# computed using the dredge function. Instead, the QAIC value needs to be computed. 
# With the function "QAIC" (Quasi AIC) of the package MuMIn, calculated a modification 
# of Akaike’s Information Criterion for overdispersed count data (e.g. quasipoisson 
# model) can be calculated, but it also may also be used to compute QAIC for binomial 
# data (e.g. quasipoisson).


####### QAIC - Quasipoisson ########

Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

# install.packages("MuMIn")
library(MuMIn)
options(na.action = "na.fail")


# Fist, we execute the function for "x.quasipossion" constructor, which allows
# for quasipoisson family objects, for ML estimation. 


x.quasipoisson <- function(...) { 
  res <- quasipoisson(...) 
  res$aic <- poisson(...)$aic 
  res 
}

riqueza.glm2<- glm(Species.richness ~ P_mg_100g + Land.use + K_mg_100g 
                + Hogweed.height, data=Hogweed, family="x.quasipoisson")

# To compute the variance inflation factor based on residuals of Pearson, where it is 
# necessary to use the "dredge" function to compute the QAIC values and its parameters.
# There is a discussion among statisticians regarding whether residuals of Pearson are 
# suitable for models (e.g. Poisson and Binomial) based on Chisq, but there are various 
# examples where it provides a  robust and reliable fit.

chat<-sum(residuals(riqueza.glm2,"pearson")^2)/riqueza.glm2$df.residual 
chat 

riqueza.QAIC <- dredge(riqueza.glm2, rank = "QAIC", chat = chat)
riqueza.QAIC

# We can compare the QAIC selection by using Pearson's residuals to fit variance 
# inflation factor,  by dividing the deviance by residual degrees of freedom.

chat2=(deviance(riqueza.glm2) / df.residual(riqueza.glm2))

riqueza2.QAIC <- dredge(riqueza.glm2, rank = "QAIC", chat = chat2)
riqueza2.QAIC

# We see that the outcomes of QAIC selection using the two approaches above, to 
# compute 'chat' as the variance inflation factors, results in virtually the 
# same outcomes.


####### QAIC - Quasibinomial ########

# install.packages("MuMIn")
library(MuMIn)


Hogweed <- read.table("Hogweed.csv", header=T, sep=",", dec=".")
head(Hogweed)

#  First, let's execute a 'hacked' constructor for quasibinomial family 
# object for ML estimation

x.quasibinomial  <- function(...) {
  res <- quasibinomial(...)
  res$aic <- binomial(...)$aic
  res
}



coverHog.glm1 <- glm(Hogweed.cover/100 ~ Veg.height + N_perc + 
                       Habitat.type + Land.use + P_mg_100g, 
                     data=Hogweed,family="x.quasibinomial")

chat3=(deviance(coverHog.glm1) / df.residual(coverHog.glm1))

options(na.action = "na.fail")
AICselectCover <- dredge(coverHog.glm1)
AICselectCover


################################ MODEL AVERAGING ######################################

# With the function "model.avg" of the package MuMIn we can compute model averaging
# based on an information criterion (e.g. AIC (Akaike) or BIC (Baysean)).


BIOVEG <- read.table("BIOVEG3.csv", header=TRUE, sep=",", dec=".")
head(BIOVEG)

library(MuMIn)
options(na.action = "na.fail")
riqueza <- glm(Richness ~ Mortality + Recruitment + Diameter 
               + Temperature_mean, family=gaussian, data=BIOVEG)
riqueza.AICc <- dredge(riqueza)
riqueza.AICc 

# Average of the best modes under delta < 4

summary(model.avg(riqueza.AICc, subset = delta < 4))



#################### POST-HOC TEST TO REGRESSION MODELS ################################

# A common necessity when using GLM, LM or GLMM using categorical (factor variable) 
# predictors, is to conduct multiple comparison of means among the pairs of categories. 
# Such comparison using regression models is still being developed and improved. For 
# instance, there are currently no 100% robust, standard unbiased post-hoc tests which 
# can be applied in Mixed Models (e.g. GLMM, LME). 
# Below we provide examples of two of the most commonly used post-hoc methods applied 
# in Ecology but we highly advise you to read  about their application and carefully 
# consider whether their use is appropriate for your analysis. 

# If after you run your model, you find a significant difference of the means among 
# categories (as shown by fitting ANOVA to count data in this GLM model example with 
# species richness as target variable), then you need a post-hoc (posteriori) test to 
# show you between which categories the difference occurred. Post-hoc tests on Linear 
# Models become even more important the more predictor variables you use.  If you 
# want to run an  a posteriori comparison of the model, knowing  that there are 
# significant differences among categories as there are in this example, this 
# comparison will consider the weight of other predictor variables (e.g. k, ca, light). 



### TUKEY'S POST HOC TEST (Tukey's HSD (honest significant difference) test) ###

BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

# Usage

install.packages("multcomp")
library(multcomp)


richtype <- glm(Richness ~ type, data=BIOVEG, family=poisson)
summary(richtype)
summary(glht(richtype, linfct=mcp(type = "Tukey")))
cldm0 <- cld(glht(richtype, linfct=mcp(type = "Tukey")))
cldm0 # here show between each pairs of type there is the significant difference


### LMERTEST ###

# This test is used for mixed models (e.g. GLMM, LME) which have Pseudoreplication 
# (when your data is distributed in subplots within each site. E.g. 3 sites with 
# 50 10 m x 10 m plots) and you need to consider this effect (normally called random 
# effect, but be careful when reading about this subject since there remains some 
# disagreement between statisticians). For instance, if you need to not only consider 
# data values per site or fragment, but also data values per subplot within each site. 
# Suppose you need to know the phylogenetic diversity per subplot within each site, 
# you need to recognize that the subplots are not real replicates and so you need  
# to use a mixed model. 

# Atention: The "lmer test" from the package of the same name is one of the most 
# commonly used post-hoc tests in this case, but there remain disagreements among  
# statisticians about its robustness when homogeneity of variance and normality 
# are violated. Read about this method and carefully consider whether it is appropriate 
# for your analysis. The "lmerTest" package calculate p-values of differences of means 
# among categories.

install.packages("lmerTest")
library(lmerTest) 


BIOVEG <- read.table("BIOVEG3.csv", header=T, sep=",", dec=".")
head(BIOVEG)

## use to Continuos data

# call the function "glmer" from package "lme4"              
DIAM<-lmer(Diameter ~ type + (1|Site), REML=FALSE, data=BIOVEG)

difflsmeans(DIAM, test.effs="type")
 


############################### MULTIVARIATE ANALYSIS #################################

# This approach is indicated to analyze more than one statistical outcome variable  
# at a time. Usually, multivariate analysis to establish relationships among variables,
# where they were sampled by multiple measurements in sample or experimental units, 
# and their sampling design.

####################### NMDS - Non-metric multidimensional scaling ####################


# NMDS is used to robustly demonstrate the original position of data in 
# multidimensional space using a reduced number of dimensions that can be 
# easily plotted and visualized.

# use "vegan" package:

install.packages("vegan")
library("vegan") 

# Adapted NMDS template option:

matrix<- read.table("matrixspp.txt", header=T, sep="\t")
matrix

# NMDS

NMDS_bray<- metaMDS(matrix, k=2, distance="bray")
NMDS_bray$stress ## < 0.3:  check if  it was lower like this in your analysis

# If stress < 0.3 the program returns “TRUE” which indicates acceptable 
# levels of stress within your model. It implies how robust your analysis' 
# result was.

plot(NMDS_bray) ### ordinary graph, just distinction in 2 colors and symbols

# Save scores

write.table(NMDS_bray$points, "sitesR.txt", row.names=T, sep="\t")
write.table(NMDS_bray$species, "species.txt", row.names=F, sep="\t")


## WARNING: The data sets used here are for training purposes only.  The 
# two matrices from test dataset don't represent a realistic  species 
# composition of the vegetation types described here.  We collated species 
# data from our personal datasets and shuffled them into the columns of the 
# species matrix. For example, you mind find that species common in Semideciduous
# forest (SEMI) might instead be listed amongst Rain Forest (OMB) species. 

# Better NMDS graph: distinction between vegetation types (as in the following 
# example: ssm, as, fo). You must to use 2 matrix files, one with species and 
# another one with classification in vegetation types (e.g. if you just want to 
# check species distribution between types) or environmental variables (e.g. soil:
# ph, Al...). In this example we just want to check the species distribution among 
# vegetation types. Thus, you need the species matrix file (here it is the object 
# "matrix") and another one with category (here it is "matrix4"). Each row number 
# of "matrix4" corresponds to the same row number in "matrix". Therefore both matrices
# must have the same sequence and number of rows. See both matrices in the example 
# below:

# matrix                   

# abarema_sp casearia_arborea 
# 0           1
# 0           0   
# 1           2
# 0           1 

# matrix 4 

# OMB         fo
# OMB
# ALUV        as

# matrix 5 -> example of soil matrix

# 0.4       pH
# 0.76      Al.
# 0.25      Ca


# Note 1: again all of the these 3 matrices must have the same sequence and 
# number of rows. Within each row you classify the soil variable value and 
# element, or in case of vegetation the category.


### Better Graph (separation among types: OMB, GAL, SEMI)

matrix<- read.table("matrixspp.txt", header=T, sep="\t")
matrix

matrix5<- read.table("matrixtype.txt", header=T, sep="\t")
matrix5

dune.nmds <- metaMDS(matrix, distance="bray")
grupo<- matrix5$codigos
grupo.nivel <- levels(grupo)
grupo
grupo.nivel
escores.nmds <- scores(dune.nmds)
escores.nmds
dune.nmds$stress 


p <- ordiplot(escores.nmds, type="n", main="Vegetation Types")
grupo.nivel
points(escores.nmds[grupo=="OMB",], pch=(18), cex=2, col="black")
points(escores.nmds[grupo=="GAL",], pch=(19), cex=2, col="darkgrey")      
points(escores.nmds[grupo=="SEMI",], pch=(20), cex=2, col="red") 



###################### Analysis of similarities (ANOSIM) ############################

# Basic template by "vegan"

# The two matrices used here are the same as those used previously in NMDS analysis

# usage:           

library("vegan")

# Data

matrix<- read.table("matrixspp.txt", header=T, sep="\t")
matrix

matrix5<- read.table("matrixtype.txt", header=T, sep="\t")
matrix5

# Analysis-ANOSIM

matrix
matrix5
attach(matrix5)
dune.dist <- vegdist(matrix, "bray")
dune.ano <- anosim(matrix, matrix5$codigo)
summary(dune.ano)
plot(dune.ano) 
dune.nmds$stress 

# The graph described in NMDS analysis can also be used to illustrate the 
# result of  ANOSIM.           
           

# Important: This ANOSIM using R just provides the Global R statistics and p-value, 
# i.e. whether there are differences or not. If you need to know between which pairs 
# there are significant differences, the most pratical and acessible way is to do 
# ANOSIM in the "PAST" software. To do this you insert the names of the categories 
# in the first column (grey), which corresponds to the rows , and the species' names
# in the horizontal grey row, which corresponds to the columns. Your columns will 
# contain the species' names. In the rest of the cells are filled with the species 
# abundance data per plot. The last step is to mark the first column with type names
# in different colors. Just go in edit and you'll find how to insert different colors
# to specific lines. "Past" requires at least 2 different color groups to be defined.
# Lastly, just click on "Multivariatve" and select "ANOSIM-One way". The output will 
# return the Global and the pairwise values.



##################### CCA (Canonical Correspondence Analysis) ########################


# CCA, also refered as constrained correspondence analysis, is one of the most 
# commonly used multivariate analyses to model and to test effects of  environmental 
# variables (e.g. soil properties (ph, N, P, K etc) on species  composition of plant 
# communities.

install.packages("vegan") 
library(vegan)

# Usage:

# Species matrix
matrix<- read.table("matrixspp.txt", header=T, sep="\t")
matrix[1:5,1:4]

# Environmental matrix
env<- read.table("matrixenv.txt", header=T, sep="\t")
head(env)

forest.cca <- cca(matrix ~ pH +	P + K + Ca + Mg + Al, env, na.action = na.fail)
forest.cca
anova(forest.cca, by="margin")
summary(forest.cca) # scores of CCA
# ordinary plot
x11(width=12, height=12)
plot(forest.cca, las=1)
dev.off()


######################## PCA (Principal Component Analysis ##########################

# PCA is a method commonly used to separate variables (e.g environmental data, such 
# as soil conditions) in ordination before submitting them to further analysis. 

# To execute the PCA in the example below, we use the function rda, though other 
# functions such as prcomp or princomp can also be used to calculate PCA. It's 
# important to correctly interpret the PCA output before trying to understand values 
# returned by the function "PCAsignificance" and the command "PC1.exp", which are 
# explained below. The command PC1.exp shows the percentage of explanation of each 
# axis. In the function "PCAsignificance" the result named "percentage of variance"
# shows the percentage of variation explained by axis 1. You can access further 
# information about this, and other values, by requesting R help for this function 
# by writing: ?PCAsignificance.

# The output of rda, v, contains the scores of the "species" (here: soil variables) 
# and sites on the ordination axes. These are not the same as the correlations of 
# variables with the axes.

# pca.env$CA$v    # Relationship of each environmental variable with the axes.         
                  # The closer the result is to ±1, the larger the association 
                  # with the positive or negative axis, respectively.



install.packages("vegan") 
install.packages("BiodiversityR") 

library(vegan) # rda function for PCA analysis
library(BiodiversityR) # function PCAsignificance

# Environmental matrix
env<- read.table("matrixenv.txt", header=T, sep="\t")
head(env)

# Usage:

pca.env <- rda(env, scale=TRUE) # scale = TRUE - to standardize the data
pca.env
summary(pca.env) 
PCAsignificance(pca.env, axes=6) 


(av <- pca.env$CA$eig) # av - inspects the eigenvalues. The sum of eigenvalues 
                       # represent the total percentage of explanation assigned 
                       # to each axis

sum.av <- sum(av)
PC1.exp <- 100*(av[1]/sum(av))  # calculate the percentage of explanation of 
PC1.exp                         # each axis. Here it shows the percentage of the 
                                # of the first axis (PCA 1).

calculates the percentage of explanation of 
PC1.exp

round(PC1.exp, 2) # You can choose how many decimals after the point you would 
                  # like to keep. Here it is 2.


# Simple gragh

x11(width=12, height=12)
plot(pca.env, type="n")
points(pca.env, display = "sites", col=1, pch=21, cex=1, bg="grey")

# Biplot graph

par(mfrow=c(1,2))
plot1 <- biplot(pca.env, scaling=1, type="text", las=1) # sites scaled by 
                                                        # eigenvalues
plot2 <- biplot(pca.env, scaling=2, type="text", las=1) # species scaled by 
                                                        # eigenvalues
dev.off()


# Screeplots - select axes with values from the ordination that are larger 
# than the expected value under the Broken stick criteria.

screeplot(pca.env, bstick = TRUE, type = "lines")




################### Analysis of Indicator Species (INDVAL) ###########################


# The association index named IndVal computed with the function "multipatt" of the 
# package "indicspecies" is generated by combinations of the input clusters and 
# comparison of each combination with the matrix species. Thus, each species receives 
# the combination with the highest association value and the best matching patterns 
# are tested for statistical significance of the associations. Thus this analysis is 
# useful to indicate species that can be used to gather information about specific 
# characteristics of habitats or environment (e.g. forest types; soil types; drought 
# or wetter sites).

# The name "type" within the command below refers to the categorical variable in the 
# forest type matrix (matrixtype.txt). "codigos" are labels given to abbreviate  
# forest type names: OMB = Rain (Ombrophilous) Forest; GAL = Gallery Forest and   
# SEMI = Semideciduos forest.

install.packages("indicspecies")

library(indicspecies)

# Species matrix
matrixspp<- read.table("matrixspp.txt", header=T, sep="\t")
matrixspp[1:5,1:4]

# Vegetation type
type <- read.table("matrixtype.txt", header=T)
head(type)

wetpt = multipatt(matrixspp, type$codigos) # Runs the combination analysis  
                                  # using IndVal.g as statistic
wetpt

summary(wetpt) # Lists those species with significant association to one combination



############# Permutational Analysis of variance (PERMANOVA) ########################

# PERMANOVA is a non-parametric multivariate statistical test used to compare groups 
# of objects based on any measure of distance. Commonly this analysis is used to 
# compare taxon in samples where the groups of samples are compared.

library(vegan)
library(permute)
library(lattice)

spp <- read.table("matrixspp.txt", header=T)
spp[1:5,10:13]

log(spp+1)->spp2
spp2[1:5,10:13]


# Usage: 

# Testing for species composition over different forest types: the name "type" 
# within the command below is the categorical variable in the forest type matrix 
# (matrixtype.txt). "codigos" are just labels given to abbreviate forest typenames:
# OMB = Rain (Ombrophilous) Forest; GAL = Gallery Forest and 
# SEMI = Semideciduos forest.


type<- read.table("matrixtype.txt", header=T, sep="\t")
head(type)

adonis(spp2 ~ codigos, data=type, permutations=999) -> perma.result1
perma.result1

 
# Testing for soil influence in the tree species abundance: the names in the right  
# side of the formula regard the soil variables  contained in the environmental 
# matrix ("matrixenv.txt"). The category "types" regard the forest type (OMB, GAL 
# and SEMI).

env <- read.table("matrixenv.txt", header=T)
head(env)

adonis(spp2 ~ pH + P + K + Ca + Mg + Al, data=env, permutations=999)-> perm.result
perm.result

## Result to permanova pseudo-F, R2 e p.
#?adonis


############################### RECOMMENDED READING #################################


##### Highlight:


# GOTELLI, N.J.; ELLISON, A.M. A Primer of Ecological Statistics. Sinauer Associates 
# Publishers, 2004.

# Version in Portuguese: # GOTELLI, N.J.; ELLISON, A.M. Princípios de Estatística em 
# Ecologia. Porto Alegre: Artmed, 2011. 527 p.

# CRAWLEY, M.J. The R Book: Second Edtion. Chichester: Imperial College London
# at Silwood Park, John Wiley & Sons Ltd., 2012. 1051p.

# BORCARD, D.; GILLET, F.; LEGENDRE, P. Numerical ecology with R. New York: 
# Springer, 2011. 306 p.

# FORTIN, M. J.; DALE, M. Spatial analysis: a guide for ecologists. Cambridge:
# Cambridge University, 2005. 365 p.

# HINKLE, D.E.; WIERSMA, W. JURS, S.G. Applied Statistics for the Behavioral 
# Sciences. 5th ed. Boston: Houghton Mifflin, 2003. 756 p.

# LEGENDRE, P.; LEGENDRE, L. Numerical ecology. Third English edition. 
# Amsterdam: Elsevier Science, 2012. 1006 p.
                                  
###### Complementary:

# BEGON'S text book is highly recommendable reading on sampling design and statistical
# analysis in ecological studies (e.g. distribution, population growth, behaviour) of
# plants and animals.

# BEGON, M.; TOWNSEND, C.R.;  HARPER, J.L. Ecology: From Individuals to Ecosystems, 
# 4th Edition.  Wiley-Blackwell, 2005. 750 pages.

# Version in Portuguese: BEGON, M.; TOWNSEND, C.R.;  HARPER, J.L. Ecologia: de 
# indivíduos a ecossistemas. 4ª ed. Porto Alegre: Artmed, 2007. 752 p.

#####################################################################################