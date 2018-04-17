###'
###'
###'    With the LGI and TH LHF as seperated from the focal trees flowering intensity Full script of the analysis 
###

rm(list=ls())
### Reed in the data frame 
data <- read.table("data.txt", header = T)
str(data)
### abundnace as a frequency 
n.indv <- with(subset(data, DBH > 50), tapply(sp,sp,length))
n.dt <- with(data, tapply(sp,sp,length))
n.indv[which(is.na(n.indv))] <- 0
data$ABN <- rep(n.indv, n.dt)

######################################################################################################################
### Creating a table of the important information for teh methods 
AB <- with(data, tapply(sBA, fullname, mean))
sp <- levels(data$fullname)
minLM <- with(data,tapply(mlim, fullname, mean))
pF <- with(data, tapply(FL2010, fullname, mean, na.rm= T))*100
ABN <- with(data, tapply(seedlings, fullname, length))
msed <- with(data[is.na(data$FL),], tapply(seedlings, fullname, mean, na.rm = T))
nqd <- with(data[is.na(data$FL),], tapply(seedlings, fullname, function(x){length(x[!is.na(x)])}))
FH_msed <- with(data[!is.na(data$FL),], tapply(seedlings, fullname, mean, na.rm = T))
FH_nqd <- with(data[!is.na(data$FL),], tapply(seedlings, fullname, function(x){length(x[!is.na(x)])}))
dtable <- data.frame(sp, minLM, ABN, AB, pF, nqd, msed, FH_nqd, FH_msed)
#write.table(dtable, file = "Table of species.txt")

#######################################################################################################################
### Select on the trees with seedlings, LGI's  
data <- data[!is.na(data$seedlings),] # all data with seedlings
data <- data[which(!is.na(data$GI)),] # There are 4 tree species where we did not know where there NFC was no NFC in the plot - these were dropped 
data[which(is.na(data$sBA)),"sBA"] <- 0 # There is one species that has no individuals above the 50 cm limit Vmic
data <- subset(data, DBH >= 30)
### LOG GI 
data$LGI <- log(data$GI_II)
data$LHF <- log(data$HF_II)
### Select data without flowering history 
fdata <- data[is.na(data$FL),] # to train the data set 
fdata <- droplevels(fdata)
dim(fdata)
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
### Using MASS glm.nb calculate the expdential term for the gravity analysis 
require(MASS)
require(ggplot2)
### LOG GI 
plot(data$LHF, data$LGI) # not colinear 

### creating the full model 
model <- glm.nb(seedlings ~ (sBA + log(DBH) + LGI + I(LGI^2) + LHF)^3 , fdata, na.action = "na.pass")
summary(model)
require(MuMIn)
require(doSNOW)
#clust <- makeCluster(8, "SOCK")
#clusterExport(clust, c("model","fdata","glm.nb"))
#d2 <- pdredge(model, rank = "AIC", trace = 3, cluster = clust)
#save(d2, file = "D2 DBH sep.R")
load(file = "D2 DBH sep.R")
stopCluster(clust)
head(d2, 20)
dim(fdata)
#save(d1, file = "D1 on con,ba,hetero.R")
model2 <- get.models(d2, df == 5)[[1]] ### Retrive the best model. best model also has a weight of 0.56 which is quite high 
summary(model2)
AIC(model2)
4275.61
### The model with 7 df has a higher AIC by 2.1, compared to a model with 11 df

sum(with(fdata,tapply(sp,sp,length)))
length(levels(fdata$sp))
### Predicting the models results 
pred1 <- expand.grid(sBA = round(quantile(fdata$sBA, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1),
                     LGI = seq(min(fdata$LGI),max(fdata$LGI),length=50))

pred1$seedlings <- predict(model2, pred1)
pred1$se <- predict(model2, pred1, se.fit = T)$se.fit * 1.96
pred1$Abundance <- factor(pred1$sBA)
### Plotting the results 
require(scales)
ggplot(pred1,aes(x=LGI, exp(seedlings), fill = Abundance))  + geom_line()  + 
  geom_line() + geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75) + 
  geom_line() + theme_classic() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  scale_fill_manual(values =c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')) + 
  theme(legend.position = c(0.15,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))  + 
  scale_x_continuous(breaks=pretty_breaks(6))

#######################################################################################################################
#########################################################################################################################
### informaiton on the model for the paper 
### what was the differance in AIC for a model where the quadratic terms was not pressent 
summary(model2)
modeli <- update(model2, .~. - I(LGI^2):sBA  - I(LGI^2))
summary(modeli)
diff(AIC(model2, modeli)[,2])
### What was the peek recruitment for the most abundant species 
predi <- data.frame(sBA = max(data$sBA), 
                    LGI = seq(min(data$LGI), max(data$LGI), length = 100))

predi$seedlings <- predict(model2, predi)
predi$se <- predict(model2, predi, se.fit = T)$se.fit
predi[which(predi$seedlings==max(predi$seedlings)),]

exp(3.186822)          
exp(3.186822 + 0.3528094)
exp(3.186822 - 0.3528094)

### What does the conspecific flowering density correspond too
mDBH_ptom <- mean(data[data$sp == "Ptom",]$DBH, na.rm = T)
sqrt(1/((exp(4.33)/5)/(mDBH_ptom^2)))
#######################################################################################################################
#######################################################################################################################
# Try the same model with abundance as a frequency 
model5 <- update(model2, .~. - I(LGI^2):sBA - sBA + I(LGI^2):ABN + ABN)
summary(model5)
AIC(model5, model2)
diff(AIC(model5, model2)[,2])

summary(lm(ABN~sBA, fdata))
with(fdata, cor(ABN,sBA))



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#########################################################################################################################
### prepare the History data frame from analysis 
hist <- data[!is.na(data$FL),]
hist <- droplevels(hist)
str(hist)

### 
require(INLA)
### Setting up the formula for INLA 
head(hist)
mat <- with(hist, model.matrix(seedlings ~ (FL + LGI + FL + I(LGI^2) + sBA)^3))
head(mat)
hist$int1 <- mat[,6]
hist$int2 <- mat[,7]
hist$int3 <- mat[,8]
hist$int4 <- mat[,9]
hist$int5 <- mat[,12]
hist$int6 <- mat[,13]
hist$int7 <- mat[,14]

summary(model2)
# formula with the priors 
formula <- seedlings ~ FL + 
  f(I(LGI^2), model = "linear", mean.linear = -0.0265434, prec.linear = 1/( 0.0093357^2)) + 
  f(sBA, model = "linear", mean.linear = 0.0003095, prec.linear = 1/(0.0009597^2)) + 
  f(I(sBA*LGI^2), model = "linear", mean.linear = -0.0010615, prec.linear = 1/(0.0001774^2)) 


# Running the model 
model <- inla(formula, data = hist, family = "nbinomial",
              control.fixed = list(mean.intercept=3.1995690,prec.intercept = 1/(0.1268740^2)), 
              control.compute = list(cpo=T, dic = T))

-2*sum(log(model$cpo$cpo))

summary(model)
model$dic$dic
# compared to the null model Fl reduced the DIC by 5.5573 
#
#FL = 987.2463
#FL + int1  =   989.2
#FL + int2  =   989.2
#FL + int3  =   980.77
#FL + int1 + int2 =  991.0
#FL + int1 + int3 =  976.0233
#FL + int2 + int3 = 978.37
#FL + int1 + int2 + int3 =  977.8
#
# FLPrior:I(LGI^2) + FLPrior:sBA these two terms interact with the flowering history 

# the inclusion of interaction terms flowering intensity or species abundance did not reduced the AIC - 4 points per parameter

# Pit is relatively flat (probability intergral transformation)
hist(model$cpo$pit)

####
#### I'd like to predict the effects from the model 
####


pred2 <- expand.grid(FL = levels(hist$FL), LGI = mean(data$LGI), sBA = mean(data$sBA), seedlings = NA)
hist2 <- hist[,c("FL","LGI","sBA","seedlings")]

model.pred <- inla(formula, data = rbind(hist2, pred2) , family = "nbinomial",
                   control.fixed = list(mean.intercept=3.1995690,prec.intercept = 1/(0.1268740^2)), 
                   control.compute = list(cpo=T, dic = T, config = T))

### Create a funciton to bootstrap the Credible intervals 
require(foreach)
FH_boot <- function(){
  s = inla.posterior.sample(1, result=model.pred)
  s = s[[1]]$latent
  tail(s, 7)[1:2]}

FH_boots <- foreach(i = 1:1000, .combine = rbind) %do% FH_boot()
#save(FH_boots, file = "boots_mean_data_n1000.R")
CI <- apply(FH_boots, 2, quantile, c(0.025,0.975))

### making the data frame to plot the [predicted values and confidence intervals
pred3 <- data.frame(FL = levels(hist$FL), CI025 = CI[1,], CI975 = CI[2,])
pred3$seedlings <- model.pred$summary.fitted.values$mean[147:148]

ggplot(pred3,aes(x=FL, exp(seedlings)))  + geom_point(size = 3) +
  geom_errorbar(aes(ymin=exp(CI025),ymax=exp(CI975)), width = 0.1) + 
  theme_classic() + xlab("Flowering history") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + 
  theme(legend.position = c(0.8,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))

1-(exp(1.741848)/exp(2.692223)) # drop of 61.3%

# so at mean flowering intensities and mean abundances for thethe whole data set, we find that there is a diffe 44.7 reducution
summary(model)
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#### predicing the curvs for the inla model at Once flowered 

pred3 <- expand.grid(FL = levels(hist$FL), 
                     LGI = seq(min(data$LGI),max(data$LGI), length = 100),
                     sBA = round(quantile(data$sBA, c(0.1,0.25,0.5,0.75,0.9)),1), 
                     seedlings = NA)
head(pred3)
dim(pred3)
hist3 <- hist[,c("FL","LGI","sBA","seedlings")]
head(hist3)
model.pred2 <- inla(formula, data = rbind(hist3, pred3) , family = "nbinomial",
                    control.fixed = list(mean.intercept=3.1995690,prec.intercept = 1/(0.1268740^2)), 
                    control.compute = list(cpo=T, dic = T, config = T))
summary(model.pred2)
### Create a funciton to bootstrap the Credible intervals 
require(foreach)
FH_boot <- function(){
  s = inla.posterior.sample(1, result=model.pred2)
  s = s[[1]]$latent
  s[147:1146]}

FH_boots2 <- foreach(i = 1:100, .combine = rbind) %do% FH_boot()
#save(FH_boots, file = "boots_mean_fdata_int_n5000.R")
CI <- apply(FH_boots2, 2, quantile, c(0.025,0.975))

### making the data frame to plot the [predicted values and confidence intervals
pred3$CI025 <- CI[1,] 
pred3$CI975 <- CI[2,]
pred3$seedlings <- model.pred2$summary.fitted.values$mean[147:1146]
pred3$Abundance <- factor(pred3$sBA)
head(pred3)

require(scales)
ggplot(pred3,aes(x=LGI, exp(seedlings), fill = Abundance))  + geom_line() + facet_wrap(~FL) + 
  geom_line() + geom_ribbon(aes(ymin=exp(CI025),ymax=exp(CI975)), alpha = 0.75) + 
  geom_line() + theme_classic() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  scale_fill_manual(values =c('#9ecae1','#6baed6','#4292c6','#2171b5','#084594')) + 
  theme(legend.position = c(0.15,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))  + 
  scale_x_continuous(breaks=pretty_breaks(6))

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#### Now just the parashorea tom

ptom <- subset(hist, sp == "Ptom")
ptom <- droplevels(ptom)
str(ptom)

### modelling the data 
model3 <- glm.nb(seedlings ~ FL_inst, ptom, na.action = "na.pass")
summary(model3)
anova(model3)

d5 <- dredge(model3, trace = 2, rank = "AIC")
head(d5)
model6 <- get.models(d5, 1)[[1]]
summary(model6)
AIC(model6)

plot(model6, which = 2)
plot(model6, which = 1)

### predicting from the model
pred4 <- expand.grid(FL_inst = levels(ptom$FL_inst), LGI = mean(ptom$LGI))
pred4$seedlings <- predict(model3, pred4)
pred4$se <- predict(model3, pred4, se.fit = TRUE)$se.fit


ggplot(pred4,aes(x=reorder(FL_inst, - seedlings), y = exp(seedlings))) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=exp(seedlings-se), ymax=exp(seedlings+se)), width = 0.2) + 
  theme_classic() + xlab("Flowering history") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + 
  theme(legend.position = c(0.8,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))

exp(3.590439)/exp(2.431418)
### predicting from the model
pred4 <- expand.grid(FL_inst = levels(ptom$FL_inst)[1], LGI = seq(min(ptom$LGI), max(ptom$LGI), length = 100))
pred4$seedlings <- predict(model3, pred4)
pred4$se <- predict(model3, pred4, se.fit = TRUE)$se.fit


ggplot(pred4,aes(x=LGI, y = exp(seedlings))) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin=exp(seedlings-se), ymax=exp(seedlings+se)), width = 0.2) + 
  theme_classic() + xlab("Flowering history") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + 
  theme(legend.position = c(0.8,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))

timesF <- as.character(ptom$FL_inst)
timesF[ptom$FL_inst == "Twice07"] <- "Twice"
timesF[ptom$FL_inst == "Twice09"] <- "Twice"
ptom$timesF <- timesF


model7 <- glm.nb(seedlings ~ FL + timesF + FL_inst, ptom, na.action = "na.pass")
summary(model7)
anova(model7)
plot(model7, which = 2)
0.1818 /59.120*100
