################################################################################
#'@author James Margrove
#'@title seedling recruitment analysis 

Sys.setenv(LANG="en")

# Clear workspace 
rm(list=ls())

# Import packages 
require(MuMIn)
require(ggplot2)
require(MASS)
require(car)
require(arm)
require(doSNOW)

# Import data 
data <- read.table("./data/SeedlingData.txt", header = T)
data <- subset(data, DBH > 50)
data <- droplevels(data)

# sBA50 that are unknown make 0
data$sBA50[which(is.na(data$sBA50))] <- 0

data$pCF <- arm::invlogit(data$sAC)
data$pHF <- arm::invlogit(data$AC)
data <- subset(data, pCF > 0.4)

formula5 <- as.formula(seedlings ~ (sBA50 * I(LGI) + I(LGI^2))^3) # all 
model1 <- glm.nb(formula5, data, na.action = "na.pass")
summary(model1)
require(doSNOW)
clust <- makeCluster(8, "SOCK")
clusterExport(clust, c("data","glm.nb","model1"))
d1 <- pdredge(model1, trace = 2, rank = "AIC", cluster = clust)
model2 <- get.models(d1, 1)[[1]]
AIC(model2)
summary(model2)
Anova(model2)
par(mfrow=c(2,2))
plot(model2)
par(mfrow=c(1,1))

pred1 <- expand.grid(sBA50 = round(quantile(data$sBA50, 
                                            c(0.025,0.25,0.5,0.75,0.975)),0),
                     LGI = seq(min(data$LGI),max(data$LGI),length=50))

pred1$seedlings <- predict(model2, pred1)
pred1$se <- predict(model2, pred1, se.fit = T)$se.fit
pred1$sBA50 <- factor(pred1$sBA50)
head(pred1)

ggplot(pred1, aes(x=LGI, y = exp(seedlings), fill = sBA50))  + geom_line() + 
  geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.333) + 
  geom_line() + theme_classic() + xlab("Conspecific flowering density") + ylab("Seedlings 16m^2") + scale_alpha(guide = 'none') + theme(legend.position = c(0.15,0.8))

