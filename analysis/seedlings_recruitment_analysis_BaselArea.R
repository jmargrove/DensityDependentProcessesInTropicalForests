################################################################################
#' @title seedling recruitment analysis Basel Area 
#' @author James MArgrove

# Clear workspace 
rm(list=ls())

# Import packages 
require(MASS)
require(ggplot2)
require(MuMIn)
require(doSNOW)
require(scales)

# Import data 
data <- read.table("./data/data.txt", header = T)

data <- data[!is.na(data$seedlings),] # All data with seedlings
data <- data[which(!is.na(data$GI)),] # There are 4 tree species where we did not know where there NFC was no NFC in the plot - these were dropped 
data[which(is.na(data$sBA)),"sBA"] <- 0 # There is one species that has no individuals above the 50 cm limit Vmic
data <- subset(data, DBH >= 30)

# Log the flowering intensities 
data$LGI <- log(data$GI_II)
data$LHF <- log(data$HF_II)
# Plot the flowering intensities to ensure that they are not colinear 
plot(data$LHF, data$LGI) # not colinear 

# Creating the full model 
model <- glm.nb(seedlings ~ (sBA + log(DBH) + LGI + I(LGI^2) + LHF)^3 , data, na.action = "na.pass")
summary(model)

  #clust <- makeCluster(8, "SOCK")
  #clusterExport(clust, c("model","data","glm.nb"))
  #d2 <- pdredge(model, rank = "AIC", trace = 3, cluster = clust)
  #save(d2, file = "./dredged/glmnb_baselArea.R")
  #stopCluster(clust)

#load the dredge 
load(file = "D2 DBH sep.R")

head(d2, 20)
dim(fdata)
#save(d1, file = "D1 on con,ba,hetero.R")
model2 <- get.models(d2, df == 5)[[1]]
summary(model2)
AIC(model2)

# Predicting the models results 
pred1 <- expand.grid(sBA = round(quantile(fdata$sBA, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1),
                     LGI = seq(min(fdata$LGI),max(fdata$LGI),length=50))

pred1$seedlings <- predict(model2, pred1)
pred1$se <- predict(model2, pred1, se.fit = T)$se.fit * 1.96
pred1$`Basel Area` <- factor(pred1$sBA)

# Plotting the results 
ggplot(pred1,aes(x=LGI, exp(seedlings), fill = `Basel Area`))  + geom_line()  + 
  geom_line() + geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75) + 
  geom_line() + theme_classic() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  scale_fill_manual(values =c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')) + 
  theme(legend.position = c(0.15,0.8), 
        axis.text = element_text(size = 15),
        text = element_text(size=16))  + 
  scale_x_continuous(breaks=pretty_breaks(6))
