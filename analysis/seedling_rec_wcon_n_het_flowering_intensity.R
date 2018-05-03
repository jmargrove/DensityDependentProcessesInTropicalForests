################################################################################
#' @title seedling recruitment analysis Abundance
#' @author James MArgrove

# Clear workspace 
rm(list=ls())

# Import packages 
require(MASS)
require(ggplot2)
require(MuMIn)
require(doSNOW)
require(scales)
require(RColorBrewer)

# Import data 
data <- read.table("./data/data.txt", header = T)

# Calculate the abundance 
n.indv <- with(subset(data, DBH > 50), tapply(sp,sp,length))
n.dt <- with(data, tapply(sp,sp,length))
n.indv[which(is.na(n.indv))] <- 0
data$ABN <- rep(n.indv, n.dt)

data <- data[!is.na(data$seedlings),] # All data with seedlings
data <- data[which(!is.na(data$GI)),] # There are 4 tree species where we did not know where there NFC was no NFC in the plot - these were dropped 
data[which(is.na(data$sBA)),"sBA"] <- 0 # There is one species that has no individuals above the 50 cm limit Vmic
data <- subset(data, DBH >= 30)

dim(data)
# Log the flowering intensities 
data$LGI <- log(data$GI_II)
data$LHF <- log(data$HF_II)
# Plot the flowering intensities to ensure that they are not colinear 
plot(data$LHF, data$LGI) # not colinear 
model <- glm.nb(seedlings ~ (ABN + DBH + LGI + I(LGI^2) + I(LHF^2))^3 , data, na.action = "na.pass")
summary(model)

plot(model, which = 1)
plot(model, which = 2)
?glm.nb

#clust <- makeCluster(8, "SOCK")
#clusterExport(clust, c("model","data","glm.nb"))
#d_abn <- pdredge(model, rank = "AIC", trace = 3, cluster = clust)
#save(d_abn2, file = "./dredged/glmnb_2Abundance.R")
#stopCluster(clust)

#load the dredge 
load(file = "./dredged/glmnb_Abundance.R")
model2 <- get.models(d_abn, 1)[[1]]
summary(model2)
AIC(model2)
plot(model2, which = 1)

# Remove the quadratic term and replace with linear term for AIC diff calculation 
model3 <- update(model2, .~. - I(LGI^2) + LGI - 
                   LHF:I(LGI^2) + LHF:LGI - 
                   ABN:I(LGI^2) + ABN:LGI)
summary(model3)
LGIdif <- AIC(model2, model3)
diff(LGIdif[, "AIC"])

# Calculate the log(DBH) to DBH value 
model4 <- update(model2, .~. - log(DBH) + DBH)


# Predicting the models results 
pred1 <- expand.grid(ABN = round(quantile(data$ABN, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = round(quantile(data$LHF, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1),
                     LGI = seq(min(data$LGI),max(data$LGI),length=50))

pred1$seedlings <- predict(model2, pred1)
dfr <- model2$df.residual

pred1$se <- predict(model2, pred1, se.fit = T)$se.fit * qt(0.95, dfr)
pred1$Abundance <- factor(pred1$ABN)
head(pred1)
#dev.off() # problems with plotting - use dev.off: close a ploting device 
cols <- c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')

# Plotting the results 
p2 <- ggplot(pred1, aes(x=LGI, exp(seedlings), fill = Abundance)) + geom_line()  + 
  geom_line() + geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75) + 
  geom_line() + theme_bw() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  scale_fill_manual(values = cols) + 
  facet_grid(~factor(LHF)) + 
  theme(legend.position = "top") + 
  scale_x_continuous(breaks=pretty_breaks(6))


p2

ggsave(p2, file = "./graphs/het_flower_ints_Abundance.png", 
       width = 10, 
       height = 4)


# Numbers from the paper
# Predicting the models results 
data[which(data$ABN == 0 ), "sp"]
with(data, tapply(ABN, sp, mean))
pred2 <- expand.grid(ABN = quantile(data$ABN, c(0.025, 0.975)),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = c(min(data$LHF), max(data$LHF)),
                     LGI = c(min(data$LHF), max(data$LHF)))

pred2$seedlings <- predict(model2, pred2)

dfr <- model2$df.residual
pred2$se <- predict(model2, pred1, se.fit = T)$se.fit * qt(0.95, dfr)

pred2$Abundance <- factor(pred2$ABN)
#


# Calculate the numbers for the species 
df <- exp(with(pred1, tapply(seedlings, list(LHF, Abundance), max)))
(1 - df[5,1]/df[1,1]) * 100# rare species 
(1 - df[5,5]/df[1,5]) * 100# abundent species 

# Which species are these at the
data[which(data$ABN == 621), "sp"]
