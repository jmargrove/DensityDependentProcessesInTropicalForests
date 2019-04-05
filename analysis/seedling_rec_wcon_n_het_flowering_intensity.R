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
str(data)
# Calculate the abundance 
n.indv <- with(subset(data, DBH > 50), tapply(sp,sp,length))
n.dt <- with(data, tapply(sp,sp,length))
n.indv[which(is.na(n.indv))] <- 0
data$ABN <- rep(n.indv, n.dt) / 160 

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
#d_abn <- pdredge(model, rank = "AIC", trace = 3, subset = dc(LGI, I(LGI^2)), cluster = clust)

head(d_abn)

?pdredge

#save(d_abn2, file = "./dredged/glmnb_2Abundance.R")
#stopCluster(clust)

#load the dredge 
#load(file = "./dredged/glmnb_Abundance.R")
#model2 <- get.models(d_abn, 1)[[1]]
#summary(model2)

model2 <- glm.nb(seedlings ~ ABN + I(LGI^2) +  LGI + LHF + log(DBH) +
                   ABN:I(LGI^2) + ABN:LHF + I(LGI^2):LHF + 1, data)

summary(model2)
anova(model2)
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
                                          c(0.1, 0.5, 0.9)),1),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = round(quantile(data$LHF, 
                                          c(0.1, 0.9)),1),
                     LGI = seq(min(data$LGI),max(data$LGI),length=50))

pred1$seedlings <- predict(model2, pred1)
dfr <- model2$df.residual


Sys.setenv(LANG = "en")

pred1$se <- predict(model2, pred1, se.fit = T)$se.fit * qt(0.95, dfr)
pred1$`Community wide density` <- factor(pred1$ABN)
head(pred1)
pred1$LHF[which(pred1$LHF == 0.9)] <- "Low heterspecific flowering intensity (0.9)"
pred1$LHF[which(pred1$LHF == 2.2)] <- "High heterspecific flowering intensity (2.2)"
pred1$LHF <- as.factor(pred1$LHF)
pred1$LHF <- relevel(pred1$LHF, ref = "Low heterspecific flowering intensity (0.9)")
?relevel
#dev.off() # problems with plotting - use dev.off: close a ploting device 
cols <- c('#B3001B', '#018E42', '#053C5E', 'yellow')

# 
#install.packages("viridis")
#require(viridis)


str(pred1)

pred1

# Plotting the results 
p2 <- ggplot(pred1, aes(x=LGI, exp(seedlings), linetype = `Community wide density`, fill = `Community wide density`)) + geom_line()  + 
  geom_line() + geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75) + 
  geom_line() + 
  theme_light() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  facet_grid(~factor(LHF), scales = "free_y") + 
  theme(legend.position = "top") + 
  scale_fill_viridis(discrete = T) + 
  #scale_fill_manual(values = cols) + 
  scale_x_continuous(breaks=pretty_breaks(6))

p2

ggsave(p2, file = "./graphs/het_flower_ints_Abundance.png", 
       width = 8, 
       height = 4.2)


# Numbers from the paper
# Predicting the models results 
data[which(data$ABN == 0 ), "sp"]
with(data, tapply(ABN, sp, mean))


pred2 <- expand.grid(ABN = c(0.09375, 3.8812),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = quantile(data$LHF, c(0.025, 0.975)),
                     LGI = quantile(data$LGI, c(0.025, 0.975)))

pred2$seedlings <- predict(model2, pred2, type = "response")

with(data, tapply(ABN, sp, mean, na.rm = T))
# Calculate the numbers for the species 
df <- (with(pred2, tapply(seedlings, list(LHF, ABN), max)))


(1 - df[1,3]/df[1,1]) * 100# rare species 


41/6.5

(1 - df[2,3]/df[1,2]) * 100# abundent species 

# Which species are these at the
data[which(data$ABN == 621), "sp"]

## maybe a desnity plot 


ggplot(pred1, aes(x = Abundance, y = LGI, fill = exp(seedlings))) + geom_tile()
head(pred2)





# Predicting the models results 
pred2 <- expand.grid(ABN = seq(11, 600, length = 50),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = round(quantile(data$LHF, 
                                          c(0.025,0.975)),1),
                     LGI = round(quantile(data$LGI, 
                                          c(0.025,0.975)),1))

pred2$seedlings <- predict(model2, pred2)
dfr <- model2$df.residual

pred2$se <- predict(model2, pred2, se.fit = T)$se.fit * qt(0.95, dfr)
pred2$Abundance <- (pred2$ABN)
head(pred2)
#dev.off() # problems with plotting - use dev.off: close a ploting device 
cols <- c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')

# Plotting the results 
p2 <- ggplot(pred2, aes(x=(Abundance), exp(seedlings), fill = factor(LHF))) +
  facet_grid(~factor(LGI)) + 
  geom_line() + 
  geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75, width = 30) + 
  theme_classic() + xlab("Abundance") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
#  scale_fill_manual(values = cols) + 
  theme(legend.position = "top")# + 
  #scale_x_continuous(breaks=pretty_breaks(6))


p2



ggsave(p2, file = "./graphs/het_flower_ints_AbundanceII.png", 
       width = 6, 
       height = 4)





# Predicting the models results 
pred3 <- expand.grid(ABN = seq(11, 600, length = 2),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = seq(min(data$LHF), max(data$LHF), length = 100),
                     LGI = round(quantile(data$LGI, 
                                          c(0.025,0.975)),1))

pred3$seedlings <- predict(model2, pred3)
dfr <- model2$df.residual

pred3$se <- predict(model2, pred3, se.fit = T)$se.fit * qt(0.95, dfr)
pred3$Abundance <- (pred3$ABN)
head(pred3)
#dev.off() # problems with plotting - use dev.off: close a ploting device 
cols <- c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')

# Plotting the results 
p3 <- ggplot(pred3, aes(x=LHF, exp(seedlings), fill = factor(Abundance))) +
  facet_grid(~factor(LGI)) + 

  geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.75, width = 30) + 
  theme_classic() + xlab("Heterospecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  geom_line() + 
  #  scale_fill_manual(values = cols) + 
  theme(legend.position = "top")# + 
#scale_x_continuous(breaks=pretty_breaks(6))
p3



ggsave(p3, file = "./graphs/het_flower_ints_AbundanceIII.png", 
       width = 6, 
       height = 4)


