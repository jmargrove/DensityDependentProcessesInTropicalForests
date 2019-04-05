################################################################################
#' @author James Margrove 
#' @title Fitting the spline models

# Clear workspace 
rm(list=ls())

# Import packages 
require(splines)

# Import data 
data <- read.table("./data/data.txt", header = TRUE)
str(data)

# Calculate the abundance 
n.indv <- with(subset(data, DBH > 50), tapply(sp,sp,length))
n.dt <- with(data, tapply(sp,sp,length))
n.indv[which(is.na(n.indv))] <- 0
data$ABN <- rep(n.indv, n.dt)

data <- data[!is.na(data$seedlings),] # All data with seedlings
data <- data[which(!is.na(data$GI)),] # There are 4 tree species where we did not know where there NFC was no NFC in the plot - these were dropped 
data[which(is.na(data$sBA)),"sBA"] <- 0 # There is one species that has no individuals above the 50 cm limit Vmic
data <- subset(data, DBH >= 30)

# Log the flowering intensities 
data$LGI <- log(data$GI_II)
data$LHF <- log(data$HF_II)

# Import model 
model3 <- glm.nb(formula = seedlings ~ ABN * ns(LGI, df = 2) + 1, data = data)
summary(model3)
load("./models/bestOverallModel.R")
summary(model2)
model4 <- update(model2, . ~ . - I(LGI^2):LHF - log(DBH) - ABN:LHF - LHF)
diff(AIC(model4, model3)[,2])
summary(model4)

# Predicting the models results 
pred1 <- expand.grid(ABN = round(quantile(data$ABN, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = round(quantile(data$LHF, 
                                          c(0.1,0.25,0.5,0.75,0.9)),1)[3],
                     LGI = seq(min(data$LGI),max(data$LGI),length=50))

pred1$seedlings <- predict(model3, pred1)
dfr <- model2$df.residual

pred1$se <- predict(model3, pred1, se.fit = T)$se.fit * qt(0.95, dfr)
pred1$Abundance <- factor(pred1$ABN)
head(pred1)
#dev.off() # problems with plotting - use dev.off: close a ploting device 
cols <- c('#a1d99b','#74c476','#41ab5d','#238b45','#005a32')

# Plotting the results 
p2 <- ggplot(pred1, aes(x=LGI, exp(seedlings), fill = Abundance)) + geom_line()  + 
  geom_line() + geom_ribbon(aes(ymin=exp(seedlings-se),ymax=exp(seedlings+se)), alpha = 0.5) + 
  geom_line() + theme_bw() + xlab("Conspecific flowering intensity") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + scale_alpha(guide = 'none')  + 
  scale_fill_manual(values = cols) + 
  theme(legend.position = "top") + 
  scale_x_continuous(breaks=pretty_breaks(6))

p2


ggsave(p2, file = "./graphs/splines_LCFI.png", 
       width = 4, height = 4)
