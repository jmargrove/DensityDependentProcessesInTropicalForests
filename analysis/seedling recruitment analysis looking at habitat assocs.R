################################################################################
#' @title seedling recruitment analysis looking at habitat assocs 
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
data <- droplevels(data)
dim(data)
# Log the flowering intensities 
data$LGI <- log(data$GI_II)
data$LHF <- log(data$HF_II)

# OK 
load(file = "./dredged/glmnb_Abundance.R")
model2 <- get.models(d_abn, 1)[[1]]
summary(model2)

sp_median <- with(data, tapply(elev, sp, median, na.rm = T))
sp_length <- with(data, tapply(elev, sp, length))
data$sp_median <- rep(sp_median, sp_length)
data$dif_elev <- with(data, sp_median - elev)
# check for colinarly between the elev and sp_median
with(data, plot(sp_median, elev))
with(data, plot(sp_median, dif_elev^2))
with(data, plot(x = elev, dif_elev^2))
with(data, plot(x = elev, dif_elev))
with(data, plot(x = sp_median, dif_elev))

quantile(data$dif_elev, c(0.025, 0.975))


# ok update the models with the hypothesis 
# 1.  that there are less at lower elevations 
# 2. that there is a quadratice shap about the median value


model3 <- update(model2, . ~ . + elev +  I(dif_elev^2))
summary(model3)


# Predicting the models results quadratic 
q_dif <- quantile(data$dif_elev, c(0.025, 0.975)) ## extremes have outliers 
str(q_dif)
pred1 <- expand.grid(ABN = mean(data$ABN, na.rm = T),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = mean(data$LHF, na.rm = T),
                     LGI = mean(data$LGI, na.rm = T), 
                     elev = mean(data$elev, na.rm = T), 
                     dif_elev =  seq(q_dif[1], q_dif[2],length=50))

pred1$seedlings <- predict(model3, pred1)
dfr <- model3$df.residual
pred1$se <- predict(model3, pred1, se.fit = T)$se.fit * qt(0.95, dfr, lower.tail = F)

p1 <- ggplot(pred1, aes(x = dif_elev, y = exp(seedlings))) + 
  geom_ribbon(aes(ymax = exp(seedlings + se), ymin = exp(seedlings - se), fill = "dude"), alpha = 0.5) + 
  theme_bw() + 
  geom_line() + 
  geom_vline(xintercept = 0, linetype = 2, color = "#CC1C2B") +
  xlab("Elevation Differance (m)") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + 
  scale_fill_manual(values = "#004C75") + 
  theme(legend.position = "none")
  
p1  

ggsave(p1, file = "./graphs/elevation_differance_seedlings.png", 
       width = 4, height = 4)



# Predicting the models results quadratic 
pred2 <- expand.grid(ABN = mean(data$ABN, na.rm = T),
                     DBH = mean(data$DBH, na.rm = TRUE),
                     LHF = mean(data$LHF, na.rm = T),
                     LGI = mean(data$LGI, na.rm = T), 
                     dif_elev = mean(data$dif_elev, na.rm = T), 
                     elev =  seq(min(data$elev),max(data$elev),length=50))

pred2$seedlings <- predict(model3, pred2)
dfr <- model3$df.residual
pred2$se <- predict(model3, pred2, se.fit = T)$se.fit * qt(0.95, dfr)

p2 <- ggplot(pred2, aes(x = elev, y = exp(seedlings))) + 
  geom_ribbon(aes(ymax = exp(seedlings + se), ymin = exp(seedlings - se), fill = "#1E1E24"), 
              alpha = 0.5) + 
  theme_bw() + 
  geom_line() + 
  #geom_vline(xintercept = 0, linetype = 2, color = "red") +
  xlab("Elevation (m)") + 
  ylab(bquote('Seedlings 16'~ m^-2)) + 
  scale_fill_manual(values = "#004C75") + 
  theme(legend.position = "none")

p2

ggsave(p2, file = "./graphs/elevation_effect_seedlings.png", 
       width = 4, height = 4)
