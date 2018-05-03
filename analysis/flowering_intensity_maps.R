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

#' @title Flowering distribution map 
#' @discription Function that returns a ggplot of the flowering intensity distribution 
#' @param name takes a species name as a string 

flDist <- function(name){
  ggplot(data[data$sp == name & data$FL2010 == 0,], aes(x = X, y = Y)) + 
    geom_point() + 
    theme_bw() + 
    theme(legend.position = "null") + 
    geom_point(data = data[data$sp == name & data$FL2010 == 1,], alpha = 0.75,
               inherit.aes = F, aes(x = X, y = Y, size = log(GI + 1), fill = "p", color = "p", shape = "p")) + 
    scale_color_manual(values = c("black")) + 
    coord_fixed(ratio = 1) + 
    ylab("Latitude (m)") + 
    xlab("Longitude (m)") + 
    scale_shape_manual(values = c(21)) + 
    scale_fill_manual(values = "#6B0F1A") + 
    xlim(c(min(data$X), max(data$X))) + 
    ylim(c(min(data$Y), max(data$Y))) 
}

# vector of species names that have flowered
sp <- levels(droplevels(data[which(data$FL2010 == 1),])$sp)

#' A for loop to save all the different graphs 
for(i in 1:length(sp)){
  p1 <- flDist(sp[i])
  p1  
  ggsave(p1, file = paste("./graphs/flower_dist/", sp[i] ,"_fl_dist.png", sep = ""), 
         width = 8, height = 4)
}

# Using the same function as above but some added arguments
#' @param col colour of the plots 
#' @param al alpha level for the flowering points 
flDist <- function(name, col, al = 1){
  ggplot(data[data$sp == name & data$FL2010 == 0,], aes(x = X, y = Y)) + 
    geom_point() + 
    theme_bw() + 
    theme(legend.position = "null") + 
    geom_point(data = data[data$sp == name & data$FL2010 == 1,], alpha = al,
               inherit.aes = F, aes(x = X, y = Y, size = log(GI + 1), fill = "p", color = "p", shape = "p")) + 
    scale_color_manual(values = c("black")) + 
    coord_fixed(ratio = 1) + 
    ylab("") + 
    xlab("") + 
    scale_shape_manual(values = c(21)) + 
    scale_fill_manual(values = col) + 
    xlim(c(min(data$X), max(data$X))) + 
    ylim(c(min(data$Y), max(data$Y))) 
  
}

# Some colours 
cols = c("#86DEB7", "#6B0F1A", "#6C6EA0", "#F4D58D")
# so the preferance is for 
spp <- "Ssmi" # select species name 
p2 <- flDist(spp, col = cols[2], al = 0.75)
p2

# Save the graphic 
ggsave(p2, file = paste("./graphs/", spp ,"_fl_dist.png", sep = ""), 
       width = 6, height = 3)
