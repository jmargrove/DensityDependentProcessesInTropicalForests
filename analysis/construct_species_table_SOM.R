################################################################################
#' @title construct species table SOM
#' @author James Margrove

# Clear workspace 
rm(list=ls())

# Import data 
data <- read.table("./data/data.txt", header = T)
str(data)

# Abundnace as a frequency 
n.indv <- with(subset(data, DBH > 50), tapply(sp,sp,length))
n.dt <- with(data, tapply(sp,sp,length))
n.indv[which(is.na(n.indv))] <- 0
data$ABN <- rep(n.indv, n.dt)

# Creating a table of the important information 
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

write.table(dtable, file = "./data/table_of_species_information.txt")