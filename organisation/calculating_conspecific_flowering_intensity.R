##########################################################################################
#' @title Calculating the conspecific flowering intensity 
#' @author James Margrove 

# Clear workspace 
rm(list=ls())

# Import data 
spdata <- read.table("./data/spdata.txt", header = T)
spdata <- subset(spdata, dead != "DEAD2010")
spdata$seedlings <- apply(spdata[,26:29], 1, sum, na.rm = T)
spdata[which(spdata$DBH < 20),"DBH"] <- NA
str(spdata)


# Distance function 
distance <- function(x1,x2, y1, y2){
  sqrt((x1-x2)^2+(y1-y2)^2)}

finaldata <- spdata[0,]

##### the function to calculate the force of flowering :) using gravetry function 
for(species in levels(spdata$sp)){
  # Select a dataframe with only conspecific individuals 
  dt_species <- subset(spdata, sp == species)
  # Assign species with NA DBH the species mean 
  dt_species$DBH[which(is.na(dt_species$DBH))] <- mean(dt_species$DBH, na.rm = T)
  
  if(dim(dt_species)[1] == 1){
    Sum_FF <- 0 ### incase there is only one speceis the force value will equal 0!
    }else{ 
    mat_size <- dim(dt_species)[1]
    dist_matrix <- matrix(rep(NA, mat_size*mat_size), nrow = mat_size)
    
    coords <- as.matrix(dt_species[,c(23,24)])
    
    for(j in 1:mat_size){
      sp1 <- coords[j,c(1,2)]
      for(i in 1:mat_size){
        spN <- coords[i,c(1,2)]
        dist_matrix[i,j] <- distance(sp1[1], spN[1], sp1[2], spN[2])
      }}
    
    rownames(dist_matrix) <- dt_species[1:mat_size,1]
    colnames(dist_matrix) <- dt_species[1:mat_size,1]
    
    dim(dist_matrix)
    #dist_matrix[1:5,1:5]
    F_matrix <- matrix(rep(NA, mat_size*mat_size), nrow = mat_size)
    tree_size <- dt_species$DBH # the tree sizes (dbh)
    tree_size[1:10]
    for(j in 1:mat_size){
      indv_dist_vec <- dist_matrix[,j] # distance col 
      for(i in 1:mat_size){
        F_matrix[i,j] <-    1/indv_dist_vec[i]^2}} # Force matrix ... 
    
    rownames(F_matrix) <- dt_species[1:mat_size,1]
    colnames(F_matrix) <- dt_species[1:mat_size,1]
    #F_matrix[1:10,1:10]
    F_matrix[which(is.infinite(F_matrix))] <- 0 # this is changing the inf values to zero (1/0 = inf) F_matrix[which(is.infinite(F_matrix))]
    ######################################################################
    flowering_trees <- which(dt_species$FL2010 == 1)# which are the flowering trees
    F_matrix_fl <- F_matrix[flowering_trees,] # matrix with rows of flowering trees, col of all trees
    
    #F_matrix_fl[1:5,1:5]
    #######################################################################
    #the sum of all the Forces is equal to the closeness of the nearest flowering neighbor 
    
    # if statement as cannot use teh apply function to sum over a vector. there is no need
    if(class(F_matrix_fl) == "matrix") {
      Sum_FF <- apply(F_matrix_fl, 2, sum, na.rm = T)}else{
        Sum_FF <- F_matrix_fl
      } # this is the sum of the flowering force on all trees
  }
  dt_species$Sum_FF <- Sum_FF # add to the final dataframe
  
  finaldata <- rbind(finaldata, dt_species)}

dim(finaldata)
head(finaldata)


result_data$Sum_FF3 <- finaldata$Sum_FF
head(result_data)

write.table(result_data, file = "./data/calc_data.txt")
