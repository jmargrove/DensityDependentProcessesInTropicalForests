################################################################################
#' @author James margrove 
#' @title hypothesis graphics 
#' 
#' 

rm(list=ls())



x1 <- runif(100)
require(ggplot2)


quad <- function(b1, b2, alpha, x){
  x * b1 + x^2 * b2 + alpha 
}




test <- data.frame(x1 - 0.5, y = quad(12, -1 * 12, 0, x1))


ggplot(test, aes(x = x1, y = exp(y))) + geom_line() + theme_bw()


# x1 > 
y2 <- test$y
y2[which(test$x1 > 0)] <- max(test$y)
test$y2 <- y2

ggplot(test, aes(x = x1, y = y)) + geom_line() + 
  geom_line(aes(x = x1, y = y2 + 0.05), color = "red", linetype = 2)


# x1 > 
y3 <- test$y
y3[which(test$x1 < 0)] <- max(test$y)
test$y3 <- y3   

ggplot(test, aes(x = x1, y = y)) + geom_line() + 
  theme_bw() + 
  geom_line(aes(x = x1, y = y2 + 0.01), color = "red", linetype = 2) + 
  geom_line(aes(x = x1, y = y3 + 0.02), color = "blue", linetype = 2)


# Invert
y4 <- test$y
y4 <- y4 * -1 + max(test$y)
test$y4 <- y4   


p1 <- ggplot(test, aes(x = x1, y = exp(y))) + geom_line( color = "black", linetype = 2) + 
  theme_bw() + 
  geom_line(aes(x = x1, y = exp(y2 + 0.01)), color = "red", linetype = 2) + 
  geom_line(aes(x = x1, y = exp(y3 + 0.02)), color = "blue", linetype = 2) + 
  geom_line(aes(x = x1, y = exp(y4 + 0.0)), color = "green", linetype = 2) + 
  ylab("Recruitment success") + 
  xlab("Conspecific flowering intensity")


p1

### now abundance 
linear <- function(beta, alpha, x){
  x * beta + alpha 
}
max(linear(4, 0, x1))
n = max(test$y)
abn_df <- data.frame(x1 = x1 - 0.5, y = linear(n, 0, x1), y1 = linear(-1*n, max(linear(n, 0, x1)), x1))

p1 <- ggplot(abn_df, aes(x = x1, y = exp(y))) + geom_line() + theme_bw() + 
  geom_line(aes(x1, y = exp(y1))) + 
  geom_line(data = test, inherit.aes = F, aes(x = x1 - 0.5, y = exp(y - 0.01)), 
            color = "green") + 
  geom_line(data = test, inherit.aes = F, aes(x = x1 - 0.5, y = exp(y2 + 0.03)), 
            color = "red", linetype = 2) + 
  geom_line(data = test, inherit.aes = F, aes(x = x1 - 0.5, y = exp(y3 + 0.03)), 
            color = "blue", linetype = 2) +
  ylab("Seedlings") +
  xlab("Conspecific flowering intensity")

p1

# This is the tree 
ggsave(p1, file = "./graphs/hyp_LGI.png", 
       width = 4, height = 4)

#### additional with abundance 

quad <- function(b2, b1, alpha, x){
  (x^2 * b2) + (x * b1 ) + alpha 
}


x1 <- seq(-10, 10, length = 100)



res = c();
b <- c(0.02, 0.05, 0.1, 0.75, 2) * -1
test2 <- data.frame(x1 = rep(x1, times = 5), qt = rep(b, each = 100)) 
for(i in b){
  res <- append(quad(i,0, log(10), x1), res)
}

test2$y2 <- res
test2$qt <- as.factor(test2$qt)

colnames(test2)[2] <- "Quadratic term"

p2 <- ggplot(test2, aes(x = x1, y = exp(y2), linetype = `Quadratic term`)) + geom_line() + 
  theme_classic() + 
  xlim(-10, 0) + 
  theme(legend.position = "none", axis.text = element_blank()) + 
  ylab("") + 
  xlab("") 

p2

ggsave(p2, file = "./graphs/quad term hypothesis.png", 
       width = 4, height = 4)


p3 <- ggplot(test2, aes(x = x1, y = exp(y2), linetype = `Quadratic term`)) + geom_line() + 
  theme_classic() + 
  xlim(0, 10) + 
  theme(legend.position = "none", axis.text = element_blank()) + 
  ylab("") + 
  xlab("") 

p3

ggsave(p3, file = "./graphs/quad term hypothesis_rev.png", 
       width = 4, height = 4)





#######

quad <- function(b2, b1, alpha, x){
  (x^2 * b2) + (x * b1 ) + alpha 
}


x1 <- seq(-10, 10, length = 100)



res = c();
b <- c(12.5, 10, 7.5, 5, 2.5)
test2 <- data.frame(x1 = rep(x1, times = 5), qt = rep(b, each = 100)) 
for(i in b){
  res <- append(quad(-0.2,0, log(i), x1), res)
}

test2$y2 <- res
test2$qt <- as.factor(test2$qt)

colnames(test2)[2] <- "Heterospecific flowering intensity"

p4 <- ggplot(test2, aes(x = x1, y = exp(y2), linetype = `Heterospecific flowering intensity`)) + geom_line() + 
  theme_classic() + 
  theme(legend.position = "null", axis.text = element_blank()) + 
  ylab("") + 
  xlim(-5,5) + 
  xlab("") 

p4

ggsave(p4, file = "./graphs/quad term hypothesis.png", 
       width = 4, height = 4)




p5 <- ggplot(test2, aes(x = (x1), y = exp(y2)*-1, linetype = `Heterospecific flowering intensity`)) + geom_line() + 
  theme_classic() + 
  theme(legend.position = "top", axis.text = element_blank()) + 
  ylab("") + 
  xlim(-5,5) + 
  xlab("") 

p5

ggsave(p5, file = "./graphs/quad term hypothesis_rev.png", 
       width = 4, height = 4)




test3 <- test2
test3$`Heterospecific flowering intensity` <-  rev(test3$`Heterospecific flowering intensity`)
test4 <- rbind(test2, test3)
test4$abn <-  rep(c("high", "low"), each = dim(test2)[1])

p3 <- ggplot(test4, aes(x = x1, y = exp(y2), color = `Heterospecific flowering intensity`)) + 
  geom_line() + 
  facet_grid(~abn) + 
  theme_classic() + 
  theme(legend.position = "top", axis.text = element_blank()) + 
  ylab("Reproductive output") + 
  xlim(-5,0) + 
  xlab("Heterospecific flowering intensity") 
  

p3

ggsave(p3, file = "./graphs/quad term hypothesis.png", 
       width = 4, height = 4)




### Davids hypothesis

#1 
x <- 1:100








