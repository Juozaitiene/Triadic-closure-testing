library(igraph)
library(data.table)
library(ggplot2)
source("functions.R")

#------------Initialization of parameters------------
# number of iterations of simulation
ni<-50

# number of events
ne<-100

# number of sampling moments
ns<-50

#model parameters:
lambda_0 <- 2
mu_0 <- 2.5
lambda <- 0.5
mu <- -0.3

#number of actors:
n <- 50

variance <- rep(0,4)
answer <- rep(0,4)
#-------------------Simulation-------------------
for(z in 1:ni)
{

#vertices 
V <- 1:n
#links
data <- as.data.frame(matrix(0, n* (n-1)/2, 3))
data[,1:2] <- t(combn(V,2))
colnames(data) <- c("i", "j",0)
#number of links
m <- nrow(data)
#parameters
lambda_t <- rep(0,m)
mu_t <- rep(0,m)

#initial states
data[,3] <- rbinom(m,1,1/10)

#time intervals
#t <- rep(0,T)
t <- 0
i <- 3
repeat{
  G <- graph.data.frame(data[data[,i]!=0,c(1,2,i)], directed=F, vertices=V)
  adjSM <- as(get.adjacency(G), "dgTMatrix")
  adjDT <- data.table(V1=adjSM@i+1, V2=adjSM@j+1)
  res <- adjDT[adjDT, nomatch=0, on="V2", allow.cartesian=TRUE
               ][V1 < i.V1, .(Neighbours=paste(V2, collapse=",")),
                 by=c("V1","i.V1")][order(V1)]
  
  res[,3] <- apply(res[,3], 1, function(x) length(unlist(strsplit(as.character(x), ","))))
  neighbors <- as.data.frame(matrix(0, n* (n-1)/2, 2))
  neighbors[,1:2] <- t(combn(n,2))
  rez <- as.numeric(names(V(G)[unlist(res[,1])]))
  rez <- cbind(rez, as.numeric(names(V(G)[unlist(res[,2])])))
  rez <- cbind(rez,res[,3])
  colnames(neighbors) <- colnames(rez)[1:2]
  neighbors <- merge(neighbors, rez, by = c("rez", "V2"), all = TRUE)   
  neighbors[is.na(neighbors[,3]),3] <- 0
  neighbors <- neighbors[,3]
  #update parameters
  lambda_t <- lambda_0 + lambda * neighbors
  mu_t <- mu_0 + mu * neighbors
  #max(0, mu_0 + mu * neighbors)
  mu_t <- ifelse(mu_t > 0 , mu_t, 0)
  lambda_t <- ifelse(lambda_t > 0 ,lambda_t, 0)
  
  parameters <- ifelse(data[,i] == 0 , lambda_t , mu_t)
  
  #smallest time
  t <- c(t,rexp(1, sum(parameters)))
  
  #which event
  event <- rmultinom(1,1, prob = c(parameters)/sum(parameters))
  
  id <- which(event==1)
  data <- cbind(data,data[,i])
  data[id,i+1] <- ifelse(data[id,i]==1, 0, 1)

  
  #if(sum(t) > 1)
  if(ncol(data)==(ne+2))
  {
    break
  }
  i <- i+1
}
time <- rep(0, length(t))
for(i in 1:length(t))
{
  time[i] <- sum(t[1:i])
}

colnames(data)[3:ncol(data)] <- time

o <- sample(3:(ne+2), ns)
o <- o[order(o)]
time_t <- time[o-2]
dataset <- data[,c(1:2,o)]

answ <- estimation(dataset, time_t)
variance <- rbind(variance,answ[[2]])
answer <- rbind(answer,unlist(answ[1]))
}
colnames(answer) <- c("m_0","l_0","m","l")
answer <- answer[-1,]
variance <- variance[-1,]

#---------Parameters estimates for each simulation---------
head(answer)
#---------Variance estimates for each simulations---------
head(variance)

#---------Results visualization---------
dat <- data.frame(dens = c(answer[,1], answer[,2],
                           answer[,3], answer[,4]),
                  lines = rep(c("a", "b","c","d"), each = ni))
p <- ggplot(dat, aes(x = dens, fill = lines, color =lines))+ 
  geom_density(alpha = 0.5)+
  scale_x_continuous(breaks=c(-0.3, 0.5, 2, 2.5)) +
  labs(x="") 
p <- p + guides(color =FALSE)  
p <- p+ scale_fill_discrete(name = "", labels =c(expression(mu[0]),
                                                 expression(lambda[0]), 
                                                 expression(mu),
                                                 expression(lambda))) +
  theme_grey(base_size = 22) +  theme(legend.text.align = 0, 
                                      # Remove panel border
                                      panel.border = element_blank(),  
                                      # Remove panel grid lines
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      # Remove panel background
                                      panel.background = element_blank(),
                                      axis.line = element_line(colour = "black")) 
p <- p+ geom_vline(xintercept = mu)
p <- p+ geom_vline(xintercept = lambda)
p <- p+ geom_vline(xintercept = lambda_0)
p <- p+ geom_vline(xintercept = mu_0)

print(p)

