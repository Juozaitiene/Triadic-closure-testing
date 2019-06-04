#parameters estimation
estimation <- function(dataset, laikas)
{
  #number of time points
  T <- ncol(dataset)
  
  delta_t <- diff(laikas) 
  S1 <- sum(dataset[,3:(T-1)])
  S0 <- nrow(dataset)*(T-3) - S1
  
  S10 <- 0
  S01 <- 0
 
  
  #each time point
  for(i in 3:(T-1))
  {
    S10 <- S10 + sum(apply(dataset[,i:(i+1)], 1, function(x) ifelse(x[1]==1 && x[2]==0, 1, 0)))
    S01 <- S01 + sum(apply(dataset[,i:(i+1)], 1, function(x) ifelse(x[1]==0 && x[2]==1, 1, 0)))
  }
  #initial values
  l_0 <- -log(1-S01/S0)/mean(delta_t)
  m_0 <- -log(1-S10/S1)/mean(delta_t)
  l <- 0
  m <- 0
  
  previous <- rbind(m_0,l_0,m,l)
  ats <- functions(delta_t, dataset, c(m_0,l_0,m,l))
  answer <- previous - ats[[1]] %*% ats[[2]]
  
  repeat{
    previous <- answer
    ats <- functions(delta_t, dataset, c(previous))
    answer <- previous - ats[[1]] %*% ats[[2]]
    if(sum(abs(answer-previous)) < 0.001)
    {
      break
    }
    
  }
  
  d <- functions(delta_t, dataset, c(answer))
  
  var <- d[[1]] %*% diag(d[[3]]) %*% t(d[[1]])
  
  return(list(answer,diag(var), d[[3]]))
  
}

functions <- function(c_t, info, initial)
{
  #number of time points
  T <- ncol(info)
  
  #number of links
  links <- nrow(info)
  
  f_1_m_0 <- 0
  f_1_l_0 <- 0
  f_1_m <- 0
  f_1_l <- 0
  
  f_2_m_0 <- 0
  f_2_l_0 <- 0
  f_2_m <- 0
  f_2_l <- 0
  
  f_3_m_0 <- 0
  f_3_l_0 <- 0
  f_3_m <- 0
  f_3_l <- 0
  
  f_4_m_0 <- 0
  f_4_l_0 <- 0
  f_4_m <- 0
  f_4_l <- 0
  
  
  f_1 <- 0
  f_2 <- 0
  f_3 <- 0
  f_4 <- 0
  
  R0 <- 0
  R1 <- 0
  U0 <- 0
  U1 <- 0
  
  #List of graph nodes
  V <- unique(c(unlist(info[,1]), unlist(info[,2])))
  
  #for each time point
  for(i in 4:T)
  {
    #graph at time t-1
    g <- graph.data.frame(info[info[,i-1]!=0,c(1,2,i-1)], directed=F, vertices=V)
    #graph at time t
    g1 <- graph.data.frame(info[info[,i]!=0,c(1,2,i)], directed=F, vertices=V)
    
    #lines where first element is 0 at time t-1
    cl <- ifelse(info[,i-1]==0, 1,0) 
    a <- which(cl==1)
    #lines where first element is 1 at time t-1
    b <- which(cl==0)
    #lines where first element is 0 at time t
    cl <- ifelse(info[,i]==0, 1,0) 
    c <- which(cl==1)
    #lines where first element is 1 at time t
    d <- which(cl==0)
    
    #procedure to find number of common neighbours between each pair of nodes at time t-1
    adjSM <- as(get.adjacency(g), "dgTMatrix")
    adjDT <- data.table(V1=adjSM@i+1, V2=adjSM@j+1)
    res <- adjDT[adjDT, nomatch=0, on="V2", allow.cartesian=TRUE
                 ][V1 < i.V1, .(Neighbours=paste(V2, collapse=",")),
                   by=c("V1","i.V1")][order(V1)]
    
    res[,3] <- apply(res[,3], 1, function(x) length(unlist(strsplit(as.character(x), ","))))
    
    kaimynai <- as.data.frame(matrix(0, length(V)* (length(V)-1)/2, 2))
    kaimynai[,1:2] <- t(combn(V,2))
    
    rez <- as.numeric(names(V(g)[unlist(res[,1])]))
    rez <- cbind(rez, as.numeric(names(V(g)[unlist(res[,2])])))
    rez <- cbind(rez,res[,3])
    colnames(kaimynai) <- colnames(rez)[1:2]
    kaimynai <- merge(kaimynai, rez, by = c("rez", "V2"), all = TRUE)   
    kaimynai[is.na(kaimynai[,3]),3] <- 0
    
    #procedure to find number of common neighbours between each pair of nodes at time t
    adjSM <- as(get.adjacency(g1), "dgTMatrix")
    adjDT <- data.table(V1=adjSM@i+1, V2=adjSM@j+1)
    res1 <- adjDT[adjDT, nomatch=0, on="V2", allow.cartesian=TRUE
                  ][V1 < i.V1, .(Neighbours=paste(V2, collapse=",")),
                    by=c("V1","i.V1")][order(V1)]
    
    res1[,3] <- apply(res1[,3], 1, function(x) length(unlist(strsplit(as.character(x), ","))))
    
    kaimynai_1 <- as.data.frame(matrix(0, length(V)* (length(V)-1)/2, 2))
    kaimynai_1[,1:2] <- t(combn(V,2))
    
    rez <- as.numeric(names(V(g)[unlist(res1[,1])]))
    rez <- cbind(rez, as.numeric(names(V(g)[unlist(res1[,2])])))
    rez <- cbind(rez,res1[,3])
    colnames(kaimynai_1) <- colnames(rez)[1:2]
    kaimynai_1 <- merge(kaimynai_1, rez, by = c("rez", "V2"), all = TRUE)   
    kaimynai_1[is.na(kaimynai_1[,3]),3] <- 0
    
    #vector of numbers of common neighbours for each pair at time t-1
    N <- kaimynai[,3]
    #vector of numbers of common neighbours for each pair at time t
    N1 <- kaimynai_1[,3]
    
    #----------------------------------------------------
    #SAMPLE STATISTICS
    #----------------------------------------------------
    
    #link appeared
    index <- intersect(a,d)
    #not a part of open non-triangle
    U1 <- U1 + (length(index) - sum(N[index]!=0))
    #part of open triangle
    R1 <- R1 + sum(N[index]!=0)
    
    #link disappeared
    index <- intersect(b,c)
    #part of closed non-triangle
    U0 <- U0 + (length(index) - sum(N[index]!=0))
    #part of closed triangle
    R0 <- R0 + sum(N[index]!=0)
    
    #----------------------------------------------------
    #EXPECTED VALUES OF STATISTICS
    #----------------------------------------------------
    #nodes without neighbours
    q <- which(N==0)
    
    index <- intersect(a,q)
    #not a part of open triangle
    temp1 <- initial[2] 
    #m_0 + m * N(t_(k+1))
    temp2 <- initial[1] + initial[3] * N1[index]
    #indicator whether mu_0 + N_ij * mu is equal to 0
    temp3 <- ifelse(temp2 > 0 , 1, 0)
    temp2 <- ifelse(temp2 > 0 , temp2, 0)
    
    laik1 <- (-(temp1)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] )  * temp3
    f_2_m_0 <- f_2_m_0 + sum(laik1)
    
    laik2 <- (-(temp1)/(temp1 + temp2)^2 * N1[index] * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N1[index]) * temp3
    f_2_m <- f_2_m + sum(laik2)
    
    laik3 <- (temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] 
    f_2_l_0 <- f_2_l_0 + sum(laik3)
    
    f_2_l <- 0
    
    laik5 <- temp1/(temp1 + temp2) * (1-exp(-(temp1 + temp2)*c_t[i-3]))
    f_2 <- f_2 + sum(laik5)
    
    #nodes having neigbours
    q <- which(N!=0)
    index <- intersect(a,q)
    #part of open triangle
    #l_0 + l * N(t_k)
    temp1 <- initial[2] + initial[4] * N[index]
    #m_0 + m * N(t_(k+1))
    temp2 <- initial[1] + initial[3] * N1[index]
    #indicator whether mu_0 + N_ij * mu is equal to 0
    temp3 <- ifelse(temp2 > 0 , 1, 0)
    #indicator whether lambda_0 + N_ij * lambda is equal to 0
    temp4 <- ifelse(temp1 > 0 , 1, 0)
    
    temp1 <- ifelse(temp1 > 0 , temp1, 0)
    temp2 <- ifelse(temp2 > 0 , temp2, 0)
    
    laik1 <- (-(temp1)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3]) * temp3
    f_1_m_0 <- f_1_m_0 + sum(laik1)
    
    laik2 <- (-(temp1)/(temp1 + temp2)^2 * N1[index] * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N1[index]) * temp3
    f_1_m <- f_1_m + sum(laik2)
    
    laik3 <- ((temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] ) * temp4
    f_1_l_0 <- f_1_l_0 + sum(laik3)
    
    laik4 <- ((N[index] * temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp1/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N[index]) *temp4
    f_1_l <- f_1_l + sum(laik4)
    
    laik5 <- temp1/(temp1 + temp2) * (1-exp(-(temp1 + temp2)*c_t[i-3]))
    f_1 <- f_1 + sum(laik5)
    
    q <- which(N==0)
    index <- intersect(b,q)
    #not a part of closed triangle
    #l_0 + l * N(t_(k+1))
    temp1 <- initial[2] + initial[4] * N1[index]
    #m_0 
    temp2 <- initial[1] 
    #indicator whether lambda_0 + N_ij * lambda is equal to 0
    temp4 <- ifelse(temp1 > 0 , 1, 0)
    
    temp1 <- ifelse(temp1 > 0 , temp1, 0)
    
    laik1 <- (temp1)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3]
    f_4_m_0 <- f_4_m_0 + sum(laik1)
    
    f_4_m <- 0
    
    laik3 <- (-(temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] )*temp4
    f_4_l_0 <- f_4_l_0 + sum(laik3)
    
    laik4 <- (-(N1[index] * temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N1[index])*temp4
    f_4_l <- f_4_l + sum(laik4)
    
    laik5 <- temp2/(temp1 + temp2) * (1-exp(-(temp1 + temp2)*c_t[i-3]))
    f_4 <- f_4 + sum(laik5)
    
    q <- which(N!=0)
    index <- intersect(b,q)
    #part of closed triangle
    #l_0 + l * N(t_(k+1))
    temp1 <- initial[2] + initial[4] * N1[index]
    #m_0 + m * N(t_k)
    temp2 <- initial[1] + initial[3] * N[index]
    #indicator whether mu_0 + N_ij * mu is equal to 0
    temp3 <- ifelse(temp2 > 0 , 1, 0)
    #indicator whether lambda_0 + N_ij * lambda is equal to 0
    temp4 <- ifelse(temp1 > 0 , 1, 0)
    
    temp1 <- ifelse(temp1 > 0 , temp1, 0)
    temp2 <- ifelse(temp2 > 0 , temp2, 0)
    
    laik1 <- ((temp1)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3]) * temp3
    f_3_m_0 <- f_3_m_0 + sum(laik1)
    
    laik2 <- ((temp1 * N[index])/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N[index]) * temp3
    f_3_m <- f_3_m + sum(laik2)
    
    laik3 <- (-(temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] ) * temp4
    f_3_l_0 <- f_3_l_0 + sum(laik3)
    
    laik4 <- (-(N1[index] * temp2)/(temp1 + temp2)^2 * (1-exp(-(temp1 + temp2)*c_t[i-3])) + temp2/(temp1 + temp2) * exp(-(temp1 + temp2)*c_t[i-3]) * c_t[i-3] * N1[index]) * temp4
    f_3_l <- f_3_l + sum(laik4)
    
    laik5 <- temp2/(temp1 + temp2) * (1-exp(-(temp1 + temp2)*c_t[i-3]))
    f_3 <- f_3 + sum(laik5)
  }
  
  l <- matrix(c(f_1_m_0, f_1_l_0, f_1_m, f_1_l,
                f_2_m_0, f_2_l_0, f_2_m, f_2_l,
                f_3_m_0, f_3_l_0, f_3_m, f_3_l,
                f_4_m_0, f_4_l_0, f_4_m, f_4_l), nrow = 4,ncol = 4, byrow = T)
  l <- solve(l)
  
  f <- rbind(f_1 - R1,f_2 - U1,f_3 - R0,f_4 - U0)
  S <- c(R1, U1, R0, U0)
  
  return(list(l,f, S))
}

Hotelling_test <- function(theta_hat, sigma,n)
{
  T2 <- t(theta_hat) %*% solve(sigma) %*% theta_hat
  T2 <- (n-2)/(2*(n-1)) * T2
  p <- 1 - pf(T2,2,n-2,lower.tail = T)
  return(p)
}