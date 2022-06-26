library("matrixcalc")
library("markovchain")
library("Rcpp")
library("Matrix")


#Power without prior info
power1<-matrix(data = c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0),byrow = FALSE, nrow = 5, dimnames = list(c("100","500","1000","5000","10000"), c("0","1/10","1/50","1/100")))
for (c in c(0,1/10,1/50,1/100)) {
  powerv <- c()
  for (n in c(100,500,1000,5000,10000)) {
    
    sample_size <- n
    
    c
    
    States <- c("1", "2", "3")
    
    byRow <- TRUE
    
    TransMatrix <- matrix(data = c(1/3, 1/3, 1/3,
                                   2/6, 1/6, 3/6,
                                   (1+6*c)/6, (2+6*c)/6, (3-12*c)/6), byrow = byRow, nrow = 3, 
                          dimnames = list(States, States))
    
    w <- matrix(data = c(1, 1, 1,
                         1, 1, 1,
                         1, 1, 1), nrow = 3, ncol = 3,byrow = T)
    
    
    stationary <- matrix.power(TransMatrix,1000) #stationary distribution
    
    sum_wp <- function(w, TransMatrix,stationary){
      sum <- 0
      for (i in 1:length(States)) {
        for (j in 1:length(States)) {
          sum <- sum + w[i,j]*stationary[1,i]*TransMatrix[i,j]
        }
      }
      return(sum)
    }
    
    C <- matrix(1:81, nrow =9, ncol = 9)
    
    C[1,]<-c(((w[1,1]*stationary[1,1])/(TransMatrix[1,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0,0)
    C[2,]<-c(0,((w[1,2]*stationary[1,1])/(TransMatrix[1,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0)
    C[3,]<-c(0,0,((w[1,3]*stationary[1,1])/(TransMatrix[1,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0)
    C[4,]<-c(0,0,0,((w[2,1]*stationary[1,2])/(TransMatrix[2,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0)
    C[5,]<-c(0,0,0,0,((w[2,2]*stationary[1,2])/(TransMatrix[2,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0)
    C[6,]<-c(0,0,0,0,0,((w[2,3]*stationary[1,2])/(TransMatrix[2,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0)
    C[7,]<-c(0,0,0,0,0,0,((w[3,1]*stationary[1,3])/(TransMatrix[3,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0)
    C[8,]<-c(0,0,0,0,0,0,0,((w[3,2]*stationary[1,3])/(TransMatrix[3,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0)
    C[9,]<-c(0,0,0,0,0,0,0,0,((w[3,3]*stationary[1,3])/(TransMatrix[3,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))))
    
    Sigma <- matrix(1:81, nrow =9, ncol = 9)
    
    Sigma[1,]<-c((1/stationary[1,1])*(TransMatrix[1,1]*(1-TransMatrix[1,1])),(-TransMatrix[1,1]*TransMatrix[1,2])/stationary[1,1],(-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[2,]<-c((-TransMatrix[1,2]*TransMatrix[1,1])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,2]*(1-TransMatrix[1,2])),(-TransMatrix[1,2]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[3,]<-c((-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],(-TransMatrix[1,3]*TransMatrix[1,2])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,3]*(1-TransMatrix[1,3])),0,0,0,0,0,0)
    Sigma[4,]<-c(0,0,0,(1/stationary[1,2])*(TransMatrix[2,1]*(1-TransMatrix[2,1])),(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],0,0,0)
    Sigma[5,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,2]*(1-TransMatrix[2,2])),(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],0,0,0)
    Sigma[6,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,3]*(1-TransMatrix[2,3])),0,0,0)
    Sigma[7,]<-c(0,0,0,0,0,0,(1/stationary[1,3])*(TransMatrix[3,1]*(1-TransMatrix[3,1])),(-TransMatrix[3,1]*TransMatrix[3,2])/stationary[1,3],(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3])
    Sigma[8,]<-c(0,0,0,0,0,0,(-TransMatrix[3,2]*TransMatrix[3,1])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,2]*(1-TransMatrix[3,2])),(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3])
    Sigma[9,]<-c(0,0,0,0,0,0,(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3],(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,3]*(1-TransMatrix[3,3])))
    
    CSigma <- C%*%Sigma
    
    ev <- eigen(CSigma)
    
    betas <- ev$values
    betas
    
    b1 <- betas[1]
    b2 <- betas[2]
    b3 <- betas[3]
    b4 <- betas[4]
    b5 <- betas[5]
    b6 <- betas[6]
    b7 <- betas[7]
    b8 <- betas[8]
    b9 <- betas[9]
    
    distr <- c()
    for(i in 1:10000){
      
      r.chi <- b1*rchisq(n=10000, df=1) + b2*rchisq(n=10000, df=1) + b3*rchisq(n=10000, df=1)+ b4*rchisq(n=10000, df=1)+ b5*rchisq(n=10000, df=1) + b6*rchisq(n=10000, df=1) + b7*rchisq(n=10000, df=1) + b8*rchisq(n=10000, df=1)+ b9*rchisq(n=10000, df=1)
      
      distr[i]<- r.chi
    }
    
    q.distr <- quantile(distr, probs =seq(0, 1, 0.01))
    
    TransMatrix.alt <- matrix(data = c(1/3, 1/3, 1/3,
                                       2/6, 1/6, 3/6,
                                       1/6, 2/6, 3/6), byrow = byRow, nrow = 3, 
                              dimnames = list(States, States))
    
    mc <- new("markovchain", states = States, byrow = byRow,
              transitionMatrix = TransMatrix.alt, name = "Markov Model")
    
    CWKL1 <-c()
    for (i in 1:1000) {
      
      sample <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est <- function(X, prob=T){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est(X=sample,prob = T)
      
      CWKL.1 <- function(TransMatrix, w, TransMatrix.est, stationary){
        sum <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum <- sum + (stationary[1,i]*w[i,j]*(TransMatrix.est[i,j]*log((TransMatrix.est[i,j])/(TransMatrix[i,j]))-TransMatrix.est[i,j]+TransMatrix[i,j]))
          }
        }
        return(sum)
      } 
      
      CWKL.2 <- function(TransMatrix, w, stationary){
        sum1 <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum1 <- sum1 + (w[i,j]*stationary[1,i]*TransMatrix.est[i,j])
          }
        }
        return(sum1)
      } 
      
      CWKL.est.1 <- CWKL.1(TransMatrix=TransMatrix, w=w, TransMatrix.est=TransMatrix.est, stationary=stationary)
      
      CWKL.est.2 <- CWKL.2(TransMatrix=TransMatrix, w=w, stationary=stationary)
      
      CWKL.est <- (CWKL.est.1)/(CWKL.est.2) #CWKL divergence
      
      CWKL <- (2*sample_size*CWKL.est)
      
      CWKL1[i] <- CWKL
      
      CWKL2 <- CWKL1[!is.na(CWKL1)]
      
      u <- CWKL2 >= q.distr[96]
      power <- length(u[u==TRUE])/length(u)
    }
    powerv[n] <- power
    powerv<-powerv[!is.na(powerv)]
  }
  if(c==0){
    power1[,1]<-powerv
  }
  
  if(c==1/10){
    power1[,2]<-powerv
  }
  
  if(c==1/50){
    power1[,3]<-powerv
  }
  
  if(c==1/100){
    power1[,4]<-powerv
  }
}  




#Power with Prior info 1
power2<-matrix(data = c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0),byrow = FALSE, nrow = 5, dimnames = list(c("100","500","1000","5000","10000"), c("0","1/10","1/50","1/100")))
for (c in c(0,1/10,1/50,1/100)) {
  powerv <- c()
  for (n in c(100,500,1000,5000,10000)) {
    
    sample_size <- n
    
    c
    
    States <- c("1", "2", "3")
    
    byRow <- TRUE
    
    TransMatrix <- matrix(data = c(1/3, 1/3, 1/3,
                                   2/6, 1/6, 3/6,
                                   (1+6*c)/6, (2+6*c)/6, (3-12*c)/6), byrow = byRow, nrow = 3, 
                          dimnames = list(States, States))
    
    w <- matrix(data = c(1, 1, 1,
                         1, 1, 1,
                         2, 2, 2), nrow = 3, ncol = 3,byrow = T)
    
    
    stationary <- matrix.power(TransMatrix,1000) #stationary distribution
    
    sum_wp <- function(w, TransMatrix,stationary){
      sum <- 0
      for (i in 1:length(States)) {
        for (j in 1:length(States)) {
          sum <- sum + w[i,j]*stationary[1,i]*TransMatrix[i,j]
        }
      }
      return(sum)
    }
    
    C <- matrix(1:81, nrow =9, ncol = 9)
    
    C[1,]<-c(((w[1,1]*stationary[1,1])/(TransMatrix[1,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0,0)
    C[2,]<-c(0,((w[1,2]*stationary[1,1])/(TransMatrix[1,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0)
    C[3,]<-c(0,0,((w[1,3]*stationary[1,1])/(TransMatrix[1,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0)
    C[4,]<-c(0,0,0,((w[2,1]*stationary[1,2])/(TransMatrix[2,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0)
    C[5,]<-c(0,0,0,0,((w[2,2]*stationary[1,2])/(TransMatrix[2,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0)
    C[6,]<-c(0,0,0,0,0,((w[2,3]*stationary[1,2])/(TransMatrix[2,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0)
    C[7,]<-c(0,0,0,0,0,0,((w[3,1]*stationary[1,3])/(TransMatrix[3,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0)
    C[8,]<-c(0,0,0,0,0,0,0,((w[3,2]*stationary[1,3])/(TransMatrix[3,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0)
    C[9,]<-c(0,0,0,0,0,0,0,0,((w[3,3]*stationary[1,3])/(TransMatrix[3,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))))
    
    Sigma <- matrix(1:81, nrow =9, ncol = 9)
    
    Sigma[1,]<-c((1/stationary[1,1])*(TransMatrix[1,1]*(1-TransMatrix[1,1])),(-TransMatrix[1,1]*TransMatrix[1,2])/stationary[1,1],(-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[2,]<-c((-TransMatrix[1,2]*TransMatrix[1,1])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,2]*(1-TransMatrix[1,2])),(-TransMatrix[1,2]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[3,]<-c((-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],(-TransMatrix[1,3]*TransMatrix[1,2])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,3]*(1-TransMatrix[1,3])),0,0,0,0,0,0)
    Sigma[4,]<-c(0,0,0,(1/stationary[1,2])*(TransMatrix[2,1]*(1-TransMatrix[2,1])),(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],0,0,0)
    Sigma[5,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,2]*(1-TransMatrix[2,2])),(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],0,0,0)
    Sigma[6,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,3]*(1-TransMatrix[2,3])),0,0,0)
    Sigma[7,]<-c(0,0,0,0,0,0,(1/stationary[1,3])*(TransMatrix[3,1]*(1-TransMatrix[3,1])),(-TransMatrix[3,1]*TransMatrix[3,2])/stationary[1,3],(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3])
    Sigma[8,]<-c(0,0,0,0,0,0,(-TransMatrix[3,2]*TransMatrix[3,1])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,2]*(1-TransMatrix[3,2])),(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3])
    Sigma[9,]<-c(0,0,0,0,0,0,(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3],(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,3]*(1-TransMatrix[3,3])))
    
    CSigma <- C%*%Sigma
    
    ev <- eigen(CSigma)
    
    betas <- ev$values
    betas
    
    b1 <- betas[1]
    b2 <- betas[2]
    b3 <- betas[3]
    b4 <- betas[4]
    b5 <- betas[5]
    b6 <- betas[6]
    b7 <- betas[7]
    b8 <- betas[8]
    b9 <- betas[9]
    
    distr <- c()
    for(i in 1:10000){
      
      r.chi <- b1*rchisq(n=10000, df=1) + b2*rchisq(n=10000, df=1) + b3*rchisq(n=10000, df=1)+ b4*rchisq(n=10000, df=1)+ b5*rchisq(n=10000, df=1) + b6*rchisq(n=10000, df=1) + b7*rchisq(n=10000, df=1) + b8*rchisq(n=10000, df=1)+ b9*rchisq(n=10000, df=1)
      
      distr[i]<- r.chi
    }
    
    q.distr <- quantile(distr, probs =seq(0, 1, 0.01))
    
    TransMatrix.alt <- matrix(data = c(1/3, 1/3, 1/3,
                                       2/6, 1/6, 3/6,
                                       1/6, 2/6, 3/6), byrow = byRow, nrow = 3, 
                              dimnames = list(States, States))
    
    mc <- new("markovchain", states = States, byrow = byRow,
              transitionMatrix = TransMatrix.alt, name = "Markov Model")
    
    CWKL1 <-c()
    for (i in 1:1000) {
      
      sample <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est <- function(X, prob=T){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est(X=sample,prob = T)
      
      CWKL.1 <- function(TransMatrix, w, TransMatrix.est, stationary){
        sum <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum <- sum + (stationary[1,i]*w[i,j]*(TransMatrix.est[i,j]*log((TransMatrix.est[i,j])/(TransMatrix[i,j]))-TransMatrix.est[i,j]+TransMatrix[i,j]))
          }
        }
        return(sum)
      } 
      
      CWKL.2 <- function(TransMatrix, w, stationary){
        sum1 <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum1 <- sum1 + (w[i,j]*stationary[1,i]*TransMatrix.est[i,j])
          }
        }
        return(sum1)
      } 
      
      CWKL.est.1 <- CWKL.1(TransMatrix=TransMatrix, w=w, TransMatrix.est=TransMatrix.est, stationary=stationary)
      
      CWKL.est.2 <- CWKL.2(TransMatrix=TransMatrix, w=w, stationary=stationary)
      
      CWKL.est <- (CWKL.est.1)/(CWKL.est.2) #CWKL divergence
      
      CWKL <- (2*sample_size*CWKL.est)
      
      CWKL1[i] <- CWKL
      
      CWKL2 <- CWKL1[!is.na(CWKL1)]
      
      u <- CWKL2 >= q.distr[96]
      power <- length(u[u==TRUE])/length(u)
    }
    powerv[n] <- power
    powerv<-powerv[!is.na(powerv)]
  }
  if(c==0){
    power2[,1]<-powerv
  }
  
  if(c==1/10){
    power2[,2]<-powerv
  }
  
  if(c==1/50){
    power2[,3]<-powerv
  }
  
  if(c==1/100){
    power2[,4]<-powerv
  }
}  




#Power with prior info 2
power3<-matrix(data = c(0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0),byrow = FALSE, nrow = 5, dimnames = list(c("100","500","1000","5000","10000"), c("0","1/10","1/50","1/100")))
for (c in c(0,1/10,1/50,1/100)) {
  powerv <- c()
  for (n in c(100,500,1000,5000,10000)) {
    
    sample_size <- n
    
    c
    
    States <- c("1", "2", "3")
    
    byRow <- TRUE
    
    TransMatrix <- matrix(data = c(1/3, 1/3, 1/3,
                                   2/6, 1/6, 3/6,
                                   (1+6*c)/6, (2+6*c)/6, (3-12*c)/6), byrow = byRow, nrow = 3, 
                          dimnames = list(States, States))
    
    TransMatrix.alt <- matrix(data = c(1/3, 1/3, 1/3,
                                       2/6, 1/6, 3/6,
                                       1/6, 2/6, 3/6), byrow = byRow, nrow = 3, 
                              dimnames = list(States, States))
    
    mc <- new("markovchain", states = States, byrow = byRow,
              transitionMatrix = TransMatrix.alt, name = "Markov Model")
    
    
    #w <-abs(TransMatrix.alt-TransMatrix)
    
    w <- array(NA,dim=c(3,3,1000))
    for (i in 1:1000) {
      
      sample1 <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est1 <- function(X, prob=T){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est1(X=sample1,prob = T)
      
      w[,,i] <-abs(TransMatrix.est-TransMatrix)
    }
    
    sum11 <- matrix(0,3,3)
    for (i in 1:1000) {
      sum11 <- sum11 + w[,,i]
    }
    w <- sum11*(1/1000)
    
    
    stationary <- matrix.power(TransMatrix,1000) #stationary distribution
    
    sum_wp <- function(w, TransMatrix,stationary){
      sum <- 0
      for (i in 1:length(States)) {
        for (j in 1:length(States)) {
          sum <- sum + w[i,j]*stationary[1,i]*TransMatrix[i,j]
        }
      }
      return(sum)
    }
    
    C <- matrix(1:81, nrow =9, ncol = 9)
    
    C[1,]<-c(((w[1,1]*stationary[1,1])/(TransMatrix[1,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0,0)
    C[2,]<-c(0,((w[1,2]*stationary[1,1])/(TransMatrix[1,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0,0)
    C[3,]<-c(0,0,((w[1,3]*stationary[1,1])/(TransMatrix[1,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0,0)
    C[4,]<-c(0,0,0,((w[2,1]*stationary[1,2])/(TransMatrix[2,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0,0)
    C[5,]<-c(0,0,0,0,((w[2,2]*stationary[1,2])/(TransMatrix[2,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0,0)
    C[6,]<-c(0,0,0,0,0,((w[2,3]*stationary[1,2])/(TransMatrix[2,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0,0)
    C[7,]<-c(0,0,0,0,0,0,((w[3,1]*stationary[1,3])/(TransMatrix[3,1]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0,0)
    C[8,]<-c(0,0,0,0,0,0,0,((w[3,2]*stationary[1,3])/(TransMatrix[3,2]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))),0)
    C[9,]<-c(0,0,0,0,0,0,0,0,((w[3,3]*stationary[1,3])/(TransMatrix[3,3]*sum_wp(w=w, TransMatrix = TransMatrix, stationary = stationary))))
    
    Sigma <- matrix(1:81, nrow =9, ncol = 9)
    
    Sigma[1,]<-c((1/stationary[1,1])*(TransMatrix[1,1]*(1-TransMatrix[1,1])),(-TransMatrix[1,1]*TransMatrix[1,2])/stationary[1,1],(-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[2,]<-c((-TransMatrix[1,2]*TransMatrix[1,1])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,2]*(1-TransMatrix[1,2])),(-TransMatrix[1,2]*TransMatrix[1,3])/stationary[1,1],0,0,0,0,0,0)
    Sigma[3,]<-c((-TransMatrix[1,1]*TransMatrix[1,3])/stationary[1,1],(-TransMatrix[1,3]*TransMatrix[1,2])/stationary[1,1],(1/stationary[1,1])*(TransMatrix[1,3]*(1-TransMatrix[1,3])),0,0,0,0,0,0)
    Sigma[4,]<-c(0,0,0,(1/stationary[1,2])*(TransMatrix[2,1]*(1-TransMatrix[2,1])),(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],0,0,0)
    Sigma[5,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,2]*(1-TransMatrix[2,2])),(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],0,0,0)
    Sigma[6,]<-c(0,0,0,(-TransMatrix[2,1]*TransMatrix[2,3])/stationary[1,2],(-TransMatrix[2,3]*TransMatrix[2,2])/stationary[1,2],(1/stationary[1,2])*(TransMatrix[2,3]*(1-TransMatrix[2,3])),0,0,0)
    Sigma[7,]<-c(0,0,0,0,0,0,(1/stationary[1,3])*(TransMatrix[3,1]*(1-TransMatrix[3,1])),(-TransMatrix[3,1]*TransMatrix[3,2])/stationary[1,3],(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3])
    Sigma[8,]<-c(0,0,0,0,0,0,(-TransMatrix[3,2]*TransMatrix[3,1])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,2]*(1-TransMatrix[3,2])),(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3])
    Sigma[9,]<-c(0,0,0,0,0,0,(-TransMatrix[3,1]*TransMatrix[3,3])/stationary[1,3],(-TransMatrix[3,3]*TransMatrix[3,2])/stationary[1,3],(1/stationary[1,3])*(TransMatrix[3,3]*(1-TransMatrix[3,3])))
    
    CSigma <- C%*%Sigma
    
    ev <- eigen(CSigma)
    
    betas <- ev$values
    betas
    
    b1 <- betas[1]
    b2 <- betas[2]
    b3 <- betas[3]
    b4 <- betas[4]
    b5 <- betas[5]
    b6 <- betas[6]
    b7 <- betas[7]
    b8 <- betas[8]
    b9 <- betas[9]
    
    distr <- c()
    for(i in 1:10000){
      
      r.chi <- b1*rchisq(n=10000, df=1) + b2*rchisq(n=10000, df=1) + b3*rchisq(n=10000, df=1)+ b4*rchisq(n=10000, df=1)+ b5*rchisq(n=10000, df=1) + b6*rchisq(n=10000, df=1) + b7*rchisq(n=10000, df=1) + b8*rchisq(n=10000, df=1)+ b9*rchisq(n=10000, df=1)
      
      distr[i]<- r.chi
    }
    
    q.distr <- quantile(distr, probs =seq(0, 1, 0.01))
    
    
    w <- array(NA,dim=c(3,3,1000))
    for (i in 1:1000) {
      
      sample1 <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est1 <- function(X, prob=T){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est1(X=sample1,prob = T)
      
      w[,,i] <-abs(TransMatrix.est-TransMatrix)
    }
    
    sum11 <- matrix(0,3,3)
    for (i in 1:1000) {
      sum11 <- sum11 + w[,,i]
    }
    w <- sum11*(1/1000)
    
    CWKL1 <-c()
    for (i in 1:1000) {
      
      sample <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est <- function(X, prob=T){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est(X=sample,prob = T)
      
      #w <-abs(TransMatrix.est-TransMatrix)
      
      CWKL.1 <- function(TransMatrix, w, TransMatrix.est, stationary){
        sum <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum <- sum + (stationary[1,i]*w[i,j]*(TransMatrix.est[i,j]*log((TransMatrix.est[i,j])/(TransMatrix[i,j]))-TransMatrix.est[i,j]+TransMatrix[i,j]))
          }
        }
        return(sum)
      } 
      
      CWKL.2 <- function(TransMatrix, w, stationary){
        sum1 <- 0
        for (i in 1:length(States)) {
          for (j in 1:length(States)) {
            sum1 <- sum1 + (w[i,j]*stationary[1,i]*TransMatrix.est[i,j])
          }
        }
        return(sum1)
      } 
      
      CWKL.est.1 <- CWKL.1(TransMatrix=TransMatrix, w=w, TransMatrix.est=TransMatrix.est, stationary=stationary)
      
      CWKL.est.2 <- CWKL.2(TransMatrix=TransMatrix, w=w, stationary=stationary)
      
      CWKL.est <- (CWKL.est.1)/(CWKL.est.2) #CWKL divergence
      
      CWKL <- (2*sample_size*CWKL.est)
      
      CWKL1[i] <- CWKL
      
      CWKL2 <- CWKL1[!is.na(CWKL1)]
      
      u <- CWKL2 >= q.distr[96]
      power <- length(u[u==TRUE])/length(u)
    }
    powerv[n] <- power
    powerv<-powerv[!is.na(powerv)]
  }
  if(c==0){
    power3[,1]<-powerv
  }
  
  if(c==1/10){
    power3[,2]<-powerv
  }
  
  if(c==1/50){
    power3[,3]<-powerv
  }
  
  if(c==1/100){
    power3[,4]<-powerv
  }
} 




#Chi-Squared power
distr <- c()
for(i in 1:10000){
  
  r.chi <- rchisq(n=10000, df=(length(States)*(length(States)-1)))
  
  distr[i]<- r.chi
}
q.distr.chi <- quantile(distr, probs =seq(0, 1, 0.01))

power.chi<-matrix(data = c(0,0,0,0,
                           0,0,0,0,
                           0,0,0,0,
                           0,0,0,0,
                           0,0,0,0),byrow = FALSE, nrow = 5, dimnames = list(c("100","500","1000","5000","10000"), c("0","1/10","1/50","1/100")))
for (c in c(0,1/10,1/50,1/100)) {
  powerv <- c()
  for (n in c(100,500,1000,5000,10000)) {
    
    sample_size <- n
    
    c
    
    States <- c("1", "2", "3")
    
    byRow <- TRUE
    
    TransMatrix <- matrix(data = c(1/3, 1/3, 1/3,
                                   2/6, 1/6, 3/6,
                                   (1+6*c)/6, (2+6*c)/6, (3-12*c)/6), byrow = byRow, nrow = 3, 
                          dimnames = list(States, States))
    
    
    
    TransMatrix.alt <- matrix(data = c(1/3, 1/3, 1/3,
                                       2/6, 1/6, 3/6,
                                       1/6, 2/6, 3/6), byrow = byRow, nrow = 3, 
                              dimnames = list(States, States))
    
    mc <- new("markovchain", states = States, byrow = byRow,
              transitionMatrix = TransMatrix.alt, name = "Markov Model")
    
    
    
    values <-c()
    for(l in 1:10000){
      
      sample <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est1 <- function(X, prob){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est1(X=sample,prob = F)
      
      chi.sq <- function(TransMatrix.est, TransMatrix){
        sum <- 0
        for (i in length(States)) {
          for (j in length(States)) {
            sum <- sum + (((TransMatrix.est[i,j] - sum(TransMatrix.est[i,])*TransMatrix[i,j])^2)/(sum(TransMatrix.est[i,])*TransMatrix[i,j]))
          }
        }
        return(sum)
      } #Chi-Square test for markov chains
      
      values[l] <- chi.sq(TransMatrix.est = TransMatrix.est, TransMatrix = TransMatrix)
      
      u1 <- values[!is.na(values)]
      
      u <- u1 >= q.distr.chi[96]
      
      power <- length(u[u==TRUE])/length(u)
    }
    powerv[n] <- power
    powerv<-powerv[!is.na(powerv)]
  }
  if(c==0){
    power.chi[,1]<-powerv
  }
  
  if(c==1/10){
    power.chi[,2]<-powerv
  }
  
  if(c==1/50){
    power.chi[,3]<-powerv
  }
  
  if(c==1/100){
    power.chi[,4]<-powerv
  }
}  




#LRT power 
power.lrt<-matrix(data = c(0,0,0,0,
                           0,0,0,0,
                           0,0,0,0,
                           0,0,0,0,
                           0,0,0,0),byrow = FALSE, nrow = 5, dimnames = list(c("100","500","1000","5000","10000"), c("0","1/10","1/50","1/100")))
for (c in c(0,1/10,1/50,1/100)) {
  powerv <- c()
  for (n in c(100,500,1000,5000,10000)) {
    
    sample_size <- n
    
    c
    
    States <- c("1", "2", "3")
    
    byRow <- TRUE
    
    TransMatrix <- matrix(data = c(1/3, 1/3, 1/3,
                                   2/6, 1/6, 3/6,
                                   (1+6*c)/6, (2+6*c)/6, (3-12*c)/6), byrow = byRow, nrow = 3, 
                          dimnames = list(States, States))
    
    
    
    TransMatrix.alt <- matrix(data = c(1/3, 1/3, 1/3,
                                       2/6, 1/6, 3/6,
                                       1/6, 2/6, 3/6), byrow = byRow, nrow = 3, 
                              dimnames = list(States, States))
    
    mc <- new("markovchain", states = States, byrow = byRow,
              transitionMatrix = TransMatrix.alt, name = "Markov Model")
    
    
    values <-c()
    for(l in 1:10000){
      
      sample <- rmarkovchain(sample_size, object = mc) #Markov chain realization
      
      TransMatrix.est1 <- function(X, prob){
        tt <- table(c(X[-length(X)]), c(X[-1]))
        if(prob) tt <- tt / rowSums(tt)
        tt
      } #Transition Probability estimation
      
      TransMatrix.est <- TransMatrix.est1(X=sample,prob = F)
      
      lrt.markov.order1 <- function(States, sample, TransMatrix){
        nbstates <- length(States)
        
        sample_size <- length(sample)
        
        TransMatrix.est <- function(X, prob=T){
          tt <- table(c(X[-length(X)]), c(X[-1]))
          if(prob) tt <- tt / rowSums(tt)
          tt
        } 
        
        initial <- function(sample){
          init <- table(sample)/length(sample)
          return(init)
        } #initial distribution estimator
        
        proba1<-initial(sample)
        Pest1<-TransMatrix #hypothetical trans matrix
        
        proba2<-initial(sample)
        Pest2<-TransMatrix.est(X=sample,prob = T) #estimated trans matrix
        
        log.like.markov.order1 <- function(seq,alphabet,proba,Pest){ # Attention aux zéros dans Pest
          Nij<-matrix(count(seq,2,alphabet = s2c(paste(States, collapse = ""))), byrow=T, ncol=nbstates)
          lV<- as.numeric(log(proba[which(alphabet==seq[1])])) + sum(Nij*log(Pest))
          return(lV)
        } #log likelihood
        
        L0<- log.like.markov.order1(seq=sample,alphabet=States,proba1,Pest1)
        L1<- log.like.markov.order1(seq=sample,alphabet=States,proba2,Pest2)
        
        K1<-length(States)*(length(States)-1)
        K0<-0
        
        pvalue <- pchisq(-2*(L0-L1),K1-K0,lower.tail=FALSE)
        
        return(pvalue)
      } ##lrt function
      
      values[l] <- lrt.markov.order1(States = States, sample=sample,TransMatrix = TransMatrix)
      
      u1 <- values[!is.na(values)]
      
      u <- u1 < 0.05
      
      power <- length(u[u==TRUE])/length(u)
    }
    powerv[n] <- power
    powerv<-powerv[!is.na(powerv)]
  }
  if(c==0){
    power.lrt[,1]<-powerv
  }
  
  if(c==1/10){
    power.lrt[,2]<-powerv
  }
  
  if(c==1/50){
    power.lrt[,3]<-powerv
  }
  
  if(c==1/100){
    power.lrt[,4]<-powerv
  }
}  

