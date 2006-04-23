#########################################################################
#
# function to solve the Mixed Model Equations
# Part of the R/maanova package
#
# Written by Hao Wu 2003
# Modified May 2004 for new algorithms
#
# Currently method of scoring is used to solve REML and ML
#
# For performance purpose, I will skip all data checking in function
#
#########################################################################


mixed <- 
  function(y, X, Z, XX, XZ, ZZ, Zi, ZiZi, dimZ, s20,
           method=c("noest","MINQE-I","MINQE-UI","ML", "REML"), maxiter=100)
{
  
 # define macro
  EPS <- 1e-4
  
  method <- match.arg(method)
  # local variables
  yy <- sum(y*y)
  Xy <- t(X) %*% y
  Zy <- t(Z) %*% y
  a <- c(Xy, Zy)
  n <- length(y)
  k <- dim(X)[2]
  m <- dim(Z)[2]
  rx <- matrank(XX)
  r <- length(s20)
  crit <- 1
  loops <- 0
  
  #===============================================
  # calculate b and u using initial sigma
  #===============================================
  result <- solveMME(s20, dimZ, XX, XZ, ZZ, a)

  #======================================================================
  # "noest" : No estimation of variance components 
  # we are done!
  #======================================================================
  if(method=="noest") {
    # we are done, return
    result$s2 <- s20
    result$loops <- loops
    return(result)
  }

  #=====================================================
  # do one step of MINQUE. The result will be used as
  # the initial values for REML and ML
  #=====================================================
  loops <- loops + 1
  # make H matrix and q vector
  if(method == "ML") #ML
  Hq <- makeHq(s20, y, X, Z, Zi, ZiZi, dim, result$b, "MINQE-I")
  else if(method == "REML") #REML
    Hq <- makeHq(s20, y, X, Z, Zi, ZiZi, dim, result$b, "MINQE-UI")
  else
    Hq <- makeHq(s20, y, X, Z, Zi, ZiZi, dim, result$b, method)

  #================================
  # re-estimated s20 and solve MME
  #================================
  s20 <- as.vector(pinv(Hq$H) %*% Hq$q)
  result <- solveMME(s20, dimZ, XX, XZ, ZZ, a)
  # if method is MINQUE, stop here
  if (method=="MINQE-I" | method=="MINQE-UI") {
    result$s2 <- s20
    result$loops <- loops
    return(result)
  }

  #========================================
  # solve MME iteratively for ML and REML
  #========================================
  # loop until converge
  while ( (crit > EPS) && (loops<maxiter) ) {
    loops <- loops+1
    sigaux <- s20
    # make H and q 
    Hq <- makeHq(s20, y, X, Z, Zi, ZiZi, dim, result$b, method)
    if(method=="ML") { # ML
      # update s2 and b 
      tmp <- pinv(Hq$H) %*% Hq$q
      result$b <-  result$b + tmp[1:k]
      s20 <- s20 + tmp[(k+1):length(tmp)]
    }
    else { # REML
      # update s20
      s20 <- s20 + as.vector(pinv(Hq$H) %*% Hq$q)
    }
    # s20 cannot be negative
    s20[s20<0] <- 1e-10
    # calculate crit
    crit <- norm(sigaux-s20)
  }

  # solve MME using converged s2
  result <- solveMME(s20, dimZ, XX, XZ, ZZ, a);
  result$s2 <- s20
  result$loops <- loops

  result
}


#=====================================================================
#This is the function to make H and q matrices for different methods
# For MINQUE, H*s20 = q
# For ML/REML, H and q are Hessian and gradient
#=====================================================================
makeHq <- function(s20, y, X, Z, Zi, ZiZi, dim, b, method)
{
  # some local variables
  # number of observations
  N <- length(y)
  # number of random variables
  nr <- length(s20)
  
  # make V matrix, V=var(y), V=sum(Zi*Zi'*sigma_i^2) + sigma_e^2*I(N)
  V <- matrix(0, N, N)
  for(i in 1:nr)
    V <- V + ZiZi[[i]]*s20[i]
  invV <- solve(V)
  
  # init output
  H <- matrix(0, nr, nr)
  q <- rep(0, nr)
  
  #===================================
  # for ML or MINQUE0
  #===================================
  if(method=="ML" | method=="MINQE-I") {
    # First make H
    # for sigma
    Is <- matrix(0, nr,nr)
    for (i in 1:(nr-1) ) {
      tmp <- t(Zi[[i]]) %*% invV %*% Zi[[i]]
      tmp <- sum(tmp^2)
      Is[i,i] <- tmp
      for(j in (i+1):nr) {
        tmp <- t(Zi[[j]]) %*% invV %*% Zi[[i]]
        tmp <- sum(tmp^2)
        Is[i,j] <- tmp
        Is[j,i] <- tmp
      }
    }
    # last element of Is
    Is[nr,nr] <- sum(invV^2)
    
    if(method == "MINQE-I") # for MINQUE-I
      H <- Is
    else { # for ML
      # for b
      Ib <- t(X) %*% invV %*% X
      H <- blkdiag(Ib, 0.5*Is)
    }
    
    # make q
    tmp <- y - X%*%b
    if(method=="ML") { # for ML
      # for b
      lb <- t(X) %*% invV %*% tmp
      # for sigma
      ls <- rep(0,nr)
      for(i in 1:nr) {
        tmp2 <- t(tmp) %*% invV %*% ZiZi[[i]] %*% invV %*% tmp
        tmp3 <- invV %*% ZiZi[[i]]
        ls[i] <- (tmp2-sum(diag(tmp3)))/2
      }
      q <- c(lb, ls)
    }
    else { # For MINQUE
      for(i in 1:nr)
        q[i] <- t(tmp) %*% invV %*% ZiZi[[i]] %*% invV %*% tmp
    }
  }
  

  #===================================
  # for MINQUE(U,I) or REML
  #===================================
  else if(method=="REML" | method=="MINQE-UI") {
    # make P matrix
    invVX <- invV %*% X
    P <- invV - invVX %*% pinv(t(X)%*%invVX) %*% t(invVX)
    # make H
    for(i in 1:(nr-1)) {
      tmp <- t(Zi[[i]]) %*% P %*% Zi[[i]]
      tmp <- sum(tmp^2)
      H[i,i] <- tmp
      for(j in (i+1):nr) {
        tmp <- t(Zi[[j]]) %*% P %*% Zi[[i]]
        tmp <- sum(tmp^2)
        H[i,j] <- tmp
        H[j,i] <- tmp
      }
    }
    # last element of H
    H[nr,nr] <- sum(P^2)
    
    if(method == "REML") # for REML
      H <- H * 0.5
    
    # make q
    if(method == "REML") { # for REML
      for(i in 1:nr) { # loop for random terms
        tmp <- P %*% ZiZi[[i]]
        tmp2 <- t(y) %*% tmp %*% P %*% y
        q[i] <- (tmp2-sum(diag(tmp)))/2
      }
    }
    else { # for MINQUE(I,U)
      for(i in 1:nr) 
        q[i] <- t(y) %*% P %*% ZiZi[[i]] %*% P %*% y
    }
  }

  # return
  list(H=H, q=q)
}


#===============================================
# function to solve MME given sigma^2
#===============================================
solveMME <- function(s20, dim, XX, XZ, ZZ, a)
{
  k <- dim(XZ)[1]
  m <- dim(XZ)[2]
  errors2 <- s20[length(s20)]
  # make D matrix
  D <- makeD(s20, dim)
  # make V matrix
  V <- diag(errors2,m,m) + ZZ%*%D
  # make A matrix
  A <- rbind(cbind(XX,XZ%*%D), cbind(t(XZ),V))
  # calculate estimates
  bb <- pinv(A)%*%a
  b <- bb[1:k]
  v <- bb[(k+1):(k+m)]
  u <- D%*%v

  list(b=b, u=as.vector(u))
}


#==============================
# sub functions for mixed
# make D matrix
#==============================
makeD <- function(s20, dimZ)
{
  r <- length(dimZ)
  d <- NULL
  for(i in 1:r) 
    d <- c(d, s20[i]*rep(1, dimZ[i]))

  result <- diag(d)

  result
}


#===================================================
# function to make Zi and Zi*Zi' for each part of Z
#===================================================
makeZiZi <- function(Z, dimZ)
{
  # number of observations
  N <- dim(Z)[1]
  # number of random variables (exclude error)
  nr <- length(dimZ);

  # divide Z into several parts for each random term
  Zi <- NULL
  tmp <- c(0, cumsum(dimZ))
  for (i in 1:nr) {
    colstart <- tmp[i] + 1
    colend <- tmp[i+1]
    Zi[[i]] <- Z[,colstart:colend]
  }

  # calculate Zi*Zi' for each Zi
  ZiZi <- NULL
  for (i in 1:nr)
    ZiZi[[i]] <- Zi[[i]] %*% t(Zi[[i]])

  Zi[[nr+1]] <- diag(1, N, N)

  ZiZi[[nr+1]] <- diag(1, N, N)

  list(Zi=Zi, ZiZi=ZiZi)

}
