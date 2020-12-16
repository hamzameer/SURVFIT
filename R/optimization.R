### SOCP OPTIMIZATION ###

solve.socp.cplex <- function(rdata, nvar, gamma, lambda1, lambda2, groups = NULL, trace=trace){

  ## complete observations
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  ## censored observations

  # number of beta variables
  nb = 2*(ncol(rdata.complete)-2)
  # number of z variables
  nz = nrow(rdata.cens)
  # number of variables # nvar

  A = rdata.complete[,3:ncol(rdata.complete)]
  Atild = t(as.matrix(A))%*%as.matrix(A)
  t = rdata.complete$time
  B = as.matrix(unname(rdata.cens[,3:ncol(rdata)]))
  bvec = -as.matrix(rdata.cens$time)

  Qmat = 2*rbind(cbind(rbind(cbind(Atild, -Atild, matrix(0,nb/2,nz)),
                             cbind(-Atild, Atild, matrix(0,nb/2,nz)),
                             cbind(matrix(0,nz,nb), gamma*diag(nz))),
                       matrix(0,nb+nz,nvar)),matrix(0,nvar,nb+nz+nvar))


  Amat = cbind(-B,B, diag(nz),matrix(0,nz,nvar))
  # rm(rdata, rdata.complete, rdata.cens, Atild, B)


  ## with quadratic constraints and linear part :)
  QC <- list()
  # QC$QC$L <- NULL
  QC$dir <- rep("L",nvar)
  QC$b <- rep(0,nvar)
  Qi_def = matrix(0,nb+nz+nvar,nb+nz+nvar)
  for(i in 1:nvar){
    ci <- c(0,groups[[i]])
    Qi <- Qi_def
    Qi[1:(nb/2),1:(nb/2)] <- diag(ci)
    Qi[((nb/2)+1):nb,((nb/2)+1):nb] <- diag(ci)
    Qi[1:(nb/2),((nb/2)+1):nb] <- -diag(ci)
    Qi[((nb/2)+1):nb,1:(nb/2)] <- -diag(ci)
    Qi[nb+nz+i,nb+nz+i] <- -1
    QC$QC$Q[[i]] <- Qi
    QC$QC$L[[i]] <- rep(0,nb+nz+nvar)
  }

  ## lower bounds
  lb = c(rep(0,nb),rep(-Inf,nz),rep(0,nvar))
  # upper bounds
  ub = c(rep(Inf,nb),rep(0,nz),rep(Inf,nvar))
  # rm(Qi_def)

  cy = c()
  for(i in 1:nvar){
    cy = c(cy,lambda2*(1/sqrt(sum(groups[[i]]))))
  }
  cy[which(is.nan(cy))] <- 0
  cy[which(cy == Inf)] <- 0

  cvec = c(c(rep(lambda1,nb/2) - 2*t%*%as.matrix(A),
             rep(lambda1,nb/2) + 2*t%*%as.matrix(A), rep(0,nz)),cy)


  res <- Rcplex::Rcplex_solve_QCP(cvec, Amat, bvec, Qmat = Qmat, QC=QC,
                          lb = lb, ub = ub, sense = "L", objsense = c("min"),
                          vtype= NULL, control = list(preind = 1,trace=trace))

  beta_coefs = res$xopt[1:(nb/2)]-res$xopt[((nb/2)+1):nb]
  z_coefs = res$xopt[(nb+1):(nb+nz)]
  return(list(beta = beta_coefs, z = z_coefs))
}


### QUADRATIC PROGRAMMING OPTIMIZATION ###
solve.qp.cplex <- function(rdata, nvar, gamma, lambda1, trace = trace){

  ## complete observations
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  ## censored observations

  # number of beta variables
  nb = 2*(ncol(rdata.complete)-2)
  # number of z variables
  nz = nrow(rdata.cens)
  # number of variables # nvar

  A = rdata.complete[,3:ncol(rdata.complete)]
  Atild = t(as.matrix(A))%*%as.matrix(A)
  t =rdata.complete$time
  B = as.matrix(unname(rdata.cens[,3:ncol(rdata)]))
  bvec = -as.matrix(rdata.cens$time)

  Qmat = 0.5*rbind(cbind(Atild, -Atild, matrix(0,nb/2,nz)),
                   cbind(-Atild, Atild, matrix(0,nb/2,nz)),
                   cbind(matrix(0,nz,nb), gamma*diag(nz)))
  cvec = c(rep(lambda1,nb/2) - 2*t%*%as.matrix(A),
           rep(lambda1,nb/2) + 2*t%*%as.matrix(A), rep(0,nz))

  Amat = cbind(B,-B, -diag(nz))
  lb = c(rep(0,nb),rep(-Inf,nz))
  ub = c(rep(Inf,nb),rep(1,nz))
  # quadratic program with no quadratic constraints
  res <- Rcplex::Rcplex(cvec, Amat, bvec, Qmat = Qmat,lb = lb,
                ub = ub, sense = "L", control = list(trace = trace), objsense = c("min"))

  beta_coefs = res$xopt[1:(nb/2)]-res$xopt[((nb/2)+1):nb]
  z_coefs = res$xopt[(nb+1):(nb+nz)]

  return(list(beta = beta_coefs, z = z_coefs))
}

