cv.gamma <- function(rdata,nvar, nfolds=5, lambda1 = 0){
  gamma.values <- c(0,1e-06,5e-06, 1e-05, 5e-05, 1e-04, 5e-04, 1e-03, 5e-03, 1e-02,1e-01,1,2)
  errors <- rep(0,length(gamma.values))
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  for(i in 1:length(gamma.values)){
    gamma = gamma.values[i]
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma, lambda1 = lambda1)
      beta <- coefs$beta
      pred = as.matrix(rdata.test[,3:ncol(rdata)])%*%beta
      error = mae(rdata.test$time,pred)
      errors[i] <- errors[i]+error
    }
  }
  index <- which(errors == min(errors))
  gamma.star <- gamma.values[index]
  res <- list()
  res$gamma.star <- gamma.star
  res$gamma.values <- gamma.values
  res$error.values <- errors/5
  res
}

cv.lambda <- function(rdata, nvar, gamma, nfolds=5){
  lambda.values <- c(0,1e-08,1e-06,1e-04,1e-03,1e-02, 1e-01,1,2,5,10,20)
  errors <- rep(0,length(lambda.values))
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  for(i in 1:length(lambda.values)){
    lambda = lambda.values[i]
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma, lambda1 = lambda)
      beta <- coefs$beta
      pred = as.matrix(rdata.test[,3:ncol(rdata)])%*%beta
      error = mae(rdata.test$time,pred)
      errors[i] <- errors[i]+error
    }
  }
  index <- which(errors == min(errors))
  lambda.star <- lambda.values[index]
  res <- list()
  res$lambda.star <- lambda.star
  res$lambda.values <- lambda.values
  res$error.values <- errors/5
  res
}

cv.lambda_gamma <- function(rdata, nvar, nfolds = 4){

  gamma.values <- c(0,1e-8,1e-06,5e-06, 1e-05, 5e-05,
                    1e-04, 5e-04, 1e-03, 5e-03, 1e-02,1e-01,1,2)
  lambda.values <- c(0,1e-08,1e-06,1e-04,1e-03
                     ,1e-02, 1e-01,1,2,5,10,20)
  grid <- expand.grid(gamma.values,lambda.values)
  errors <- rep(0,nrow(grid))

  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)

  for (i in 1:nrow(grid)) {
    gamma = grid[i,1]
    lambda = grid[i,2]
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma, lambda1 = lambda, trace = 0)
      beta <- coefs$beta
      pred = as.matrix(rdata.test[,3:ncol(rdata)])%*%beta
      error = mae(rdata.test$time,pred)
      errors[i] <- errors[i]+error
    }
  }
  index <- which(errors == min(errors))
  gamma.star <- grid[index,1]
  lambda.star <- grid[index,2]
  res <- list()
  res$gamma.star <- gamma.star
  res$gamma.values <- gamma.values
  res$lambda.star <- lambda.star
  res$lambda.values <- lambda.values
  res$grid <- grid
  res$error.values <- errors/nfolds
  res

}

cv.lambda12 <- function(rdata, nvar, nfolds = 4,gamma,groups){

  lambda.values <- c(0,1e-06,1e-04,1e-03
                     ,1e-02, 1e-01,1,2,5,10,20,30)
  grid <- expand.grid(lambda.values,lambda.values)
  errors <- rep(0,nrow(grid))

  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)

  for (i in 1:nrow(grid)) {
    lambda1 = grid[i,1]
    lambda2 = grid[i,2]
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      if(nrow(rdata.train) == 0 | nrow(rdata.test) == 0) next
      coefs <- solve.socp.cplex(rdata = rdata.train, nvar = nvar, gamma = gamma,
                                lambda1 = lambda1, lambda2 = lambda2,
                                groups = groups, trace = 0)
      beta <- coefs$beta
      pred = as.matrix(rdata.test[,3:ncol(rdata)])%*%beta
      error = mae(rdata.test$time,pred)
      errors[i] <- errors[i]+error
    }
  }
  index <- which(errors == min(errors))
  lambda1.star <- grid[index,1]
  lambda2.star <- grid[index,2]
  res <- list()
  res$lambda1.star <- lambda1.star
  res$lambda2.star <- lambda2.star
  res$grid <- grid
  res$lambda1_values <- lambda.values
  res$lambda2_values <- lambda.values
  res$error.values <- errors/nfolds
  res

}


