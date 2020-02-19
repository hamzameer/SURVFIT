cv.gamma <- function(rdata, nfolds=5, gamma.values = NULL, lambda1 = 0, trace = 0){
  if(is.null(gamma.values)){
  gamma.values <- c(0,1e-06,5e-06, 1e-05, 5e-05, 1e-04, 5e-04, 1e-03, 5e-03, 1e-02,1e-01,1,2)}
  errors <- rep(0,length(gamma.values))
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  lgv = length(gamma.values)
  for(i in 1:l){
    gamma = gamma.values[i]
    cat(paste0(i,'/',lgv,': CV for gamma=',gamma,'\n'))
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma, lambda1 = lambda1, trace = trace)
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

cv.lambda <- function(rdata, gamma.star, lambda.values = NULL, nfolds=5, trace = 0){
  if(is.null(lambda.values)){
    lambda.values <- c(0,1e-08,1e-06,1e-04,1e-03,1e-02, 1e-01,1,2,5,10,20)}

  errors <- rep(0,length(lambda.values))
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  llv = length(lambda.values)
  for(i in 1:length(lambda.values)){
    lambda = lambda.values[i]
    cat(paste0(i,'/',llv,': CV for lambda=',lambda,'\n'))
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma.star, lambda1 = lambda, trace = trace)
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

cv.lambda2 <- function(rdata, nvar, gamma, lambda1 = lambda1.star, groups, lambda2.values = NULL, nfolds=5, trace = 0){
  if(is.null(lambda2.values)){
    lambda2.values <- c(0,1e-08,1e-06,1e-04,1e-03,1e-02, 1e-01,1,2,5,10,20)}

  errors <- rep(0,length(lambda2.values))
  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  llv = length(lambda2.values)
  for(i in 1:llv){
    lambda2 = lambda2.values[i]
    cat(paste0(i,'/',llv,': CV for lambda2=',lambda2,'\n'))
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
                                groups = groups, trace = trace)
      beta <- coefs$beta
      pred = as.matrix(rdata.test[,3:ncol(rdata)])%*%beta
      error = mae(rdata.test$time,pred)
      errors[i] <- errors[i]+error
    }
  }
  index <- which(errors == min(errors))
  lambda2.star <- lambda2.values[index]
  res <- list()
  res$lambda2.star <- lambda2.star
  res$lambda2.values <- lambda2.values
  res$error.values <- errors/5
  res
}


cv.lambda_gamma <- function(rdata, nfolds = 4,gamma.values = NULL, lambda.values = NULL,trace = 0, max.grid = 75){

  if(is.null(gamma.values)){
    gamma.values <- c(0,1e-8,1e-06,5e-06, 1e-05, 5e-05,
                    1e-04, 5e-04, 1e-03, 5e-03, 1e-02,1e-01,1,2)}
  if(is.null(lambda.values)){
    lambda.values <- c(0,1e-08,1e-06,1e-04,1e-03
                     ,1e-02, 1e-01,1,2,5,10,20)}

  grid <- expand.grid(gamma.values,lambda.values)
  if (nrow(grid) > max.grid){
    grid <- grid[sample(1:nrow(grid),max.grid),]
  }
  lg = nrow(grid)
  errors <- rep(0,lg)

  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  for (i in 1:lg) {
    gamma = grid[i,1]
    lambda = grid[i,2]
    cat(paste0(i,'/',lg,': CV for gamma=',gamma,', lamba = ',lambda,'\n'))
    for(k in 1:nfolds){
      #test
      test <- ((k-1)*r%/%nfolds+1):(k*r%/%nfolds)
      rdata.test <- rdata.complete[test,]
      #train
      rdata.train <- rbind(rdata.complete[-test,],rdata.cens)
      #model
      coefs <- solve.qp.cplex(rdata.train, gamma = gamma, lambda1 = lambda, trace = trace)
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

cv.lambda12 <- function(rdata, nvar, lambda1.values = NULL, lambda2.values = NULL,nfolds = 4,gamma = 0,groups, trace = 0, max.grid = 75){

  if(is.null(lambda1.values)){
    lambda1.values <- c(0,1e-06,1e-04,1e-03
                     ,1e-02, 1e-01,1,2,5,10,20,30)}

  if(is.null(lambda2.values)){
    lambda2.values <- c(0,1e-06,1e-04,1e-03
                        ,1e-02, 1e-01,1,2,5,10,20,30)}

  grid <- expand.grid(lambda1.values,lambda2.values)
  if (nrow(grid) > max.grid){
    grid <- grid[sample(1:nrow(grid),max.grid),]
  }
  lv = nrow(grid)
  errors <- rep(0,lv)

  rdata.complete <- rdata[rdata$status == 1,]
  rdata.cens <- rdata[rdata$status == 0,]
  r = nrow(rdata.complete)
  for (i in 1:lv) {
    lambda1 = grid[i,1]
    lambda2 = grid[i,2]
    cat(paste0(i,'/',lv,': CV for lambda1=',lambda1,', lamba2 = ',lambda2,'\n'))
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
                                groups = groups, trace = trace)
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
  res$lambda1_values <- lambda1.values
  res$lambda2_values <- lambda2.values
  res$error.values <- errors/nfolds
  res

}


