survfit <- function(formula = formula,
                    data = data,
                    ntree = 200,
                    nodesize = null, nodedepth = 3,
                    input_rule_list = FALSE,
                    rule_list = NULL,
                    digit = 10,
                    seed = NULL,
                    family = "surv",
                    doubly.sparse = FALSE,
                    gamma = NULL, lambda1 = NULL, lambda2 = NULL,
                    crossvalidate = TRUE,
                    nfolds = 5,
                    num_toprules = 16,
                    num_totalrules = 2000,
                    ...){

  #### PRELIMNARY TESTING AND FORMATTING ####

  if(!is.null(seed)){
    set.seed(seed)
  }
  ## survival
  formulaPrelim <- parseFormula(formula, data)
  ## data cannot be missing
  if (missing(data)) stop("data is missing")
  ## conduct preliminary formula validation
  if (missing(formula) | (!missing(formula) && is.null(formula))) {
    stop("formula is missing or null")
  }
  # save the call/formula for the return object
  my.call <- match.call()
  my.call$formula <- eval(formula)

  ## finalize the formula based on the pre-processed data
  formulaDetail <- finalizeFormula(formulaPrelim, data)

  ## coherence checks on option parameters
  ntree <- round(ntree)
  if (ntree < 1) stop("Invalid choice of 'ntree'.  Cannot be less than 1.")
  nodedepth = round(nodedepth)
  if (nodedepth < 1) stop("Invalid choice of 'nodedepth'.  Cannot be less than 1.")

  if (!is.logical(doubly.sparse)) stop("Invalid choice of 'var.sparse'. Must be logical. ")
  ## initialize the seed
  # seed <- get.seed(seed)

  ## save the family for convenient access
  family <- formulaDetail$family

  ## save the names for convenient access
  xvar.names <- formulaDetail$xvar.names
  yvar.names <- formulaDetail$yvar.names
  nvar <- length(xvar.names)
  ## reality check on x and y
  ## are there any x-variables?
  if (length(xvar.names) == 0) {
    stop("something seems wrong: your formula did not define any x-variables")
  }
  ## .. are there any y-variables?
  if (length(yvar.names) == 0) {
    stop("something seems wrong: your formula did not define any y-variables")
  }

  if (family != 'surv'){
    stop("family should be 'surv'")
  }

  #### RULE GENERATION ####

  if (input_rule_list == TRUE | !is.null(rule_list)){
    ## if rules to be evaluated are already given in the data
    rule_list = rule_list

  }else{

    rule_list <- get_rules(formula = formula,
                            data = data,
                            xvar.names = xvar.names,
                            yvar.names = yvar.names,
                            ntree = ntree,
                            nodesize = NULL, nodedepth = nodedepth,
                            digit = digit,
                            seed = NULL,
                            family = "surv",
                            num_totalrules = num_totalrules,
                            ...)
  }


  #### RULE DATA GENERATION ####

  rdata <- generate_ruledata(data,rule_list)
  intercept <- rep(1,nrow(rdata))
  rdata <- cbind(data[,yvar.names],intercept,rdata)
  colnames(rdata) <- c("time","status","intercept",sapply(1:length(rule_list),function(x) paste('rule',x,sep="")))


  #### OPTIMIZATION SOLUTION AND ALGORITHM ####

  if(is.null(gamma) & is.null(lambda)){

    cv <- cv.lambda_gamma(rdata,nvar)
    gamma <- cv$gamma.star
    lambda <- cv$lambda.star
  }

  if(!doubly.sparse){

    coefs <- solve.qp.cplex(rdata,nvar,gamma,lambda)
    beta <- coefs$beta
    return(list(data = data, rule_list = rule_list, rdata = rdata,beta = beta))

    }else {

    groups <- extractgroups(rule_list, nvar, family)
    # adjust groups:
    for(i in 1:nvar){
      groups[[i]] <- (sqrt(sum(groups[[i]] == 1)))*groups[[i]]
    }

    if(is.null(lambda1) & is.null(lambda2)){

      cv <- cv.lambda12(rdata, nvar, nfolds, gamma,groups)
      lambda1 = cv$lambda1
      lambda2 = cv$lambda2
    }

    coefs <- Rcplex::solve.socp.cplex(rdata,nvar,gamma,lambda1,lambda2,groups=groups)
    beta <- coefs$beta
    return(list(data = data, rule_list = rule_list, rdata = rdata,beta = beta))

  }
}
