#' Doubly Sparse Survival Rule Extraction
#'
#' \code{SURVFIT} extracts a doubly sparse (sparse in both number of rules and
#' in number of variables in the rules) survival rule ensemble from survival
#' data
#'
#' @param formula \code{formula}. The model formula specifying time, status and
#'   dependent variables of the form \code{Surv(time, status)~ x1 + x2 + .. }
#' @param data \code{data.frame}. Training data.
#' @param rulelength \code{Integer}. Maximum length of the rule. (Default = 3)
#' @param doubly.sparse \code{Logical} for whether double sparsity required.
#'   (Default = FALSE)
#' @param gamma \code{Numeric} or \code{list}. Hyperparameter (Default = NULL)
#' @param lambda1 \code{Numeric} or \code{list}. Hyperparameter (Default = NULL)
#' @param lambda2  \code{Numeric} or \code{list}. Hyperparameter (Default =
#'   NULL)
#' @param crossvalidate \code{Logical}. Whether crossvalidation to be done to
#'   find hyperparameters. (Default = TRUE)
#' @param nfolds \code{Integer}. Number of cross validation folds. (Default = 5)
#' @param num_toprules \code{Integer}. Number of rules extracted. (Default = 16)
#' @param num_totalrules \code{Integer}. Number of rules considered. (Default =
#'   2000)
#' @param input_rule_list \code{Logical} Whether rule list supplied. (Default =
#'   FALSE)
#' @param rule_list \code{List}. List of supplied rules. (Default = NULL)
#' @param ntree \code{Integer} .Number of trees built
#' @param digit \code{Integer}. Decimal points.
#' @param seed \code{Numeric}. Seed for reproducible experiments.
#' @param nodesize \code{Integer}. (Default = NULL)
#' @param ... Other inputs

#' @return Object of class \code{list} with elements
#'   \item{\code{rules}}{List of top \code{\link{num_toprules}} rules}
#'   \item{\code{all_rules}}{List of all \code{\link{num_totalrules}} rules}
#'   \item{\code{rule_data}}{\code{Data.frame} of rules evaluated over data}
#'   \item{\code{beta}}{Coefficients of \link{all_rules} in the model}
#' @examples
#' ## For ovarian data from survival package.
#' SURVFIT(Surv(futime, fustat) ~ ., data = ovarian)
#'
#'@export
SURVFIT <- function(formula = formula,
                    data = data,
                    rulelength = 3,
                    doubly.sparse = FALSE,
                    gamma = NULL,
                    lambda1 = NULL,
                    lambda2 = NULL,
                    crossvalidate = TRUE,
                    nfolds = 5,
                    num_toprules = 16,
                    num_totalrules = 2000,
                    input_rule_list = FALSE,
                    rule_list = NULL,
                    ntree = 200,
                    digit = 10,
                    seed = NULL,
                    nodesize = NULL,
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


  ## if rules to be evaluated need to be generated from the data
  newdata <- data[,c(yvar.names,xvar.names)]
  remove(data)

  ## this step simply omits the observation with missing values in the dataset (NOTE)
  newdata <- parseMissingData(newdata)




  #### RULE GENERATION ####

  if (input_rule_list == TRUE | !is.null(rule_list)){
    ## if rules to be evaluated are already given in the data
    rule_list = rule_list

  }else{

    rule_list <- get_rules(formula = formula,
                            data = newdata,
                            xvar.names = xvar.names,
                            yvar.names = yvar.names,
                            ntree = ntree,
                            nodesize = NULL, nodedepth = rulelength,
                            digit = digit,
                            seed = NULL,
                            family = "surv",
                            num_totalrules = num_totalrules,
                            ...)
  }


  #### RULE DATA GENERATION ####

  rdata <- generate_ruledata(newdata,rule_list)
  intercept <- rep(1,nrow(rdata))
  rdata <- cbind(data[,yvar.names],intercept,rdata)
  colnames(rdata) <- c("time","status","intercept",sapply(1:length(rule_list),function(x) paste('rule',x,sep="")))

  # rules to rule_names (using colnames of newdata)




  #### OPTIMIZATION SOLUTION AND ALGORITHM ####

  if(is.null(gamma) & is.null(lambda)){

    cv <- cv.lambda_gamma(rdata,nvar)
    gamma <- cv$gamma.star
    lambda <- cv$lambda.star
  }

  if(!doubly.sparse){

    coefs <- solve.qp.cplex(rdata,nvar,gamma,lambda)
    beta <- coefs$beta

    top_ix <- sort(abs(coefs), decreasing = TRUE, index.return = TRUE)$ix -1
    top_rule_ix <- top_ix[1:num_toprules]
    top_rules <- rule_list[top_rule_ix]

    return(list(rules = top_rules, all_rules = rule_list, rdata = rdata,beta = beta))

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
