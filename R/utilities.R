parseFormula <- function(f,data){
  ## confirm coherency of the formula
  if (!inherits(f, "formula")) {
    stop("'formula' is not a formula object.")
  }
  if (is.null(data)) {
    stop("'data' is missing.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  ## pull the family and y-variable names
  fmly <- all.names(f, max.names = 1e7)[2]
  all.names <- all.vars(f, max.names = 1e7)
  yvar.names <- all.vars(formula(paste(as.character(f)[2], "~ .")), max.names = 1e7)
  yvar.names <- yvar.names[-length(yvar.names)]
  ## survival forests
  if ((fmly == "Surv")) {
    if (sum(is.element(yvar.names, names(data))) != 2) {
      stop("Survival formula incorrectly specified.")
    }
    family <- "surv"
  }else {
    ## must be survival formula
    stop('Survival formual is not specified.')
  }
  ## done: return the goodies
  return (list(all.names=all.names, family=family, yvar.names=yvar.names))
}

finalizeFormula <- function(formula.obj, data) {
  ## parse the formula object
  yvar.names <- formula.obj$yvar.names
  all.names  <- formula.obj$all.names
  index      <- length(yvar.names)
  fmly       <- formula.obj$family
  ## total number of variables should exceed number of yvars
  if (length(all.names) <= index) {
    stop("formula is misspecified: total number of variables does not exceed total number of y-variables")
  }
  ## extract the xvar names
  if (all.names[index + 1] == ".") {
    if(index == 0) {
      xvar.names <- names(data)
    }
    else {
      xvar.names <- names(data)[!is.element(names(data), all.names[1:index])]
    }
  }else {
    if(index == 0) {
      xvar.names <- all.names
    }
    else {
      xvar.names <- all.names[-c(1:index)]
    }
    not.specified <- !is.element(xvar.names, names(data))
    if (sum(not.specified) > 0) {
      stop("formula is misspecified, object ", xvar.names[not.specified], " not found")
    }
  }
  ## return the goodies
  return (list(family=fmly, yvar.names=yvar.names, xvar.names=xvar.names))
}

parseMissingData <- function(data) {
  ## if impute, use imputation techniques to fill in the missing values
  return(na.omit(data))
}

generate_ruledata <- function(data, rules){
  if(length(rules) == 0) stop("rules must be of length >= 1")
    return(NULL)
  X <- data[,3:ncol(data)]
  expr <- parse(text = paste0("cbind(", paste0(rules, collapse = ", "), ")"))
  x <- eval(expr, data)
  colnames(x) <- names(rules)
  x <- ifelse(x == TRUE,1,0)
  return(x)
}

extractgroups <- function(rule_list, nvar, family = "surv"){

  ## function to return a list mapping variables
  ## to rules involving the variables

  ## n is total number of rules
  n = length(rule_list)

  ## initiate the list of groups
  groups <- list()
  # nvar is number of variables in the data
  for(i in 1:nvar){
    groups[[i]] <- rep(0,n)
  }

  for(i in 1:n){
    rule = rule_list[i]
    r <- strsplit(rule, '&')[[1]]
    ## rn is the number of variables in the rule[i]
    rn <- length(r)
    ## get variables in the rules
    vars <- c()
    for(j in 1:rn){
      r. <- r[j]
      z = strsplit(r., "]")[[1]][1]
      vars <- c(vars,as.numeric(substr(trimws(z),4,10)))
    }
    vars <- unique(vars)
    ## set the group corrsponding to this variable to be 1
    for(v in vars){
      groups[[v]][i] = 1
    }
  }
  groups
}

norm_vec <- function(x) sqrt(sum(x^2))

# extractgroups(rule_list[1:5],4,"surv")
mae <- function(a,b) sum(abs(a-b))




