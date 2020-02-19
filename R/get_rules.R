get_rules <- function(formula,
                      data,
                      xvar.names,
                      yvar.names,
                      ntree,
                      nodesize,
                      nodedepth,
                      digits,
                      seed,
                      family = "surv",
                      num_totalrules,
                      ...){


  #### GENERATE THE RULE LIST ####

  ## build a ranger survival random forest
  rf.ranger = ranger::ranger(formula, data, num.trees = ntree,
                             max.depth = nodedepth, respect.unordered.factors = 'partition')
  ## get list of all trees from the ranger forest
  tree.list <- Ranger2List(rf.ranger)

  ## initialize list of all rules
  rule_list <- c()
  ## get rules of depth 1:nodedepth from the tree
  for(i in 1:nodedepth){
    rules <- unique(extractRules(tree.list,data, ntree = ntree, maxdepth = i, digits=digits))
    ## if node depth is 1 | retain only every other rule
    if (i == 1){
      # nrules = length(rules)
      # idx <- c(1:(nrules%/%2))*2
      # rules <- rules[idx]
      if(length(rules) > 300){
        rules <- sample(rules, 300)
      }
    } else if (i == 2){
      if(length(rules) > 1000){
        rules <- sample(rules, 1000)
      }
    } else{
      if(length(rules) > 700){
        rules <- sample(rules, 700)
      }
    }
    rule_list <- c(rule_list,rules)
  }
  ## select 2000 rules randomly from all the rules extracted (simply for convenience (NOTE))
  if(length(rule_list) > num_totalrules){
    rule_list <- sample(rule_list, num_totalrules, replace = FALSE)
    paste0('A random subset of ', num_totalrules, ' rules was selected')
  }

  #### GENERATE RULES WITH VARIABLES SUBSTITUTED ####

  rule_list
}
