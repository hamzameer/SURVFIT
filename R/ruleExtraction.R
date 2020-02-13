Ranger2List <- function(rf_ranger)
{
  formatRanger <- function(tree){
    rownames(tree) <- 1:nrow(tree)
    tree$`status` <- ifelse(tree$`terminal`==TRUE,-1,1)
    tree$`left daughter` <- tree$`leftChild` + 1
    tree$`right daughter` <- tree$`rightChild` + 1
    tree$`split var` <- tree$`splitvarID`+1
    tree$`split point` <- tree$`splitval`
    tree$`prediction` <- tree$`prediction`
    tree <- tree[,c("left daughter","right daughter","split var","split point","status")]
    tree <- as.data.frame(tree)
    return(tree)
  }
  treeList <- NULL
  treeList$ntree <- rf_ranger$num.trees
  treeList$list <- vector("list",rf_ranger$num.trees)
  for(i in 1:rf_ranger$num.trees){
    treeList$list[[i]] <- formatRanger( ranger::treeInfo(rf_ranger, tree = i) )
  }
  return(treeList)
}

extractRules <-function(treeList,X,ntree=100,maxdepth=6,random=FALSE,digits=NULL){
 if(is.numeric(digits)) digits <- as.integer(abs(digits))

 levelX = list()
 for(iX in 3:ncol(X))
 levelX <- c(levelX,list(levels(X[,iX])))
 # X <- NULL; target <- NULL
 ntree=min(treeList$ntree,ntree)
 allRulesList = list()
 for(iTree in 1:ntree){
 if(random==TRUE){max_length = sample(1:maxdepth,1,replace=FALSE)}else{
 max_length = maxdepth}
 rule = list(); count = 0; rowIx = 1;
 # tree = getTree(rf,iTree,labelVar=FALSE)
 tree <- treeList$list[[iTree]]
 if(nrow(tree)<=1) next # skip if there is no split
 ruleSet = vector("list", length(which(tree[,"status"]==-1)))
 res = treeVisit(tree,rowIx = rowIx,count,ruleSet,rule,levelX,length=0,max_length=max_length,digits=digits)
 allRulesList = c(allRulesList, res$ruleSet)
 }
allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
cat(paste(length(allRulesList)," rules (length<=",
max_length, ") were extracted from the first ", ntree," trees.","\n",sep=""))

rulesExec <- ruleList2Exec(X,allRulesList)
return(rulesExec)
}

treeVisit <-function(tree,rowIx,count,ruleSet,rule,levelX,length,max_length,digits=NULL)
{
  res <-list()
  if( tree[rowIx,"status"] == -1 | length == max_length ){
    count = count + 1
    ruleSet[[count]] = rule
    return(list(ruleSet = ruleSet, count=count))
  }
  xIx <- tree[rowIx,"split var"]
  xValue <- tree[rowIx,"split point"]
  if(is.integer(digits)) xValue <- round(as.numeric(tree[rowIx,"split point"]), digits)

  if(is.null(levelX[[xIx]])){
    lValue <- paste("X[,",xIx, "]<=",xValue,sep="")
    rValue <- paste("X[,",xIx, "]>",xValue,sep="")
  }else{
    xValue<- which(as.integer(intToBits(as.integer(xValue)))>0)
    lValue <- levelX[[xIx]][xValue]
    rValue <- setdiff(levelX[[xIx]],lValue)
    #   lValue <- paste("X[,",xIx, "]%in% '",lValue,"'",sep="")
    #   rValue <- paste("X[,",xIx, "]%in% '",rValue,"'",sep="")
  }
  xValue <- NULL
  ruleleft <- rule
  if(length(ruleleft)==0)
  {
    ruleleft[[as.character(xIx)]] <- lValue
  }else{
    if(as.character(xIx) %in% ls(ruleleft)) {
      if(!is.null(levelX[[xIx]])){
        lValue <- intersect(ruleleft[[as.character(xIx)]],lValue)
        ruleleft[[as.character(xIx)]] <- lValue
      }else{
        ruleleft[[as.character(xIx)]] <- paste(ruleleft[[as.character(xIx)]], "&", lValue)
      }
    }else{
      ruleleft[[as.character(xIx)]] <- lValue
    }
  }

  #thisItem = paste("X[,",xIx, "] %in% ", nxValue, sep="")
  ruleright <- rule
  if(length(ruleright)==0){
    ruleright[[as.character(xIx)]] <- rValue
  }else{
    if(as.character(xIx) %in% ls(ruleright)) {
      if(!is.null(levelX[[xIx]])){
        rValue <- intersect(ruleright[[as.character(xIx)]],rValue)
        ruleright[[as.character(xIx)]] <- rValue
      }else{
        ruleright[[as.character(xIx)]] <- paste(ruleright[[as.character(xIx)]], "&", rValue)
      }
    }else{
      ruleright[[as.character(xIx)]] <- rValue
    }
  }

  thisList = treeVisit(tree, tree[rowIx,"left daughter"],count,ruleSet,ruleleft,levelX,length+1,max_length,digits)
  ruleSet = thisList$ruleSet; count = thisList$count

  thisList = treeVisit(tree, tree[rowIx,"right daughter"],count,ruleSet,ruleright,levelX,length+1,max_length,digits)
  ruleSet = thisList$ruleSet; count = thisList$count

  res$ruleSet <- ruleSet
  res$count <- count
  return(res)
}

ruleList2Exec <-
  function(X,allRulesList){
    typeX = getTypeX(X)
    ruleExec <- unique(t(sapply(allRulesList,singleRuleList2Exec,typeX=typeX)))
    ruleExec <- t(ruleExec)
    colnames(ruleExec) <- "condition"
    return(ruleExec)
  }

singleRuleList2Exec <-
  function(ruleList,typeX){ #numeric: 1; categorical: 2s
    #ruleExec <- "which("
    ruleExec <- ""
    vars <- ls(ruleList)
    #ruleL <- length(unique(vars))
    vars <- vars[order(as.numeric(vars))]
    for(i in 1:length(vars)){
      if(typeX[as.numeric(vars[i])]==2){
        values <- paste("c(",paste(  paste("'",ruleList[[vars[i]]],"'",sep="")    ,collapse=","),")",sep="")
        tmp = paste("X[,",vars[i], "] %in% ", values, sep="")
      }else{
        tmp = ruleList[[vars[i]]]
      }
      if(i==1)ruleExec <- paste(ruleExec, tmp,sep="")
      if(i>1)ruleExec <- paste(ruleExec, " & ", tmp, sep="")
    }
    #ruleExec <- paste(ruleExec,")",sep="")
    return(c(ruleExec))
  }

getTypeX <-
  function(X){
    typeX = rep(0,ncol(X)-2)
    for(i in 1:ncol(X)-2){ #numeric: 1; categorical: 2s
      if(is.numeric(X[,i+2])){ typeX[i] = 1 }else{
        typeX[i] = 2
      }
    }
    return(typeX)
  }

# getTypeX <-
#   function(X){
#     typeX = rep(0,ncol(X))
#     for(i in 1:ncol(X)){ #numeric: 1; categorical: 2s
#       if(is.numeric(X[,i])){ typeX[i] = 1 }else{
#         typeX[i] = 2
#       }
#     }
#     return(typeX)
#   }

# getTree <- function(rfobj, k=1, labelVar=FALSE) {
#   if (is.null(rfobj$forest)) {
#     stop("No forest component in ", deparse(substitute(rfobj)))
#   }
#   if (k > rfobj$ntree) {
#     stop("There are fewer than ", k, "trees in the forest")
#   }
#   if (rfobj$type == "regression") {
#     tree <- cbind(rfobj$forest$leftDaughter[,k],
#                   rfobj$forest$rightDaughter[,k],
#                   rfobj$forest$bestvar[,k],
#                   rfobj$forest$xbestsplit[,k],
#                   rfobj$forest$nodestatus[,k],
#                   rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
#   } else {
#     tree <- cbind(rfobj$forest$treemap[,,k],
#                   rfobj$forest$bestvar[,k],
#                   rfobj$forest$xbestsplit[,k],
#                   rfobj$forest$nodestatus[,k],
#                   rfobj$forest$nodepred[,k])[1:rfobj$forest$ndbigtree[k],]
#   }
#
#   dimnames(tree) <- list(1:nrow(tree), c("left daughter", "right daughter",
#                                          "split var", "split point",
#                                          "status", "prediction"))
#
#   if (labelVar) {
#     tree <- as.data.frame(tree)
#     v <- tree[[3]]
#     v[v == 0] <- NA
#     tree[[3]] <- factor(rownames(rfobj$importance)[v])
#     if (rfobj$type == "classification") {
#       v <- tree[[6]]
#       v[! v %in% 1:nlevels(rfobj$y)] <- NA
#       tree[[6]] <- levels(rfobj$y)[v]
#     }
#   }
#   tree
# }


# extractRules.surv <-
#   function(treeList,X,ntree=500,maxdepth=6,random=FALSE,digits=NULL){
#     if(is.numeric(digits)) digits <- as.integer(abs(digits))
#
#     levelX = list()
#     for(iX in 2:ncol(X)){
#       levelX <- c(levelX,list(levels(X[,iX])))
#     }
#     # X <- NULL; target <- NULL
#     ntree=min(treeList$ntree,ntree)
#     allRulesList = list()
#     for(iTree in 1:ntree){
#       if(random==TRUE){max_length = sample(1:maxdepth,1,replace=FALSE)}else{
#         max_length = maxdepth}
#       rule = list(); count = 0; rowIx = 1;
#       # tree = getTree(rf,iTree,labelVar=FALSE)
#       tree <- treeList$list[[iTree]]
#       if(nrow(tree)<=1) next # skip if there is no split
#       ruleSet = vector("list", length(which(tree[,"status"]==-1)))
#       res <- treeVisit(tree,rowIx = rowIx,count,ruleSet,rule,levelX,length=0,max_length=max_length,digits=digits)
#       allRulesList = c(allRulesList, res$ruleSet)
#     }
#     allRulesList <- allRulesList[!unlist(lapply(allRulesList, is.null))]
#     cat(paste(length(allRulesList)," rules (length<=",
#               max_length, ") were extracted from the first ", ntree," trees.","\n",sep=""))
#
#     rulesExec <- ruleList2Exec.surv(X,allRulesList)
#     return(rulesExec)
#   }

