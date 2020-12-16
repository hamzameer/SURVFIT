# # install SURVFIT and RCplex
#
# # install.packages('SURVFIT')
#
#
# # load required packages
# library(SURVFIT)
# library(Rcplex)
# library(ranger)
#
#
# # load example data
# library(survival)
# data(ovarian)
# head(ovarian)
# dim(ovarian)
#
#
# # remember to type SURVFIT in caps
#
# sfit <- SURVFIT(Surv(futime, fustat) ~., data = ovarian, gamma = 0.5, lambda1 = 2,
#                 lambda2 = 2, crossvalidate = FALSE )
