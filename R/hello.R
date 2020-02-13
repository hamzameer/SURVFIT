# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# devtools::load_all() | reload all your code
# shortcut: Ctrl-Shift-L | also saves all open files saving you a keystroke

# Ctrl + . | go to file/function from within a package
# Ctrl-F9 |go back to where you were
# source('script.R') | runs the whole script | never use to load code from a file | use devtools::load_all()
# R packages can only contain ASCII format in the .R files. To include unicode characters use the special escape format
# For example:
#
# x <- "This is a bullet <weird shit>"
# y <- "This is a bullet\u2022"
# identical(x,y) # FALSE
# cat(stringi::stri_escape_unicode(x)) #

# install.packages(c("devtools","roxygent2","testthat","knitr"))
# install.packages('rstudioapi')
# rstudioapi::isAvailable("0.99.149")
#
# # library(devtools)
# # has_devel() # great.
# # devtools::session_info()
#
# # use_package("Rcplex") #Rcplex ,osqp, "dplyr","ranger","survival","ggplot2","Formula")
# library(survival)
# data(ovarian)
# use_data(ovarian)

head(ovarian)
formula = Surv(futime, fustat)~.
# rules <- survfit(formula = formula, ovarian)
