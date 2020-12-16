SURVFIT example
================

#### Applying SURVFIT to the *ovarian* data-set from the **survival** package. A minimal working example.

``` r
# load all the dependencies
library(survival)
library(ranger)
library(Rcplex)
```

    ## Loading required package: slam

``` r
library(SURVFIT)
```

``` r
# load the ovarian data-set
data(ovarian)
# run the model
model <- SURVFIT(Surv(futime, fustat)~., data = ovarian)
```

    ## 400 rules (length<=1) were extracted from the first 200 trees.
    ## 710 rules (length<=2) were extracted from the first 200 trees.
    ## 987 rules (length<=3) were extracted from the first 200 trees.
    ## 1/25: CV for gamma=1e-06, lamba = 10
    ## 2/25: CV for gamma=0, lamba = 2
    ## 3/25: CV for gamma=0.001, lamba = 1
    ## 4/25: CV for gamma=2, lamba = 0.1
    ## 5/25: CV for gamma=0.1, lamba = 2
    ## 6/25: CV for gamma=5e-06, lamba = 20
    ## 7/25: CV for gamma=1, lamba = 1e-08
    ## 8/25: CV for gamma=1e-08, lamba = 2
    ## 9/25: CV for gamma=0.005, lamba = 0.01
    ## 10/25: CV for gamma=1e-05, lamba = 2
    ## 11/25: CV for gamma=0.001, lamba = 0.001
    ## 12/25: CV for gamma=1e-05, lamba = 1e-06
    ## 13/25: CV for gamma=5e-05, lamba = 10
    ## 14/25: CV for gamma=0.01, lamba = 1e-04
    ## 15/25: CV for gamma=0, lamba = 0.1
    ## 16/25: CV for gamma=5e-06, lamba = 0.001
    ## 17/25: CV for gamma=1, lamba = 1
    ## 18/25: CV for gamma=1, lamba = 0.1
    ## 19/25: CV for gamma=0.005, lamba = 1
    ## 20/25: CV for gamma=0.1, lamba = 1e-04
    ## 21/25: CV for gamma=1e-06, lamba = 0.01
    ## 22/25: CV for gamma=5e-04, lamba = 1e-08
    ## 23/25: CV for gamma=1e-08, lamba = 10
    ## 24/25: CV for gamma=1e-08, lamba = 1e-04
    ## 25/25: CV for gamma=5e-04, lamba = 1e-06
    ## Rcplex: num variables=1766 num constraints=14

#### To get the rules extracted from SURVFIT:

``` r
# rules extracted from SURVFIT are:
survfit_rules <- model$rules
print(survfit_rules)
```

    ##  [1] "X[,1]<=64.3 & X[,1]<=56.7479 & X[,2]<=1.5"         
    ##  [2] "X[,1]<=58.263 & X[,3]>1.5 & X[,4]>1.5"             
    ##  [3] "X[,1]<=64.3 & X[,1]>56.7479 & X[,2]<=1.5"          
    ##  [4] "X[,1]<=61.537 & X[,1]>41.2041 & X[,3]<=1.5"        
    ##  [5] "X[,1]<=55.8041 & X[,2]<=1.5 & X[,4]>1.5"           
    ##  [6] "X[,1]<=67.77535 & X[,1]<=59.74245 & X[,1]>58.263"  
    ##  [7] "X[,1]>59.74245 & X[,3]>1.5"                        
    ##  [8] "X[,1]<=58.9493 & X[,2]<=1.5 & X[,3]>1.5"           
    ##  [9] "X[,1]<=59.7219 & X[,1]>56.99455 & X[,3]>1.5"       
    ## [10] "X[,1]<=64.3 & X[,1]>47.46985 & X[,1]<=56.1151"     
    ## [11] "X[,1]<=51.85205 & X[,1]<=43.8685 & X[,4]<=1.5"     
    ## [12] "X[,1]>58.9493 & X[,2]<=1.5 & X[,3]>1.5"            
    ## [13] "X[,1]<=65.32055 & X[,1]>54.54245 & X[,1]<=57.6233" 
    ## [14] "X[,1]<=54.38765 & X[,3]<=1.5 & X[,4]<=1.5"         
    ## [15] "X[,1]<=61.42465 & X[,1]<=56.99455 & X[,1]>50.22465"
    ## [16] "X[,1]>59.74245 & X[,4]>1.5"

#### To get all rules extracted and their coefficients in the model:

``` r
# rules extracted from SURVFIT are:
all_rules <- model$all_rules
beta <- model$beta

# printing rules and coefficients of first 10 rules
print(all_rules[1:10])
```

    ##  [1] "X[,1]<=68.3781"  "X[,1]>68.3781"   "X[,1]<=65.44525" "X[,1]>65.44525" 
    ##  [5] "X[,3]<=1.5"      "X[,3]>1.5"       "X[,4]<=1.5"      "X[,4]>1.5"      
    ##  [9] "X[,1]<=69.33425" "X[,1]>69.33425"

``` r
print(beta[1:10])
```

    ##  [1]  1.107661  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000
    ##  [8] -9.641074  1.063347  0.000000
