

## Generate a dataset that contains a true submodel of 'ntrue' regressors

## The data matrix X is drawn from the standard normal distribution

## The coefficients are 0 or 1. There a 'ntrue' 1 coefficients randomly
##     selected out of the range 1..nreg

## The error vector is drawn from the normal distribution with zero mean and
## standard deviation 'sd'


## nobs: numeric, number of observation in the model
## nreg: numeric, total number of regressors
## ntrue: numeric, number of "true" regressors
## intercept: logical
## sd: numeric, standard deviation of the error term

M_PI = 3.14159265358979323846

model <- function(nobs = 80, nreg = 60, ntrue = 30 , intercept = TRUE, sd = 0.01)
{
    ## independent variables
    X <- rnorm(nobs * nreg,sd=0.01)

    X = to_6_digits(X)

    dim(X) <- c(nobs, nreg)

    
    ## coefficients
    coefs <- rep_len(0, length.out = nreg)

    ##coefs_real <- runif(nreg, min=-0.5, max=1.5)

    coefs = to_6_digits(coefs)

    true <- sample.int(nreg, size = ntrue)

    coefs[true] <- 1 ##coefs_real[true]

    ## error
    error <- rnorm(nobs, sd = sd)

    error = to_6_digits(error)


    ## dependent variable
    y <- cbind(intercept, X) %*% c(1, coefs) + error

    y = to_6_digits(y)

    model = list(y = y, X = X, nobs = nobs, nreg = nreg, ntrue= ntrue, intercept = intercept, true = true, error = error)
 
}

norm_vec <- function(x) (sum(x^2))
fitness_func <- function(n, k, RSS) ((n + n * log(2 * M_PI) + n * log(RSS / n)+ 2 * (k + 2))) #k+1 without intercept

to_6_digits <- function(x)
{
    for (i in seq_along(x)){
     x[i] = as.numeric(formatC(round(x[i], 6), digits = 6, format = "f"))
        ## x[i] = as.numeric(format(round(x[i], 5), nsmall = 5))
    }

    to_6_digits = x
}

inc <- function(x)
{
 eval.parent(substitute(x <- x + 1))
}

write_to_file <-function(filename)
{
    temp_model = model()

    setwd("./generated_data")
    file.create(filename)

    ##get number of lines and columns to be written to file
    nlines = as.numeric(temp_model[3])
    ncolumns = as.numeric(temp_model[4]) + 1
    ncolumns_model = as.numeric(temp_model[5])

    y = unlist(temp_model[1])
    x = unlist(temp_model[2])

    combine = matrix(c(y,x),nrow = nlines,ncol = ncolumns,byrow = FALSE)

    write.table(combine, file=filename, row.names=FALSE, col.names=TRUE)

    best = (unlist(temp_model[7]))

    #print(temp_model[7])
    print(best)

    e = unlist(temp_model[8])
    RSS = norm_vec(e)
    AIC = fitness_func(ncolumns, ncolumns_model, RSS)

    print(AIC)
}

filename <- commandArgs(trailingOnly = TRUE)
##filename = "Gdata0.txt"
write_to_file(filename)
