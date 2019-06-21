import numpy as np
import sys

## nobs: numeric, number of observation in the model
## nreg: numeric, total number of regressors
## ntrue: numeric, number of "true" regressors
## intercept: logical
## sd: numeric, standard deviation of the error term

nobs = 100
nreg = 30
ntrue = 15
sd = 0.01
intercept = True
M_PI = 3.14159265358979323846


def generate_data():

    # generate X
    X = np.random.normal(size=nobs*nreg, scale=sd)
    X = np.reshape(X, (nobs, nreg))

    # generate coefficients
    coeffs = np.zeros(nreg)
    coeffs[:ntrue] = 1
    np.random.shuffle(coeffs)

    true = [i+1 for i, x in enumerate(coeffs) if x == 1]

    # generate error
    error = np.random.normal(size=nobs, scale=sd)

    # add intercept
    intercept = np.ones(shape=(nobs,1), dtype=np.int8)

    temp_x =      np.concatenate((intercept, X), axis=1)
    temp_coeffs = np.concatenate((np.array([1]), coeffs), axis=0)
    y = np.matmul(temp_x, temp_coeffs) + error

    return [y, X, error, nobs, nreg, ntrue, true]


def write_to_file(filename):

    temp_model = generate_data()

    y = temp_model[0].copy()
    x = temp_model[1].copy()
    err = temp_model[2].copy()
    best = temp_model[6].copy()
    n = temp_model[4]
    k = temp_model[5]

    filepath = "./generated_data/"

    file = filepath + filename

    try:
        f = open(file, "w+")

        # write first line
        for i in range(temp_model[4]+1):
            f.write('column{} '.format(i+1))
        f.write("\n")

        # write rest of lines
        for i in range(temp_model[3]):
            f.write('{} '.format(y[i]))
            for j in range(temp_model[4]):
                f.write('{} '.format(x[i][j]))
            f.write("\n")

    except Exception as err:
        print("Write to file error: ",err)

    try:
        RSS = sum([x**2 for x in err])
        AIC = (n + n * np.log(2 * M_PI) + n * np.log(RSS / n) + 2 * (k + 2)) #k+1 without intercept

    except ZeroDivisionError as err2:
        print("Compute RSS/AIC error", err2)

    if len(best) > 0:
        return [best, AIC]


generate_data()
