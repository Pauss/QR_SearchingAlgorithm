
import matplotlib.pyplot as plt
import numpy as np


# use two global list in order to have easy access of them
iterations = list()
fitness = list()
nr_columns = list()
columns = list()


# function to get GA computed data in order to show it on graphics
def get_data(filepath):

    global fitness
    global nr_columns
    global iterations
    global columns

    try:
        fd = open(filepath, 'r')
        line = fd.readline()
        cnt = 0
        while line:
            values = line.split(",")
            # get size of an individual in order to get all columns that is composed of
            size_of_individual = int(values[2])

            # for al columns of an individual set same RSS, iteration and size
            for i in range(size_of_individual):
                iterations.append(int(values[0]))
                fitness.append(float(values[1]))
                nr_columns.append(int(values[2]))
                columns.append(int(values[3+i]))
            line = fd.readline()

            # save number of instances
            cnt += 1
    finally:
        fd.close()

    return cnt


def get_graphic():

    # used for setting up style adjustments
    font = {'family': 'sans-serif',
            'color': 'black',
            'weight': 'normal',
            'size': 8,
            'ha': 'center',
            }

    # create a figure
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=False, sharey=True)

    ########################################
    # first graphic
    # scatter the data
    index_min_fitness = fitness.index(min(fitness))
    index_min = iterations[index_min_fitness]

    x1fit = [min(iterations), index_min]
    y1fit = [fitness[0], min(fitness)]
    ax1.scatter(iterations, fitness, c=iterations, s=fitness, marker=".", alpha=0.4)
    ax1.plot(x1fit, y1fit, "r-")

    # customization
    ax1.set_xlabel('Iteration', fontdict=font)
    ax1.set_ylabel('AIC', fontdict=font)
    # ax1.set_title('Visualisation of fitness and dimension of individuals', fontdict=font)
    ax1.grid(True)
    ########################################

    # second graphic
    # scatter the data
    x2fit = nr_columns[index_min_fitness]
    y2fit = min(fitness)
    ax2.scatter(nr_columns, fitness, c=nr_columns, s=fitness, marker=".", alpha=0.5)
    ax2.plot(x2fit, y2fit, "r*")

    # customization
    ax2.set_xlabel('Number of columns', fontdict=font)
    ax2.set_ylabel('AIC', fontdict=font)
    # ax2.set_title('Visualisation of progress of the fitness within iterations', fontdict=font)
    ax2.grid(True)
    ########################################

    # third graphic
    # scatter the data
    x3fit = list()
    y3fit = list()
    for i in range(index_min_fitness, (index_min_fitness + nr_columns[index_min_fitness])):
        x3fit.append(columns[i])
        y3fit.append(min(fitness))

    ax3.scatter(columns, fitness, c=columns, s=fitness, marker=".", alpha=0.4)
    ax3.plot(x3fit, y3fit, "r*")

    # customization
    ax3.set_xlabel('Individual columns', fontdict=font)
    ax3.set_ylabel('AIC', fontdict=font)
    # ax3.set_title('Visualisation of progress of the fitness within iterations', fontdict=font)
    ax3.grid(True)
    ########################################

    fig.tight_layout()
    # Add a legend
    # ax1.legend()

    plt.show()
    plt.savefig('../output_graphics/My_Graphic.png')


n = get_data("../output_individuals/Out_individuals.csv")

if n:
    get_graphic()
else:
    print("Empty data, graphic picture could not be created")
