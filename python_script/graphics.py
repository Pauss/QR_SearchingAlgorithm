
import matplotlib.pyplot as plt
import numpy as np

'''create input data for graphics'''

fitness = list()
nr_columns = list()

def get_data(filepath):

    global fitness
    global nr_columns
    try:
        file = open(filepath, 'r')
        line = file.readline()
        cnt = 0;
        while line:
            t = line.split(",")
            fitness.append(float(t[1]))
            nr_columns.append(int(t[2]))
            line = file.readline()
            cnt += 1
    finally:
        file.close()

    return cnt


def get_graphic():

    fig, ax = plt.subplots()
    ax.scatter(nr_columns[:], fitness[:],  c=nr_columns, s=fitness, alpha=0.5)

    ax.set_xlabel(r'$nr of columns$', fontsize=10)
    ax.set_ylabel(r'$fitness$', fontsize=10)
    ax.set_title('Fitness of candidates')

    ax.grid(True)
    fig.tight_layout()
    plt.savefig('../output_graphics/My_Graphic.png')

n = get_data("../output_individuals/Out_individuals.csv")
get_graphic()


'''
# Data for plotting
n = 10
t = list()

for i in range(n):
    t.append(i**2)


plt.figure()
# linear
# subplot: imparte figura in mai multe grafice: arg1 : nr de linii, arg2 nr de coloane, arg pe ce pozitie se afla graficul
#plt.subplot(221)
plt.plot(t, 'r')
plt.yscale('linear')
plt.title('linear2')
plt.grid(True)

plt.show()
'''

