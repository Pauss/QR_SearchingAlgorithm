import matplotlib.pyplot as plt
import numpy as np


# use two global list in order to have easy access of them
report_match_all_columns = list()
report_match_some_columns = list()
report_RSS = list()
average_time = list()

n_list = 4

# function to get GA computed data in order to show it on graphics
def get_data(filepath):

    global report_match_all_columns
    global report_match_some_columns
    global report_RSS
    global average_time

    try:
        fd = open(filepath, 'r')
        line = fd.readline()
        cnt = 0
        while line:
            values = line.split(",")

            report_match_all_columns.append(float(values[0]))
            report_match_some_columns.append(float(values[1]))
            report_RSS.append(float(values[2]))
            average_time.append(float(values[3]))

            line = fd.readline()

            # save number of instances
            cnt += 1
    finally:
        fd.close()

    return cnt


# function to get GA computed data in order to show it on graphics
def get_data_tuning(filepath):

    global report_match_all_columns
    global report_match_some_columns
    global report_RSS
    global average_time

    try:
        fd = open(filepath, 'r')
        line = fd.readline()
        cnt = 0
        while line:
            values = line.split(",")

            report_match_all_columns.append(float(values[0]))
            report_match_some_columns.append(float(values[1]))
            report_RSS.append(float(values[2]))
            average_time.append(float(values[3]))

            line = fd.readline()

            # save number of instances
            cnt += 1
    finally:
        fd.close()

    return cnt


def get_graphic():

    ########################################
    # scatter the data

    data = list()

    for i in range(n_list):
        temp_list = []
        temp_list.extend([report_match_all_columns[i], report_match_some_columns[i], report_RSS[i], average_time[i]])
        data.append(temp_list)

    dim = len(data)

    w = 0.75
    dimw = w / dim

    fig, ax = plt.subplots()
    x = np.arange(len(data))
    x_names = ["GA", "HC", "SA", "GA_BB"]
    y_label = ["report match all columns", "report match some columns", "average AIC", "average time"]
    for i in range(len(data[0])):
        y = [d[i] for d in data]
        b = ax.bar(x + i * dimw, y, dimw, bottom=0.001, label=y_label[i])

    plt.setp(ax, xticks=x + dimw / 2, xticklabels=x_names)

    ax.legend(fancybox=True, framealpha=0.5, loc='upper right')
    ax.set_title('Comparing the performance of heuristic algorithms')

    fig.tight_layout()

    plt.show()
    plt.savefig('./output_graphics/My_Graphic_Reports.png')


def get_graphic_tuning(alg):

    ########################################
    # scatter the data
    fig, ax = plt.subplots()

    n_config = len(report_match_all_columns)
    x = np.arange(n_config)

    lines = plt.plot(x, report_match_all_columns, x, report_match_some_columns, x, report_RSS, x, average_time)
    plt.setp(lines[0], linewidth=3)
    plt.setp(lines[1], linewidth=4)
    plt.setp(lines[2], linewidth=2)
    plt.setp(lines[3], linewidth=1)

    plt.legend(("report match all columns", "report some all columns", "report AIC", "average time"), fancybox=True, framealpha=0.5, loc='upper right', fontsize = 'x-small')
    ax.set_title('Comparing the performance of ' + alg)

    # fig.tight_layout()

    plt.show()
    plt.savefig('./output_graphics/My_Graphic_Reports_Tuning.png')


def execute():
    n = get_data("./output_individuals/Out_Algorithms.csv")

    if n:
        get_graphic()
    else:
        print("Empty data, graphic picture could not be created")


def execute_tuning(alg):
    n = get_data_tuning("./output_individuals/Out_Algorithms_Tuning.csv")

    if n:
        get_graphic_tuning(alg)
    else:
        print("Empty data, graphic picture could not be created")