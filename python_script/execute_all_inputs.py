import subprocess
import sys
import time
#########################################################
# Execute QR_SearchingAlgorithm on generated files
#########################################################
ga_dict = {
    # genetic algorithm
    'ga1': ["ga", "1", "1", "1"],
    'ga2': ["ga", "1", "2", "1"],
    'ga3': ["ga", "1", "3", "1"],
    'ga4': ["ga", "1", "4", "1"],
    'ga5': ["ga", "2", "2", "1"],
    'ga6': ["ga", "2", "3", "1"],
    'ga7': ["ga", "2", "4", "1"],
    'ga8': ["ga", "4", "3", "1"],
    'ga9': ["ga", "4", "4", "1"],
    'ga10': ["ga", "1", "1", "2"],
    'ga11': ["ga", "1", "2", "2"],
    'ga12': ["ga", "1", "3", "2"],
    'ga13': ["ga", "1", "4", "2"],
    'ga14': ["ga", "2", "2", "2"],
    'ga15': ["ga", "2", "3", "2"],
    'ga16': ["ga", "2", "4", "2"],
    'ga17': ["ga", "4", "3", "2"],
    'ga18': ["ga", "4", "4", "2"],
}
ga_sa_dict = {
    # simulated annealing
    'ga_sa1': ["ga_sa", "1", "1", "2"],
    'ga_sa2': ["ga_sa", "1", "2", "2"],
    'ga_sa3': ["ga_sa", "1", "3", "2"],
    'ga_sa4': ["ga_sa", "1", "4", "2"],
    'ga_sa5': ["ga_sa", "2", "2", "2"],
    'ga_sa6': ["ga_sa", "2", "3", "2"],
    'ga_sa7': ["ga_sa", "2", "4", "2"],
    'ga_sa8': ["ga_sa", "4", "3", "2"],
    'ga_sa9': ["ga_sa", "4", "4", "2"],
}
ga_hc_dict = {
    # hill climbing
    'ga_hc1': ["ga_hc", "1", "1", "2"],
    'ga_hc2': ["ga_hc", "1", "2", "2"],
    'ga_hc3': ["ga_hc", "1", "3", "2"],
    'ga_hc4': ["ga_hc", "1", "4", "2"],
    'ga_hc5': ["ga_hc", "2", "2", "2"],
    'ga_hc6': ["ga_hc", "2", "3", "2"],
    'ga_hc7': ["ga_hc", "2", "4", "2"],
    'ga_hc8': ["ga_hc", "4", "3", "2"],
    'ga_gc9': ["ga_hc", "4", "4", "2"],
}
ga_bb_dict = {
    # Building Blocks
    'ga_bb1': ["ga_bb", "1", "1", "2"],
    'ga_bb2': ["ga_bb", "1", "2", "2"],
    'ga_bb3': ["ga_bb", "1", "3", "2"],
    'ga_bb4': ["ga_bb", "1", "4", "2"],
    'ga_bb5': ["ga_bb", "2", "2", "2"],
    'ga_bb6': ["ga_bb", "2", "3", "2"],
    'ga_bb7': ["ga_bb", "2", "4", "2"],
    'ga_bb8': ["ga_bb", "4", "3", "2"],
    'ga_bb9': ["ga_bb", "4", "4", "2"],
}
#########################################################
# Global variables
r_list = list()
# max number of columns
n_max = 10
# number of files to be generated and verified
n = 5
# details of files to be first generated and then checked
name = "GData"
type_file = ".txt"
#########################################################


def list_prepare(out):

    temp_list = list()

    temp = out.split(" ")
    temp = [maybe_int(v) for v in temp]
    temp = [v for v in temp if v is not None]

    temp.sort()
    temp_list.append(len(temp))
    temp_list.append(temp)

    return temp_list


def maybe_int(s):
    try:
        return int(s)
    except (ValueError, TypeError):
        return None


def execute_r():
    global r_list

    # number of files to be generated and verified
    global n
    # path to output files
    global name
    global type_file

    path_script = "./r_generator/generate_data.R"
    command = 'Rscript'

    # execute R Script
    try:
        for i in range(n):
            # Define arguments
            file_path = [name + str(i)+type_file]

            # Build subprocess command
            cmd = [command, path_script] + file_path
            x = subprocess.check_output(cmd).splitlines()

            # Get output
            out = x[1].decode(sys.stdout.encoding)

            new_l = list_prepare(out)
            r_list.append(new_l)

    except FileNotFoundError as err1:
        print("In Execute_r: ", err1)
    except IndexError as err2:
        print("\nIn Execute_r: ", err2)


def execute_c(c_list, par):

    # number of files to be generated and verified
    global n
    # path to output files
    global name
    global type_file

    command = "./Debug/QR_SearchingAlgorithm.exe"

    # execute c script
    try:
        for i in range(n):

            # Define arguments
            file_path = [name + str(i) + type_file]
            args = par

            # Build subprocess command
            cmd = [command] + args + file_path
            x = subprocess.check_output(cmd, cwd='./Debug').splitlines()

            # Get output
            out = x[5].decode(sys.stdout.encoding)

            new_l = list_prepare(out)
            c_list.append(new_l)

    except FileNotFoundError as err1:
        print("In Execute_c: ", err1)
    except IndexError as err2:
        print("In Execute_c: ", err2)
    except subprocess.TimeoutExpired as err3:
        print("In Execute_c: ", err3.timeout)


def matching_all_columns(list1, list2):

    # this function find all matching solutions
    # a solution is matching if all columns are the same as columns of best solution
    n_obsv = len(list1)
    n_match = 0

    for i in range(n_obsv):
        if list1[i][0] == list2[i][0]:
            temp_n = [z for z, j in zip(list1[i][1], list2[i][1]) if z == j]
            if len(temp_n) == list1[i][0]:
                n_match += 1
    probability = (float(n_match / n_obsv)) * 100
    return probability


def matching_some_columns(list1, list2):

    # this function find all matching columns of an individual
    # the probability is computed obtaining all columns that are the same as columns of best solution
    n_obsv = len(list1)
    probability = 0

    for i in range(n_obsv):
        temp_n = 0
        for j in list1[i][1]:
            if j in list2[i][1]:
                temp_n += 1

        temp_p = float(temp_n) / list1[0][0]
        probability += temp_p

    probability = (float(probability / n_obsv)) * 100
    return probability


def compare_outputs():
    global r_list

    execute_r()

    print("==========================================")

    # combine all dictionaries into one
    dicts = list()

    dicts.append(ga_dict)
    dicts.append(ga_hc_dict)
    dicts.append(ga_sa_dict)
    dicts.append(ga_bb_dict)

    super_dict = {}
    for d in dicts:
        for k, v in d.items():
            super_dict[k] = v

    # compare results
    try:
        for k in super_dict.keys():
            c_list = list()
            print(super_dict[k])
            execute_c(c_list, super_dict[k])

            print(r_list)
            print(c_list)

            p1 = matching_all_columns(r_list, c_list)
            p2 = matching_some_columns(r_list, c_list)

            print("Matching probability by all columns: ", p1)
            print("Matching probability by some columns: ", p2)
            print("==========================================\n")
            time.sleep(0.2)

    except Exception as err:
        print(err)


compare_outputs()

