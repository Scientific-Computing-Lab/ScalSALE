import argparse
import sys
import os
import math
import json
import time
import numpy as np
import default_value as dv
import datafile_fortran as dff
import replace_words as rw
import parser_3d as p3d
import parser_2d as p2d
from create_scr_conf import create_scr_conf
from mpi4py import MPI
import valid_datafile as checker
import subprocess


description_help = """
*********************************************************************
-------- Parses Backus's datafile and runs the fortran code. --------

Usages:
1) for a specific datafile: main.py -df
...  

"""
RED   = "\033[1;31m"
BLUE  = "\033[1;34m"
CYAN  = "\033[1;36m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
BOLD    = "\033[;1m"
REVERSE = "\033[;7m"


def set_datafile(f_name):
    with open(f_name) as f:
        j_object = json.load(f)

    j_object = rw.replace_words(j_object)
    dv.set_default_2d(j_object)

    create_scr_conf(run_name=j_object['checkpoint_restart']['run_name'],
                    scr_prefix=j_object['checkpoint_restart']['scr_prefix'],
                    checkpoint_seconds=j_object['checkpoint_restart']['checkpoint_seconds'],
                    checkpoint_overhead=j_object['checkpoint_restart']['checkpoint_overhead'])

    p2d.set_datafile_2d(j_object, dff)
    if j_object["data"]["name"].__contains__("3d"):
        p3d.set_datafile_3d(j_object, dff)
        dff.f2py_module.dimension = 3
    else:
        dff.f2py_module.dimension = 2
    f.close()
    # We check the datafile after we parse it, for convenience
    err, msg = checker.check_datafile(j_object, dff.f2py_module)
    manage_checker(err, msg)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    dff.f2py_module.comm = comm.py2f()
    print ("Done parsing")
    dff.f2py_module.send_input()


def manage_checker(err, msg):
    if err == 0:
        sys.stdout.write(RED)
        print ("Error in datafile: ", msg)
        print ("... Exiting")
        sys.stdout.write(RESET)
        exit(1)
    if err == -1: #warning
        sys.stdout.write(RED)
        print ("Warning in datafile: ", msg)
        sys.stdout.write(RESET)
        print ("... Continuing")


def convert_array_XXXX(arr):
    return np.array(arr, dtype=np.double)


def convert_array_string(arr):
    return np.array(arr, dtype=np.chararray)


# Returns the reduct number of materials from the datafile
def reduction_number_material(ind_arr):
    witness_arr = [None] * (len(ind_arr) + 1)
    for i in range(len(witness_arr)):
        witness_arr[i] = -1
    number = 0
    seen = False
    for i in range(len(ind_arr)):
        seen = False
        for j in range(len(witness_arr)):
            # we reached end of the array and we didn't see the index
            if witness_arr[j] == -1 and not seen:
                witness_arr[j] = ind_arr[i]
                break
            if witness_arr[j] == ind_arr[i]:
                seen = True
                break

        if not seen and ind_arr[i] != 0:
            number = number + 1
    return number


# TODO : change this 
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description_help
                                     , epilog=author_help)
    parser.add_argument("-df", action='store', dest='df_json', help="Datafile for Backus")
    options = parser.parse_args()
    df_json = options.df_json
    if not options.df_json:
        path = os.path.curdir
        if os.path.abspath(path).__contains__("Scripts"):
            df_json = "../Datafiles/datafile.json"
        elif os.path.abspath(path).__contains__("src"):
            df_json = "Datafiles/datafile.json"
        else:
            df_json = "src/Scripts/datafile.json"

    set_datafile(df_json)


if __name__ == '__main__':
    main()
    #   dff.datafile_interface_module.nxb = 10
    #dff.datafile_interface_module.nyb = 10
    #dff.datafile_interface_module.send_input()
    #print dff.datafile_interface_module.nxb
