import argparse
import os
import json
import imp
import subprocess #subprocess.call(['./module_load.sh'])

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


def runner(df):
#    subprocess.call(['source ../Scripts/module_load.sh'])
    aa = subprocess.Popen(['source ../Scripts/module_load.sh'] ,shell=True, stdout=subprocess.PIPE)
    aa.wait()
    
    f = open(df)
    j_object = json.load(f)
    np = get_np(j_object)
    f.close()
    a = os.system('mpirun -n {0} --mca btl self,sm,openib ../exec/main '.format(np, df))


def get_np(j_obj):
    if not j_obj.has_key("parallel"):
        return 1
    else:
        return j_obj["parallel"]["np"]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description_help
                                     , epilog=author_help)
    parser.add_argument("-df", action='store', dest='df_json', help="Datafile for Backus")
    options = parser.parse_args()
    df_json = "../Datafiles/datafile.json"
    if not options.df_json:
        path = os.path.curdir
        if os.path.abspath(path).__contains__("Scripts"):
            df_json = "../Datafiles/datafile.json"
        elif os.path.abspath(path).__contains__("Main"):
                df_json = "../Datafiles/datafile.json"
        elif os.path.abspath(path).__contains__("src") :
            df_json = "Datafiles/datafile.json"
        else:
            df_json = "src/Scripts/datafile.json"
    runner(os.path.abspath(df_json))

