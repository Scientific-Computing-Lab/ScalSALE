import argparse
import os
import json
import imp
import numpy
import subprocess #subprocess.call(['./module_load.sh'])

description_help = """
*********************************************************************
-------- Parses Backus's datafile and runs the fortran code. --------

Usages:
1) for a specific datafile: main.py -df
...  

"""

author_help = "Parsing done!!!"

RED   = "\033[1;31m"
BLUE  = "\033[1;34m"
CYAN  = "\033[1;36m"
GREEN = "\033[0;32m"
RESET = "\033[0;0m"
BOLD    = "\033[;1m"
REVERSE = "\033[;7m"


def runner(df, gnu):
    f = open(df)
    j_object = json.load(f)
    np = get_np(j_object)
    threads = get_threads(j_object)
    f.close()

    sbatch_file_path = 'run_slurm.sh'
    if os.path.isfile(sbatch_file_path):
        os.remove(sbatch_file_path)
    sbatch_file = open(sbatch_file_path, "w")
    os.chmod(sbatch_file_path, 0o777)
    N = int(numpy.ceil(np / 32.0))

    #execute_command = "mpirun -n {0} -mca btl self,sm,openib ../exec/main".format(np)
    execute_command = "mpirun -n {0} ../exec/main".format(np)
    if (gnu):
        module_load = 'module load gnu/9.1.0 cmake openmpi/4.1.3-gnu \n'
    else:
        module_load = 'module load intel/18.0.1.163 openmpi/4.1.3-intel\n'
   # job_submit_str = '#!/bin/bash\n#SBATCH -n {0} -N {1} --exclusive --threads-per-core={2} -p mixedp --error=slurm-%j.err --output=slurm-%j.out\nexport OMP_NUM__THREADS={3}\nmodule load intel/18.0.1.163 openmpi/4.0.4_intel mpi/impi-intel2018 cmake anaconda2\n{4}\n'.format(np, N, threads,threads, execute_command)
    
    #print(job_submit_str)
    
    #sbatch_file.write(job_submit_str)
    sbatch_file.write(
        '#!/bin/bash\n'
##        + '#SBATCH -n {0} -N {1} --exclusive --threads-per-core={2} -p mixedp -x node0[01,06-25]  --error=slurm-%j.err --output=slurm-%j.out\n'.format(
        + '#SBATCH -n {0} -N {1} --exclusive --threads-per-core={2} -p mixedp  --error=slurm-%j.err --output=slurm-%j.out\n'.format(
            np, N, threads)
        + 'export OMP_NUM_THREADS={0}\n ml purge\n'.format(threads)
        #+ 'module load intel/18.0.1.163 openmpi/4.0.4_intel mpi/impi-intel2018 cmake anaconda2 ScientificLibraries/silo/4.11\n'
        + module_load
        + '{0}\n'.format(execute_command)
    )    

    sbatch_file.close()
    #
    sub_proc = subprocess.Popen(['sbatch', sbatch_file_path], cwd=os.path.curdir)
    sub_proc.wait()


def get_np(j_obj):
    if not j_obj.has_key("parallel"):
        return 1
    else:
        return j_obj["parallel"]["np"]


def get_threads(j_obj):
    if not j_obj.has_key("parallel"):
        return 1
    else:
        if j_obj.has_key("threads"):
            return j_obj["parallel"]["threads"]
        else:
            return 1
	  


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=description_help
                                     , epilog=author_help)
    parser.add_argument("-gnu", action='store_true', dest='gnu', help="Run GNU", default=False)

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
    runner(os.path.abspath(df_json), options.gnu)

