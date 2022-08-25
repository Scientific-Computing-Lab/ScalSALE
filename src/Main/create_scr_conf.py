import argparse
import os

author_help = """
Author: Idan Mosseri (November, 2020) 
"""
description_help = """
*********************************************************************
-------- Creates the scr_conf file. --------

Usages:
1) create_scr_cof.py -n [run_name] -o [output scr_conf path] -t [scr_conf template path]
...  

"""
default_run_name = 'backus'
default_scr_conf_path = os.path.abspath(os.path.join(os.pardir, 'CR', 'scr_conf.conf'))
default_scr_conf_template_path = os.path.abspath(os.path.join(os.pardir, 'CR', 'scr_conf_template.conf'))
default_scr_prefix = os.path.abspath(os.path.join("gpfs", "fs_16m", "sc", os.getenv("USER"), "scr", "prefix"))
default_scr_runs = 3
default_checkpoint_seconds = 10
default_checkpoint_overhead = 0.0

def create_scr_conf(run_name=default_run_name,
                    scr_conf_path=default_scr_conf_path,
                    scr_conf_template_path=default_scr_conf_template_path,
					scr_prefix=default_scr_prefix,
                    scr_runs=default_scr_runs,
					checkpoint_seconds=default_checkpoint_seconds,
					checkpoint_overhead=default_checkpoint_overhead):
	scr_conf_path = os.path.abspath(scr_conf_path)
	scr_conf_template_path = os.path.abspath(scr_conf_template_path)
	print 'Creating scr_conf file in {0}'.format(scr_conf_path)
	with open(scr_conf_template_path) as f:
		scr_conf = f.read()

	scr_conf = scr_conf.replace('$RUN_NAME', run_name)
	scr_conf = scr_conf.replace('$SCR_PREFIX', os.path.join(scr_prefix, run_name, os.getenv("USER")))
	scr_conf = scr_conf.replace('$SCR_RUNS', str(scr_runs))
	scr_conf = scr_conf.replace('$SCR_CHECKPOINT_SECONDS', str(checkpoint_seconds))
	scr_conf = scr_conf.replace('$SCR_CHECKPOINT_OVERHEAD', str(checkpoint_overhead))
	with open(scr_conf_path, 'w') as f:
		f.write(scr_conf)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=description_help,
		epilog=author_help)
	parser.add_argument("-n", action='store', dest='run_name', default=default_run_name,
	                    help='A name describing the run. This will be used as the prefix to the scr datasets.')
	parser.add_argument("-o", action='store', dest='scr_conf_path', default=default_scr_conf_path,
	                    help='Path to output scr_conf file')
	parser.add_argument("-t", action='store', dest='scr_conf_template_path',
	                    help='Path to the scr_conf template to be used.')
	parser.add_argument("-p", action='store', dest='scr_prefix',
						help='Path to scr prefix.')
	parser.add_argument("-ckpt_seconds", action='store', dest='checkpoint_seconds', type=int,
						help='Set to positive number of seconds to specify minimum time between consecutive checkpoints as guided by SCR_Need_checkpoint.')
	parser.add_argument("-ckpt_overhead", action='store', dest='checkpoint_overhead', type=int,
						help='Set to positive percentage to specify maximum overhead allowed for checkpointing operations as guided by SCR_Need_checkpoint.')
	args = parser.parse_args()
	create_scr_conf(run_name=args.run_name,
					scr_conf_path=args.scr_conf_path,
					scr_conf_template_path=args.scr_conf_template_path,
					scr_prefix=args.scr_prefix,
					checkpoint_seconds=args.checkpoint_seconds,
					checkpoint_overhead=args.checkpoint_overhead)

