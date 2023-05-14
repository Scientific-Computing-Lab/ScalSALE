import os
import shutil
import argparse


description_help = """
*********************************************************************
-------- Deletes N latest checkpoints from all stores. --------

Usages:
1) remove_latest_checkpoints.py -n [Number of checkpoints]
...  

"""

default_n = 1
default_conf_path = '../CR/scr_conf.conf'

def delete_latest_checkpoints(n=default_n, conf_path=default_conf_path):
	# find all store locations for scr checkpoints
	with open(conf_path) as f:
		scr_conf = f.read()
	lines = scr_conf.split('\n')
	store_descriptions = filter(lambda x: x.startswith("STORE="), lines)
	stores = [store.split()[0][6:] for store in store_descriptions]
	stores = [store.replace('$USER', os.getenv('USER')) for store in stores]
	ckpts = []
	for store in stores:
		for root, dirs, files in os.walk(os.path.join(store, os.getenv("USER"))):
			for dir in dirs:
				if dir.startswith('scr.dataset.'):
					ckpts.append(os.path.join(root, dir))
	ckpts = [(int(ckpt.split('.')[-1]), ckpt) for ckpt in ckpts]
	# sort checkpoints by id
	ckpts.sort(lambda a, b: a[0] - b[0])
	# remove checkpoints from all locations except one
	last = max([ckpt[0] for ckpt in ckpts])
	to_delete = range(last-n+1, last+1)
	remaining = []
	for ckpt in ckpts:
		if ckpt[0] in to_delete:
			print "Deleting {0}".format(ckpt[1])
			shutil.rmtree(ckpt[1], ignore_errors=True)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=description_help,
		epilog=author_help)
	parser.add_argument("-n", action='store', dest='n', type=int, default=default_n,
	                    help='Number of latest checkpoints to remove from all stores.')
	parser.add_argument("-conf_path", action='store', dest='conf_path', default=default_conf_path,
	                    help='Path to scr_conf.conf.')
	args = parser.parse_args()
	delete_latest_checkpoints(n=args.n, conf_path=args.conf_path)

