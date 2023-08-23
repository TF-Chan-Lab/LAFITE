# AUTOGENERATED! DO NOT EDIT! File to edit: ../00_utils.ipynb.

# %% auto 0
__all__ = ['temp_dir_creation', 'keep_tmp_file', 'bam2bed', 'cmd_execution', 'loc_distance', 'Vividict']

# %% ../00_utils.ipynb 2
import math, tempfile, os, sys, subprocess
from bisect import bisect_left

# %% ../00_utils.ipynb 3
def temp_dir_creation(output_dir):
	"""
	TemporaryDirectory creation under output directory
	"""
	tmp_folder = tempfile.TemporaryDirectory(dir=output_dir)

	return tmp_folder

# %% ../00_utils.ipynb 4
def keep_tmp_file(output, tmp_dir):
	"""
	save intermediate files
	"""

	inter_dir = f'{os.path.splitext(output)[0]}_tmp'
	return f'''
	mkdir -p {inter_dir} && cp {tmp_dir}/* {inter_dir}
	'''

# %% ../00_utils.ipynb 5
def bam2bed(bam, output_dir, bedtools):
	"""
	return a command line that covert the bam file to a bed12 format file
	"""
	out = os.path.join(output_dir, 'bam.bed')
	bam2bed_cmd = f"{bedtools} bamtobed -bed12 -cigar -i {bam}|awk '$5>=1' > {out}"
	
	return bam2bed_cmd

# %% ../00_utils.ipynb 6
def cmd_execution(command):
	"""
	shell command execution
	"""
	p = subprocess.run(command,shell=True)
	if p.returncode == 0:
		sys.stdout.write(command +'\n')
	else:
		sys.stdout.write(command +'\n')

# %% ../00_utils.ipynb 7
def loc_distance(loc_list, loc):
	"""
	return the minimum relative distance between the splicing site in reference annotation and the query one
	"""
	loc_list = list(loc_list)
	if loc_list:
		pos = bisect_left(loc_list, loc)
		if pos == 0:
			loc_dis = abs(loc_list[0]-loc)
			ref_loc = loc_list[0]
		elif pos == len(loc_list):
			loc_dis = abs(loc - loc_list[-1])
			ref_loc = loc_list[-1]
		elif loc - loc_list[pos-1] >= loc_list[pos] - loc:
			loc_dis = loc_list[pos] - loc
			ref_loc = loc_list[pos]
		else:
			loc_dis = loc - loc_list[pos-1]
			ref_loc = loc_list[pos-1]
	else:
		loc_dis = ref_loc = math.inf
	
	return loc_dis, ref_loc

# %% ../00_utils.ipynb 8
class Vividict(dict):
	"""create nest dictionary
	"""
	def __missing__(self, key):
		value = self[key] = type(self)()
		return value
