#the purpose of this script is to be an automated controller script that runs runn_nnsearch_hpc.py on the entire enamine library for a given ligand shape
#the user needs to at least give an aligned ligand molecule shape (mol2 or sdf), and can optionally give a path to a location to put the results in

#imports 
import os,sys

#get the ligand file
target_molecule_file = sys.argv[1]

#create the result output in a specified location, or use location where this was called from otherwise
working_location = "."
if len(sys.argv) > 2:
	working_location = sys.argv[2]

#iterate over every chunk, and call the search script on each subchunk within the chunk
#for i in range(0,53085):
for i in range(0,1500):
	#build the chunk string
	chunk_str = str(i)

	#append leading zeroes until the string is 5 digits long
	while len(chunk_str) < 5:
		chunk_str = "0" + chunk_str

	#pause 0.1 seconds to help avoid overload
	os.system("sleep 0.1")

	#submit a bsub job that runs the nnsearch python script
	print("bsub -q long -W 8:00 -u \"\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/shapedb/run_nnsearch_hpc.py " + chunk_str + " " +  target_molecule_file + " " + working_location + " \"")
	os.system("bsub -q long -W 8:00 -u \"\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/shapedb/run_nnsearch_hpc.py " + chunk_str + " " +  target_molecule_file + " " + working_location + " \"")