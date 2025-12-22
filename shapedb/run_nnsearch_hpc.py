#the purpose of this script is to run the shapedb NNSearch out of the shapedb container on a given chunk of the enamine library
#this will write shape similarity scores to output text files for processing to select te most similar ligands
#note, the target shape must be pre-aligned against the rest of the enamine library (ideally, aligned against suvo_shifted.sdf) or else you are almost guaranteed to get all 1 values, which is useless

#ideally, this script can be called by nnsearch_controller.py to automate running the script on the entire enamine library

#this will operate on a per-chunk basis, hitting all (up to) 10 subchunks within the chunk, so this can be compartmentalized into one chunk per job

#importantly, this script will also be a preliminary filter to remove any blacklisted ligands that are removed for having unwanted chemical fingerprints (or conformer shape, if we ever do that)

#imports
import os,sys

#get the chunk as a 5 digit code (i.e. the first chunk is 00000, not 0)
working_chunk = str(sys.argv[1])

#derive the superchunk that this chunk belongs in for safer result storage (so we don't explode a directory with 53k directories)
superchunk_str = working_chunk[0:3]

#superchunk does not have preceeding zeroes, so cut any off, doing via casting
superchunk_str = int(superchunk_str)
superchunk_str = str(superchunk_str)

#get the target shape molecule (Can have a path to it, must be mol2 or sdf, or shapedb likely won't be able to read it)
target_molecule_file = sys.argv[2]

#create the result output in a specified location, or use location where this was called from otherwise
working_location = "."
if len(sys.argv) > 3:
	working_location = sys.argv[3]

#move to the working location
os.chdir(working_location)

#make a directory to perform the analysis in, based on the superchunk name
os.system("mkdir -p " + superchunk_str)
os.chdir(superchunk_str)
os.system("mkdir -p " + working_chunk)
os.chdir(working_chunk)

#iteratively decompress the params and db data to the working chunk location
for i in range (0,10):
	os.system("tar -xzf /pi/summer.thyme-umw/enamine-REAL-2.6billion/" + superchunk_str + "/" + working_chunk + "/condensed_params_and_db_" + str(i) + ".tar.gz -C .")

	#run shapedb out of the container on the subchunk database, executed via singularity
	os.system("singularity exec --bind condensed_params_and_db_" + str(i) + "/db.db --bind " + target_molecule_file + " /pi/summer.thyme-umw/enamine-REAL-2.6billion/shapedb_container.sif /pharmit/src/build/shapedb -NNSearch -k 100000 -ligand " + target_molecule_file + " -db condensed_params_and_db_" + str(i) + "/db.db -print > " + working_chunk + "_" + str(i) + "_nn.txt")

	#run through the shapedb list and remove any ligands on the blacklist
	#TODO

	#delete the decompressed folder
	os.system("rm -drf condensed_params_and_db_" + str(i))


