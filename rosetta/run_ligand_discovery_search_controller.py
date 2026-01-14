#the purpose of this script is to be a controller that runs the run_ligand_discovery_search.py script on a database of test_params directories to perform rosetta ligand discovery search
#this script requires the following inputs, some of which require knowledge of how Rosetta reads in the target pdb to know the indices of residues as Rosetta perceives it:
#the script needs the following
#1. target pdb
#2. an anchor residue  string of all positions to be used as anchors, which is comma separated if there is more than 1 (i.e. 79 for 1 or 11,79,55,403 for multiple); this residue needs to use rosetta indexing, which needs to be determined before running the pipeline; for each anchor residue, a unique job will be made to improve runtime and avoid the risk of clobbering placement files (since files from the same ligand have a risk of being clobbered)
#3. a motifs file (by default, will use the main 1.6M motif library, but if another is specified, will be used instead)
#4. a location to look down where there are test_params folders. the script will look down at all test_params folders from the given location

#additional arguments that we will treat as mandatory (5-7)
#score cutoff overrides for fa_atr, fa_rep, and ddg

#an option to bring in a text file containing arguments to be appended to the discovery arguments for things like a space fill region or to print the space fill matrix (not recommended forr memory and speed)
#dimensions and center for a space fill matrix (which can speed up discovery by filtering away placements that are not close enough to the desired binding region)

#by default, this script will try to run individual jobs for up to 7 days. if this really needs to be changed, I can make an option or something, but jobs shouldn't be running for that long anyway

#imports
import os,sys

#take in arguments
#target pdb
target_pdb = sys.argv[1]
#anchor residue(s)
anchor_residue_string = sys.argv[2]

#break apart the string by commas into a list
anchor_residue_string_list = anchor_residue_string.split(",")

#motifs file (likely want /pi/summer.thyme-umw/enamine-REAL-2.6billion/FINAL_motifs_list_filtered_2_3_2023.motifs unless you know what you are doing)
motifs_file = sys.argv[3]

#discovery directory root, will look down from here for all test params directories and attempt discovery for each anchor residue at the same level as each test_params directory
discovery_directory_root = sys.argv[4]

#atr, rep, ddg cutoffs
atr = sys.argv[5]
rep = sys.argv[6]
ddg = sys.argv[7]

#if there is an extra args file, take it in, so we can pass it down to all discovery jobs
extra_args_file = ""
if len(sys.argv) >= 9:
	extra_args_file = sys.argv[8]

#look over the discovery directory root to identify test_params directories
for r,d,f in os.walk(discovery_directory_root):
	for dire in d:
		#only look at the test_params directories
		if dire == "test_params":

			#store the root
			tp_root = r

			#for each test_params directory, go over the anchor residue string and prepare to run discovery for each unique residue in the list
			for residue in anchor_residue_string_list:
				#go to the root
				os.chdir(tp_root)

				#make a directory for the anchor residue
				os.system("mkdir -p " + str(residue))

				#enter the directory
				os.chdir(str(residue))

				#prepare to run discovery on test params for this residue
				#start the command
				cmd = "bsub -q long -W 168:00 -o output -e error -R \"rusage[mem=8000]\" \"python /pi/summer.thyme-umw/enamine-REAL-2.6billion/umass_chan_REAL-M_platform/rosetta/run_ligand_discovery_search.py "
				#target pdb
				cmd = cmd + target_pdb + " "
				#anchor residue
				cmd = cmd + str(residue) + " "
				#motifs file
				cmd = cmd + motifs_file + " "
				#test_params directory
				cmd = cmd + tp_root + "/test_params/ "
				#atr cutoff, rep cutoff, ddg cutoff
				cmd = cmd + str(atr) + " " + str(rep) + " " + ddg + " "
				#extra args file (if any, and if not will append a blank string)
				cmd = cmd + extra_args_file
				#cap off
				cmd = cmd + "\""

				print(cmd)
				os.system(cmd)

				#add a throttle at 1500 parallel jobs to avoid overflowing the system
				#bsub job throttle to make sure we do not exceed our local limit
				#write the length of the bjobs queue to this current location
				os.system("bjobs | wc -l > bjobs_length.txt")
				job_count = 0
				with open("bjobs_length.txt") as f:
					job_count = int(f.read().strip())
				while job_count > 1500:
					#sleep for 1 second to not overburden the system
					os.system("sleep 1")
					os.system("bjobs | wc -l > bjobs_length.txt")
					with open("bjobs_length.txt") as f:
						job_count = int(f.read().strip())
				#remove the length file to avoid clutter
				os.system("rm bjobs_length.txt")				

