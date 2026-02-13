#the purpose of this script is to runthe dehydration script on all cases of there being the placements.tar.gz and raw_scores.csv in the same folder

import os,sys
from pathlib import Path

#arguments
#take in the root of the folder that you want to look down
cleaner_root = sys.argv[1]

#take in the reference pdb to make a skeleton and remove residues per pdb that did not move significantly
reference_pdb = sys.argv[2]

#move to and iterate down the cleaner root
os.chdir(cleaner_root)

for r,d,f in os.walk(cleaner_root):
	#if there are cases where there is a placements.tar.gz and a raw_Scores.csv, then we know that we completed discovery and can try to dehydrate
	for file in f:
		if file == "placements.tar.gz":
			#heck if there is a raw_scores.csv
				if Path(r + "/raw_scores.csv").exists():
					#we have a placements and raw scores directory where we can try to ruthe dehydrator
					#run it in a lsf job and throttle
					print(r + "/raw_scores.csv " + r + "/placements.tar.gz")

					