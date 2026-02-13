#the purpose of this script is to help reduce memory overhead of placements by taking placement pdbs and only keeping residues that are within 5 angstroms of the placed ligand
#this requires a python environment with pymol2 (actually, maybe just go through and get distances via raw text so we can maintain residue indexing pdb comments)
#this is an individual script, which is meant to be operated in the directory that it is run in or a single command line argument for a path to process can be used
#for the purpose of pipeline cleanup, this will look for a file named placements.tar.gz, unpack the tar, make reduced versions of the placements, and then recompress them
#will also remove extraneous folders from the folder

imports os,sys
#import pymol2

working_location = os.getcwd()

#arguments (optional path to go to operate in)
if len(sys.argv) == 2:
	#get the different location
	working_location = sys.argv[1]

#move to the working location
os.chdir(working_location)

#unzip the placements.tar.gz
os.system("tar -xzf placements.tar.gz")

#move into the directory
os.chdir("placements")

placements_location = os.getcwd()

#iterate over the directory and look for placement pdb files
for r,d,f in os.walk(placements_location):
	#iterate over each pdb file
	for file in f:
		#only work on files in the placements location
		if r == placements_location:

			print(file)

			#if it is the minipose pdb, delete it and continue
			if "minipose" in file:
				os.system("rm -drf " + file)
				continue

			#otherwise, delete the corresponding folder with the file's namesake first
			file_base = file.split(".pdb")[0]

			os.system("rm -drf " + file_base)

			#now, read the file so we can prune away distant residues
			#buffer all lines into a list so we only have to read the file once
			read_file = open(file,"r")

			file_lines = []

			#create a dictionary of residues (line starts with ATOM) with all atoms so we can see if the residue is close enough to any ligand (any res atom in the residue is within 5 angstroms of any lig atom)
			#use index as a key; values will be the atom line, which the xyz coordinates can be derived for distance
			prot_res = {}

			#make a list of ligand atom coordinates (no need to hold the lines since we will write them anyway)
			lig_res_atoms = []

			for line in read_file.readlines():
				file_lines.append(line)

				#if a residue atom, work on adding it
				if line.startswith("ATOM"):
					#add it to the protein residue dictionary
					#check if this is aready in the dicitonary or not and make a new list at the key if it does not exist already
					resnum = line.split()[5]
					if resnum not in prot_res.keys():
						prot_res[resnum] = []

					prot_res[resnum].append(line)

				#if a hetatm (ligand) take the coordinates and add to the list
				if line.startswith("HETATM"):
					#append the x,y,z coords to lig_res_Atoms
					lig_res_atoms.append([float(line.split()[6]),float(line.split()[7]),float(line.split()[8])])

			#we have now buffered in the file, determine which residues to keep based on distance
			#list of residues to keep by index
			res_to_keep = []

			for res in prot_res.keys():
				for res_atom in prot_res[res]:
					#extract the residue atom coordinates
					x = float(res_atom.split()[6])
					y = float(res_atom.split()[7])
					z = float(res_atom.split()[8])

					#get the coordinate distance of the atom to each atom in the ligand
					for lig_atom in lig_res_atoms:
						distance = (((lig_atom[0] - x) ** 2) + ((lig_atom[1] - y) ** 2) + ((lig_atom[2] - z) ** 2)) ** 0.5

						#if the distance is within 5, note the residue in residues to keep
						if distance < 5:
							res_to_keep.append(res)
							#break
							break

					#if we added the residue to the res_to_keep, break here to move to the next residue
					if res in res_to_keep:
						break

			#finally, go back over the buffered file, and when ew get to ATOM lines, determine whether or not to keep the residue based on if the residue is in the res_to_keep
			write_file = ("temp.pdb", "w")

			for line in file_lines:
				#if the line does not start with ATOM, write
				if line.startswith("ATOM") == False:
					write_file.write(line)
				else:
					#if it is a residue atom line, determine the residue index and determine whether we identified if it iwas close enough
					if line.split()[5] in res_to_keep:
						write_file.write(line)

			#close the read and write streams
			read_file.close()
			write_file.close()

			#overwrite the original with the trimmed with the exact same name
			os.system("mv temp " + file)

#now, all files should be processed, recompress
#move back up to the level where the old compressed palcements file is
os.chdir(working_location)

os.system("tar -czf placements.tar.gz placements")





