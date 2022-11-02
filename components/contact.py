import os,sys

import numpy as np

from biostruct.entry import Entry
from biostruct.components.alignment import get_alignments
from biostruct.utils.residue_refs import STD_REF

CONTACT_FILENAME_PATTERN = r"(\S+)_(\S+)\(\d+\)\.contact"

def get_contact(entry:Entry, pair:tuple, cutoff:int, save_zeros:bool=False, redo=False) -> np.ndarray:
	#make dir if not exsit
	if not os.path.isdir(entry.get_path('contacts/')):
			os.makedirs(entry.get_path('contacts/'), exist_ok=True)

	(chain1_id,chain2_id) = pair
	contact_filepath = entry.get_path(
		f'contacts/{chain1_id}_{chain2_id}({cutoff}).contact')
	map_data = None
	#load if already exist
	if os.path.isfile(contact_filepath) and not redo:
		map_data = np.loadtxt(contact_filepath)
	else:
		#calculate and save if not exist
		print(f"Calculating contacts {entry}:{pair} ...",file=sys.stdout)
		chain1_model = entry.model[chain1_id]
		chain2_model = entry.model[chain2_id]

		chain1_fasta = entry.fasta_dict[chain1_id]["sequence"]
		chain2_fasta = entry.fasta_dict[chain2_id]["sequence"]

		aln1 = get_alignments(
			entry,(chain1_id,chain1_id),("model","fasta")).splitlines()
		aln2 = get_alignments(
			entry, (chain2_id, chain2_id), ("model", "fasta")).splitlines()

		map_data = np.zeros(shape=(len(chain1_fasta), len(chain2_fasta)))
		map_data[:,:] = -1
		
		aln_index1 = 0
		map_index1 = 0
		for animo1 in chain1_model.child_list:
			#Skip non-animos,only move pdb
			if not animo1.resname in STD_REF:
				continue
			#Incase aln is at end but pdb still has some tails
			if aln_index1 >= len(aln1[2]):
				break
			#Skip gaps in fasta, not move map index
			if aln1[2][aln_index1] == '-':
				aln_index1 += 1
				continue
			#Skip gaps in pdb, not move pdb
			while aln1[0][aln_index1] == '-':
				aln_index1 += 1
				map_index1 += 1

			#Skip non-matchings
			if aln1[1][aln_index1] != '|':
				aln_index1 += 1
				map_index1 += 1
				continue

			aln_index2 = 0
			map_index2 = 0
			for animo2 in chain2_model.child_list:
				#Skip non-animos,only move pdb
				if not animo2.resname in STD_REF:
					continue
				if aln_index2 >= len(aln2[2]):
					break
				#Skip gaps in fasta, not move map index
				if aln2[2][aln_index2] == '-':
					aln_index2 += 1
					continue
				#Skip gaps in pdb, not move pdb
				while aln2[0][aln_index2] == '-':
					aln_index2 += 1
					map_index2 += 1

				#Skip non-matchings
				if aln2[1][aln_index2] != '|':
					aln_index2 += 1
					map_index2 += 1
					continue

				distmap = np.array([[atom1-atom2 for atom2 in animo2 if atom2.element != 'H']
				                   for atom1 in animo1 if atom1.element != 'H'])
				d = np.min(distmap)
				map_data[map_index1, map_index2] = 0 if d > cutoff else 1

				aln_index2 += 1
				map_index2 += 1

			aln_index1 += 1
			map_index1 += 1
			
		#save mapdata
		if np.max(map_data) > 0 or save_zeros:
			np.savetxt(contact_filepath, map_data, fmt='%+2i')
			
	return map_data
