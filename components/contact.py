import os,sys

import numpy as np

from biostruct.entry import Entry,RESIDUE_REF
from biostruct.components.alignment import get_alignments

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
		
		i = 0
		for x in range(0, len(chain1_model)):
			animo1 = chain1_model.child_list[x]
			#Skip non-animos
			if not animo1.resname in RESIDUE_REF:
				continue
			#Skip gaps
			while aln1[0][i] == '-':
				i += 1
			#Skip non-matchings
			if aln1[1][i] != '|':
				i += 1
				continue
			j = 0
			for y in range(0, len(chain2_model)):
				animo2 = chain2_model.child_list[y]
				if not animo2.resname in RESIDUE_REF:
					continue
				while aln2[0][j] == '-':
					j += 1
				if aln2[1][j] != '|':
					j += 1
					continue

				distmap = np.array([[atom1-atom2 for atom2 in animo2 if atom2.element != 'H']for atom1 in animo1 if atom1.element != 'H'])
				d = np.min(distmap)
				map_data[i, j] = 0 if d > cutoff else 1
				
				j += 1
			i += 1
		#save mapdata
		if np.max(map_data) > 0 or save_zeros:
			np.savetxt(contact_filepath, map_data, fmt='%+2i')
			
	return map_data
