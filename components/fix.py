import os

import biostruct as bs
import pandas as df

import modeller 
import modeller.automodel

from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
import Bio.SeqIO
import Bio.PDB.Chain
from Bio.PDB.MMCIFParser import MMCIFParser


def get_fix(entry:bs.Entry,chain_id:str,redo=False)-> Bio.PDB.Chain.Chain:
	#Make fix dir
	os.makedirs(entry.get_path('fix/'),exist_ok=True)
	fix_file = entry.get_path(f'fix/{chain_id}.cif')
	if os.path.isfile(fix_file) and os.path.getsize(fix_file):
		#Export sequence fasta file
		fasta_file = entry.get_path(f'fix/{chain_id}.fasta')
		record = SeqRecord(Seq(entry.fasta_dict[chain_id]['sequence']),
							id="",
							name="",
							description=entry.fasta_dict[chain_id]['description'])
		Bio.SeqIO.write(record, fasta_file, "fasta")

		print(f"Fixing cif file chain {chain_id} in {entry}")
		env = modeller.Environ()
		m = modeller.Model(env,file="mmcif.cif",model_format='MMCIF')
		aln =modeller.alignment(env)
		aln.append(file=fasta_file, align_codes='all')
		a = modeller.automodel.automodel(env,alnfile=aln,knowns=m,sequence=fasta_file)
		a.starting_model = 1
		a.ending_model = 1
		a.chains = chain_id
		a.make()
		a.restraints.write(file=fix_file) #type:ignore

	return MMCIFParser(QUIET=True).get_structure(entry.pid,fix_file)[0][chain_id]