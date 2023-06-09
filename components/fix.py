import os,re

from biostruct.entry import Entry

import pandas as df

import modeller 
import modeller.automodel

from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
import Bio.SeqIO
import Bio.PDB.Chain
from Bio.PDB.MMCIFParser import MMCIFParser


def get_fix(entry:Entry,chain_id:str,redo=False)-> Bio.PDB.Chain.Chain:
	#Make fix dir
	os.makedirs(entry.get_path('fix/'),exist_ok=True)
	fix_file = entry.get_path(f'fix/{chain_id}.cif')
	if not (os.path.isfile(fix_file) and os.path.getsize(fix_file)) or redo:
		#Export sequence fasta file
		fasta_file = entry.get_path(f'fix/{chain_id}.fasta')
		record = SeqRecord(Seq(entry.fasta_dict[chain_id]['sequence']),
							id="Sequence",
							name='',
							description='')
		Bio.SeqIO.write(record, fasta_file, "fasta")

		print(f"Fixing cif file chain {chain_id} in {entry}")
		env = modeller.Environ()

		m = modeller.Model(env,file=entry.get_path("mmcif.cif"),model_format='MMCIF', model_segment=(f'FIRST:{chain_id}',f'LAST:{chain_id}'))
		aln =modeller.Alignment(env)
		aln.append(file=fasta_file, align_codes="Sequence",alignment_format='FASTA')
		aln.append_model(mdl=m, align_codes="Structure",atom_files=entry.get_path("mmcif.cif"))
		aln.align2d()
		aln.write(file=entry.get_path(f'fix/{chain_id}.pir'),alignment_format='PIR')
		#Clean \n in pir headers
		s = ''
		with open(entry.get_path(f'fix/{chain_id}.pir'),'r') as f:
			s = f.read()
			s = ':'.join([[parts[0]] + [part if '>P1;' in part else part.replace('\n', ' ')
		                  for part in parts[1:-1]] + [parts[-1]] for parts in [s.split(':')]][0])
		with open(entry.get_path(f'fix/{chain_id}.pir'),'w') as f:
			f.write(s)

		a = modeller.automodel.AutoModel(env, alnfile=entry.get_path(f'fix/{chain_id}.pir'), knowns="Structure", sequence="Sequence",
		root_name=f"{entry}_{chain_id}")
		a.starting_model = 1
		a.ending_model = 1
		a.make()
		a.write(entry.get_path(f'fix/{chain_id}.cif'),model_format='MMCIF')

		#delete cache files
		os.remove(f"{a.root_name}.ini")
		os.remove(f"{a.root_name}.rsr")
		os.remove(f"{a.root_name}.sch")
		os.remove(f"{a.root_name}.D00000001")
		os.remove(f"{a.root_name}.B99990001.pdb")
		os.remove(f"{a.root_name}.V99990001")

	return MMCIFParser(QUIET=True).get_structure(entry.pid,fix_file)[0]["A"]#The generated cif file would change chain id to A