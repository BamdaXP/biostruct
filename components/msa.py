import os

from biostruct.entry import Entry

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

def get_msa(entry:Entry, chain_id:str, db:str, ncpu=4, redo=False):
	if not os.path.isdir(entry.get_path(f'msa/{db}')):
		os.makedirs(entry.get_path(f'msa/{db}'), exist_ok=True)

	#independent fasta
	fasta_filepath = entry.get_path(f'msa/{db}/{chain_id}.fasta')
	#if not (os.path.exists(fasta_filepath) and os.path.getsize(fasta_filepath)):
	record = SeqRecord(Seq(entry.fasta_dict[chain_id]['sequence']),
						id="",
						name="",
						description=entry.fasta_dict[chain_id]['description'])
	
	SeqIO.write(record, fasta_filepath, "fasta")
	
	#required files
	database = f"/mnt/data1/Common_databases/database/{db}/{db}.fasta"
	sto_filepath = entry.get_path(f'msa/{db}/{chain_id}.msa')
	mask_filepath = entry.get_path(f'msa/{db}/{chain_id}.msk')
	msa_filepath = entry.get_path(f'msa/{db}/{chain_id}.msa.fasta')
	#load if exist
	#not exist: calculate and clear caches
	if not(os.path.exists(sto_filepath) and os.path.getsize(sto_filepath)) or redo:
		#calculate
		print(f"Searching msa {entry}:{chain_id}")
		cmd = f"jackhmmer --cpu {ncpu}  --noali -A \"{sto_filepath}\" \"{fasta_filepath}\" \"{database}\" > /dev/null"
		#cmd = f"jackhmmer -N 5 --cpu {ncpu}  --noali --incT {len(entry.fasta_dict[chain_id]['sequence'])*0.5} -A \"{sto_filepath}\" \"{fasta_filepath}\" \"{database}\" > /dev/null"
		#run
		os.system(cmd)
	if not(os.path.exists(mask_filepath) and os.path.getsize(mask_filepath)) or redo:
		print(f"Masking msa {entry}:{chain_id}")
		cmd = f"esl-alimask --rf-is-mask \"{sto_filepath}\">\"{mask_filepath}\""
		os.system(cmd)
	if not (os.path.exists(msa_filepath) and os.path.getsize(msa_filepath)) or redo:
		print(f"Formating msa to fasta {entry}:{chain_id}")
		cmd = f"esl-reformat a2m \"{mask_filepath}\">\"{msa_filepath}\""
		os.system(cmd)
		#completed, clear caches
		#os.remove(sto_filepath)
		#os.remove(mask_filepath)

	return SeqIO.parse(msa_filepath,'fasta')