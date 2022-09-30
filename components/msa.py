import os

from biostruct.entry import Entry

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
	
def get_msa(entry:Entry, chain_id, ncpu=8):
	if not os.path.isdir(entry.get_path('msa/')):
		os.makedirs(entry.get_path('msa/'), exist_ok=True)

	#independent fasta
	fasta_filepath = entry.get_path('msa/{0}.fasta'.format(chain_id))
	if not (os.path.exists(fasta_filepath) and os.path.getsize(fasta_filepath)):
		record = SeqRecord(entry.fasta_dict[chain_id]['sequence'],
							id="", description=entry.fasta_dict[chain_id]['description'])
		SeqIO.write(record, fasta_filepath, "fasta")
	
	#required files
	sto_filepath = entry.get_path('msa/{0}.msa'.format(chain_id))
	database = "/mnt/data1/Common_databases/database/uniref100/uniref100.fasta"
	mask_filepath = entry.get_path("msa/{0}.msk".format(chain_id))
	msa_filepath = entry.get_path("msa/{0}.msa.fasta".format(chain_id))
	#load if exist
	if os.path.exists(msa_filepath) and os.path.getsize(msa_filepath):
		#exist: clear caches
		if os.path.exists(sto_filepath):
			os.remove(sto_filepath)
		if os.path.exists(mask_filepath):
			os.remove(mask_filepath)
	else:
		#not exist: calculate and clear caches
		if not(os.path.exists(sto_filepath) and os.path.getsize(sto_filepath)):
			#calculate
			cmd = "jackhmmer -N 5 --cpu " + str(ncpu) + " --noali --incT {0} -A {1} {2} {3} > /dev/null".format(
				len(entry.fasta_dict[chain_id]['sequence'])*0.5,
				sto_filepath,
				fasta_filepath,
				database)
			#run
			os.system(cmd)
		if not(os.path.exists(mask_filepath) and os.path.getsize(mask_filepath)):
			cmd = "esl-alimask --rf-is-mask {0}>{1}".format(sto_filepath, mask_filepath)
			os.system(cmd)
		cmd = "esl-reformat a2m {0}>{1}".format(mask_filepath, msa_filepath)
		os.system(cmd)
		#completed, clear caches
		os.remove(sto_filepath)
		os.remove(mask_filepath)

	return SeqIO.parse(msa_filepath,'fasta')