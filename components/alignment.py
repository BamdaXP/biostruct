import os,sys


import Bio.pairwise2
from Bio.Align import substitution_matrices

from biostruct.entry import Entry


BLOSUM62 = substitution_matrices.load("BLOSUM62")
ALIGNMENT_FILE_PATTERN = r"(\S+)\((\S+)\)_(\S+)\((\S+)\)\.aln"

def get_alignments(entry:Entry, chain_pair: tuple, source_pair: tuple ,use_matrix:bool=False,redo = False) -> str:
	"""
	entry: Target entry
	chain_pair: Tuple of chain id (A,B)
	source_pair: Tuple of chain source: (fasta,model)
		- model: from structure 
		- fasta: from fasta sequence file
	redo: Whether force to recalculate
	"""
	#make dir if not exsit
	if not os.path.isdir(entry.get_path('alignments/')):
		os.makedirs(entry.get_path('alignments/'), exist_ok=True)

	filename = entry.get_path(os.path.join("alignments/", f"{chain_pair[0]}({source_pair[0]})_{chain_pair[1]}({source_pair[1]}).aln"))
	if os.path.isfile(filename) and not redo:
		with open(filename, 'r') as f:
			return f.read()

	else:
		#calculate and save if not exist
		print(f"Calculating alignments {entry}:{chain_pair} {source_pair} ...",file=sys.stdout)
		seqs = []
		for i in range(0,len(chain_pair)):
			if source_pair[i] == 'model':
				seq = entry.get_modelseq(chain_pair[i])
			elif source_pair[i] == 'fasta':
				seq = entry.fasta_dict[chain_pair[i]]['sequence']
			else:
				raise Exception("Not supported sequence source!")
			seqs.append(seq)
			
		try:
			alignments = Bio.pairwise2.align.localds(seqs[0], seqs[1], BLOSUM62, -10, -1)  # type: ignore
		except:
			alignments = Bio.pairwise2.align.localxx(seqs[0], seqs[1])  # type: ignore

		s = Bio.pairwise2.format_alignment(*alignments[0], full_sequences=True)
		with open(filename, 'w') as f:
			f.write(s)
		return s
