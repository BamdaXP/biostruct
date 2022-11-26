import os,re
from Bio.SeqIO import write
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def cdhit_cluster(seqs:list[str],similarity_cutoff:float,thread_num:int=8)->list[int]:
	seqrecords = []
	for i in range(0,len(seqs)):
		seqrecords.append(SeqRecord(seq=Seq(seqs[i]),id=str(i),name='',description=''))
	write(seqrecords,"temp.fasta","fasta")
	#{os.path.join(os.path.dirname(__file__),'../cdhit/cd-hit')}
	cmd = f"biostruct/cdhit/cd-hit -i temp.fasta -o temp.cluster -c {similarity_cutoff} -T {thread_num} -n 2 -d 0"
	os.system(cmd)

	clusters = [None for i in range(0,len(seqs))]
	with open('temp.cluster.clstr') as f:
		cluster = None
		for line in f.readlines():
			m = re.match(r'>Cluster (\d+)',line)
			if m:
				cluster = m.group(1) 
				continue

			m = re.match(r'.+>(\d+)\.\.\.',line)
			if m:
				seq_id = int(m.group(1))
				assert cluster
				clusters[seq_id] = cluster# type: ignore 
				continue
	assert not None in clusters
	return clusters# type: ignore