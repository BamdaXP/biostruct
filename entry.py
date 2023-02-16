import os,requests,json,re,wget
from typing import Optional

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqIO import parse

from biostruct.monoer import Monoer
from biostruct.pair import Pair
from biostruct.utils.residue_refs import STD_REF

import numpy as np

class Entry:
	@staticmethod
	def get_entries(root_dir,load_structure = True):
		for pid in os.listdir(root_dir):
			yield Entry(pid,root_dir,load_structure)
	@staticmethod
	def reomve(root_dir,pathname)->None:
		for pid in os.listdir(root_dir):
			entry_dir = os.path.join(root_dir, pid)
			path = os.path.join(entry_dir, pathname) 
			if os.path.isfile(path):
				os.remove(path)
			elif os.path.isdir(path):
				for root, dirs, files in os.walk():  # type: ignore
					for name in files:
						os.remove(os.path.join(root, name))
					for name in dirs:
						os.rmdir(os.path.join(root, name))

	def __init__(self, pid: str, root_dir:str,asm:Optional[int] = None) -> None:
		self.pid = pid.upper()
		if asm:
			self.dir = os.path.join(root_dir,f"{pid}({asm})")
		else:
			self.dir = os.path.join(root_dir, pid)
		self.asm = asm

	def get_path(self, name:str)->str:
		return os.path.join(self.dir, name)

	__structure__ = None
	@property
	def structure(self)->Structure:
		if self.__structure__ == None:
			try:
				self.__structure__ = MMCIFParser(QUIET=True).get_structure(
					self.pid, self.get_path('mmcif.cif'))
			except Exception as e:
				print("{0}:Load model failed.".format(self.pid))
				self.__download_cif__()
			finally:
				self.__model__ = MMCIFParser(QUIET=True).get_structure(
					self.pid, self.get_path('mmcif.cif'))
		return self.__structure__  # type: ignore

	@property
	def model(self)->Model:
		return self.structure.child_list[0]
		
	def __download_cif__(self)->None:
		print(f"Download cif {self}...")
		os.makedirs(self.dir, exist_ok=True)

		if self.asm:
			url = f"https://files.rcsb.org/download/{self.pid}-assembly{self.asm}.cif"
		else:
			url = f"https://files.rcsb.org/download/{self.pid}.cif"
		
		wget.download(url=url,
			out=self.get_path('mmcif.cif'),
			bar=None)#type: ignore
		#os.system(f"wget -O \'{self.get_path('mmcif.cif')}\' {url}")
		'''
		#download gz file
		response = requests.get(url)
		gfilepath = self.get_path('cif.gz')
		with open(gfilepath, 'wb') as f:
			f.write(response.content)
			g_file = gzip.GzipFile(os.path.join(gfilepath))
			#extract gz
			with open(self.get_path('mmcif.cif'), 'wb') as f:
				f.write(g_file.read())
			os.remove(gfilepath)
		'''
	__fasta_dict__ = None
	@property
	def fasta_dict(self)->dict:
		'''
		{chain_id:{"description":description,"sequence":sequence}}
		'''
		if self.__fasta_dict__ == None:
			#proceed fasta
			if not os.path.isfile(self.get_path('sequences.fasta')):
				self.__download_fasta__()
			self.__fasta_dict__ = {}
			for record in parse(self.get_path('sequences.fasta'), 'fasta'):
				chain_id_s = str(record.description.split('|')[1])
				description_s = str(record.description)
				chain_id_s = re.sub(r'Chains?', '', chain_id_s)
				for id in chain_id_s.split(','):
					id = id.strip()
					match = re.match(r'(\w+)\[auth (\w+)\]', id)
					if match:
						self.__fasta_dict__[match.group(2)] = {'description': description_s, 'sequence': str(record.seq)}
					else:
						self.__fasta_dict__[id] = {'description': description_s, 'sequence': str(record.seq)}
		return self.__fasta_dict__
	def __download_fasta__(self):
		print("Download fasta {0}...".format(self))
		os.makedirs(self.dir, exist_ok=True)
		fastafilepath = self.get_path('sequences.fasta')
		wget.download(
			url=f"https://www.rcsb.org/fasta/entry/{self.pid}",
			out=fastafilepath,
			bar=None)  # type: ignore
		


	def get_monoers(self,excepts:list[str]=['nuclein','unknown','not_fasta']):
		if 'not_fasta' in excepts:
			chains = self.fasta_dict.keys()
		else:
			chains = self.model.child_dict.keys()
		for chain_id in chains:
			#except chains not in fasta
			if 'unknown' in excepts and self.has_unknown(chain_id):
				continue
			if 'nuclein' in excepts and self.is_nuclein(chain_id):
				continue

			yield Monoer(f"{self}:{chain_id}")

	def get_pairs(self,self_pairing:bool=False):
		monoers = [m for m in self.get_monoers()]
		for i in range(0,len(monoers)):
			for j in range(i, len(monoers)):
				if i == j and not self_pairing:
					continue
				else:
					yield Pair(f"{monoers[i]}_{monoers[j]}")

	def get_modelseq(self,chain_id:str,residue_ref:dict[str,str]=STD_REF)->str:
		return ''.join([STD_REF[res.resname] for res in self.model.child_dict[chain_id] if res.resname in STD_REF.keys()])

	def is_nuclein(self,chain_id:str)->bool:
		seq = self.fasta_dict[chain_id]["sequence"]
		return {seq}.issubset({"AUTGC"})

	def has_unknown(self,chain_id:str)->bool:
		seq = self.fasta_dict[chain_id]["sequence"]
		return "X" in seq
	
	__uniprot_dict__ = None
	@property
	def uniprot_dict(self)->dict:
		if self.__uniprot_dict__ == None:
			#load uniprots
			try:
				with open(self.get_path("uniprot.json"), 'r') as jf:
					self.__uniprot_dict__ = json.load(jf)
			except Exception as e:
				print("{0}:Load uniprot failed!".format(self))
				self.__download_uniprot__()
		return self.__uniprot_dict__  # type: ignore
		
	def __download_uniprot__(self):  
		print("Download uniprot {0}...".format(self))
		os.makedirs(self.dir, exist_ok=True)
		l_pid = self.pid.lower()
		response = requests.get(
				"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{0}".format(l_pid))
		data = json.loads(response.content)
		self.__uniprot_dict__ = {}
		for uid in data[l_pid]['UniProt']:
			for o in data[l_pid]['UniProt'][uid]['mappings']:
				self.__uniprot_dict__[o['chain_id']] = uid
		with open(self.get_path("uniprot.json"), 'w') as jf:
			json.dump(self.__uniprot_dict__, jf)

	#__pfam_dict__ = None
	@property
	def pfam_dict(self)->dict:
		if self.__pfam_dict__ == None:
			#load pfams
			try:
				with open(self.get_path("pfam.json"), 'r') as jf:
					self.__pfam_dict__ = json.load(jf)
			except Exception as e:
				print("{0}:Load pfam failed!".format(self))
				self.__download_pfam__()
		return self.__pfam_dict__
	def __download_pfam__(self):
		print("Download pfam {0}...".format(self))
		os.makedirs(self.dir, exist_ok=True)
		
		l_pid = self.pid.lower()
		response = requests.get(
				"https://www.ebi.ac.uk/pdbe/api/mappings/pfam/{0}".format(l_pid))
		data = json.loads(response.content)
		#print(data[id]['Pfam'])
		self.__pfam_dict__ = {}
		for pfid in data[l_pid]['Pfam']:
			for o in data[l_pid]['Pfam'][pfid]['mappings']:
				self.__pfam_dict__.setdefault(o['chain_id'], list())
				self.__pfam_dict__[o['chain_id']].append(pfid)
		with open(self.get_path("pfam.json"), 'w') as jf:
			json.dump(self.__pfam_dict__, jf)

	def __str__(self):
		if self.asm:
			return f"{self.pid}({self.asm})"
		else:
			return self.pid
