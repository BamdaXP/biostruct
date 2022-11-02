
from biostruct.monoer import Monoer


'''
{pid1}:{chain1_id}_{pid2}:{chain2_id}
'''
class Pair():
	
	def __init__(self,full_string:str):
		self.full_string:str = full_string
		
	@property
	def full_string(self)->str:  # type: ignore
		return f"{self.monoer1}_{self.monoer2}"

	@full_string.setter
	def full_string(self,fs:str)->None:  # type: ignore
		fss = fs.split('_')
		assert len(fss) == 2
		self.monoer1 = Monoer(fss[0])
		self.monoer2 = Monoer(fss[1])

	def is_same_entry(self)->bool:
		return self.monoer2.pid == self.monoer1.pid
		
	def __hash__(self):
		return hash(self.full_string)

	def __eq__(self,other):
		return isinstance(other, Pair) and (self.full_string == other.full_string or (self.monoer1 == other.monoer2 and self.monoer2 == other.monoer1))

	def __str__(self):
		return self.full_string
