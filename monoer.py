import re

class Monoer():

	PATTERN = r'(\w+)(\(\d+\))?:(\w+)'

	def __init__(self,full_string:str) -> None:
		self.full_string:str = full_string


	@property
	def full_string(self):  # type: ignore
		if self.asm:
			return f"{self.pid}({self.asm}):{self.chain_id}"
		else:
			return f"{self.pid}:{self.chain_id}"

	@full_string.setter
	def full_string(self, fs):  # type: ignore
		match = re.match(Monoer.PATTERN, fs)
		assert match
		
		group = match.groups()
		self.pid = group[0]
		self.asm = int(group[1][1:-1]) if group[1] else None
		self.chain_id = group[2]

	def __str__(self):
		return self.full_string

	def __eq__(self, other) -> bool:
		return isinstance(other,Monoer) and self.full_string == other.full_string

	def __hash__(self) -> int:
		return hash(self.full_string)