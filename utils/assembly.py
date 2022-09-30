import requests

from typing import Iterator

def check_assemblies(pid:str)->Iterator[int]:
	asm = 1
	while True:
		url = f"https://files.rcsb.org/download/{pid}-assembly{asm}.cif"
		r = requests.head(url)
		if r.status_code == 200: #OK
			yield asm
			asm += 1
		else:
			break