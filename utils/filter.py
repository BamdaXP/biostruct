import json

from biostruct.entry import Entry

from typing import OrderedDict

__filter_results_dict__ = None


@property
def filter_results_dict(self) -> OrderedDict:
		if self.__filter_results_dict__ == None:
			#load filter results if exist
			try:
				with open(self.get_path('filter_results.json'), 'r') as jf:
					self.__filter_results_dict__ = json.load(jf)
			except Exception:
				self.__filter_results_dict__ = OrderedDict()
		return self.__filter_results_dict__

def set_filters(self, *filters) -> None:
	#init filter dict
	self.filter_results_dict.clear()
	for f in filters:
		self.filter_results_dict[f.__name__] = {'results': [], 'count':0}
	chains = [c for c in self.fasta_dict.keys()]
	#filter pairs and add to dict
	for i in range(0, len(chains)):
		for j in range(i+1, len(chains)):
			pair = (chains[i], chains[j])
			flag = True
			for f in filters:
				try:
					flag = flag and f(self, pair)
				except:
					flag = False
				if flag:
					self.filter_results_dict[f.__name__]['results'].append(
						"{0} {1}".format(pair[0], pair[1]))
					self.filter_results_dict[f.__name__]['count'] = len(
						self.filter_results_dict[f.__name__]['results'])
				else:
					break

	#save filter results
	with open(self.get_path('filter_results.json'), 'w') as jf:
		json.dump(self.filter_results_dict, jf)
