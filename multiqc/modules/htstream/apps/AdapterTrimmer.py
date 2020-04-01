from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" AdapterTrimmer submodule for HTStream charts and graphs """

#################################################

class AdapterTrimmer():


	def table(self, json):

		headers = OrderedDict()

		headers["Reads in"] = {'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Reads out"] = {'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Adapters"] = {
						   'description': 'Percentage of Reads (SE and PE) with an Adapter',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Blues'
						  }

		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():
			
			adapter_reads = json[key]["Single_end"]["SE_adapterTrim"] + json[key]["Paired_end"]["PE_adapterTrim"]  

			stats_json[key] = {
							   "Reads in": json[key]["totalFragmentsInput"],
							   "Reads out": json[key]["totalFragmentsOutput"],
							   "% Adapters" : (adapter_reads / json[key]["totalFragmentsInput"]) * 100,
							   "Notes": json[key]["Notes"]
							  }

		section = {"Table": self.table(stats_json)}

		return section