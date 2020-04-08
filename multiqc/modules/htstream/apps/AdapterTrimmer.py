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
		headers["Avg. BP Trimmed"] = {'description': 'Average Number of basepairs trimmed from reads', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():
			
			adapter_reads = json[key]["Single_end"]["adapterTrim"] + json[key]["Paired_end"]["adapterTrim"] # total reads trimmed
			bp_reads = json[key]["Single_end"]["adapterBpTrim"] + json[key]["Paired_end"]["adapterBpTrim"] # total basepairs trimmed

			if adapter_reads == 0:
				perc_adapters = 0
				avg_bp_trimmed = 0

			else:
				perc_adapters = (adapter_reads / json[key]["Fragment"]["in"]) * 100
				avg_bp_trimmed = (bp_reads / adapter_reads)


			stats_json[key] = {
							   "Reads in": json[key]["Fragment"]["in"],
							   "Reads out": json[key]["Fragment"]["out"],
							   "% Adapters": perc_adapters,
							   "Avg. BP Trimmed": avg_bp_trimmed,
							   "Notes": json[key]["Program_details"]["options"]["notes"]
							  }

		section = {"Table": self.table(stats_json)}


		return section