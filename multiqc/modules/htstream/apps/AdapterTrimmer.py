from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" AdapterTrimmer submodule for HTStream charts and graphs """

#################################################

class AdapterTrimmer():

	def table(self, json):

		# Table constructor. Just like the MultiQC docs.
		
		headers = OrderedDict()

		headers["At_Reads_in"] = {'title': "Reads in", 'namespace': "Reads in", 'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["At_Reads_out"] = {'title': "Reads out", 'namespace': "Reads out", 'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["At_%_Adapters"] = {'title': "% Adapters",
									'namespace': "% Adapters",
									'description': 'Percentage of Reads (SE and PE) with an Adapter',
									'suffix': '%',
									'max': 100,
									'format': '{:,.2f}',
									'scale': 'Blues'
									}
		headers["At_Avg_BP_Trimmed"] = {'title': "Avg. BP Trimmed", 'namespace': "Avg. BP Trimmed", 'description': 'Average Number of basepairs trimmed from reads', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["At_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():
			
			# calculations for reads with adapters and bps trimmed
			adapter_reads = json[key]["Single_end"]["adapterTrim"] + json[key]["Paired_end"]["Read1"]["adapterTrim"] + json[key]["Paired_end"]["Read2"]["adapterTrim"] # total reads trimmed
			bp_reads = json[key]["Single_end"]["adapterBpTrim"] + json[key]["Paired_end"]["Read1"]["adapterBpTrim"] + json[key]["Paired_end"]["Read2"]["adapterBpTrim"] # total basepairs trimmed

			# if adapter trim is zero, so is the percentage and the avg basepair trimmed. This prevents division by zero error
			if adapter_reads == 0:
				perc_adapters = 0
				avg_bp_trimmed = 0

			else:
				perc_adapters = (adapter_reads / json[key]["Fragment"]["in"]) * 100
				avg_bp_trimmed = (bp_reads / adapter_reads)

			# sample dictionary entry
			stats_json[key] = {
							   "At_Reads_in": json[key]["Fragment"]["in"],
							   "At_Reads_out": json[key]["Fragment"]["out"],
							   "At_%_Adapters": perc_adapters,
							   "At_Avg_BP_Trimmed": avg_bp_trimmed,
							   "At_Notes": json[key]["Program_details"]["options"]["notes"]
							  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json)}


		return section