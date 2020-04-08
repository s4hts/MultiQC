from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################

class QWindowTrim():


	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["R1 L/R Ratio"] = {'description': 'Ratio of basepairs trimmed from left to basepairs trimmed from right for Read 1. Pseudocounting applied if R = 0 occurs.', 
								   'format': '{:,.2f}', 'scale': 'Blues'}
		headers["R2 L/R Ratio"] = {'description': 'Ratio of basepairs trimmed from left to basepairs trimmed from right for Read 2. Pseudocounting applied if R = 0 occurs.', 
								   'format': '{:,.2f}', 'scale': 'YlOrRd'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Avg. BP Trimmed"] = {'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def bargraph(self, json, reads):

		if reads == 0:
			return

		categories  = OrderedDict()

		categories['Left Trimmed Reads'] = {'name': 'Left Trimmed Reads'}
		categories['Right Trimmed Reads'] = {'name': 'Right Trimmed Reads'}

		return bargraph.plot(json, categories)

	def execute(self, json):

		stats_json = OrderedDict()

		total_trimmed_reads = 0

		for key in json.keys():

			lefttrimmed_reads = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			trimmed_reads = (lefttrimmed_reads + rightrimmed_reads)

			# Potential creative solution or potential issue, let's see how sam takes it.
			if json[key]["Paired_end"]["Read1"]["rightTrim"] == 0:
				r1_bp_ratio = (json[key]["Paired_end"]["Read1"]["leftTrim"] + 1) / (json[key]["Paired_end"]["Read1"]["rightTrim"] + 1)
			else:
				r1_bp_ratio = json[key]["Paired_end"]["Read1"]["leftTrim"] / json[key]["Paired_end"]["Read1"]["rightTrim"]

			if json[key]["Paired_end"]["Read2"]["rightTrim"] == 0:
				r2_bp_ratio = (json[key]["Paired_end"]["Read2"]["leftTrim"] + 1) / (json[key]["Paired_end"]["Read2"]["rightTrim"] + 1)
			else:
				r2_bp_ratio = json[key]["Paired_end"]["Read2"]["leftTrim"] / json[key]["Paired_end"]["Read2"]["rightTrim"]

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "R1 L/R Ratio": r1_bp_ratio,
							   "R2 L/R Ratio": r2_bp_ratio,
							   "SE in" : json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "Avg. BP Trimmed": trimmed_reads / json[key]["Fragment"]["in"],
							   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads
						 	  }

			total_trimmed_reads += trimmed_reads 

		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, total_trimmed_reads)
				   }

		return section
		