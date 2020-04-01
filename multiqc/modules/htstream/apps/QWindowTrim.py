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
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
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

		trimmed_reads = 0

		for key in json.keys():

			lefttrimmed_reads = json[key]["Paired_end"]["R1_leftTrim"] + json[key]["Paired_end"]["R2_leftTrim"] + json[key]["Single_end"]["SE_leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["R1_rightTrim"] + json[key]["Paired_end"]["R2_rightTrim"] + json[key]["Single_end"]["SE_rightTrim"]

			trimmed_reads += (lefttrimmed_reads + rightrimmed_reads)

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
							   "Notes": json[key]["Notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads
						 	  }

		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_reads)
				   }

		return section
		