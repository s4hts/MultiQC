from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" NTrimmer submodule for HTStream charts and graphs """

#################################################

class NTrimmer():


	def table(self, json):

		headers = OrderedDict()

		headers["Reads in"] = {'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Reads out"] = {'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Discarded"] = {
						   'description': 'Percentage of Reads (SE and PE) Discarded',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Oranges'
						  }

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

			discarded_reads = json[key]["Single_end"]["SE_discarded"] + json[key]["Paired_end"]["R1_discarded"] + json[key]["Paired_end"]["R2_discarded"]  

			lefttrimmed_reads = json[key]["Paired_end"]["R1_leftTrim"] + json[key]["Paired_end"]["R2_leftTrim"] + json[key]["Single_end"]["SE_leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["R1_rightTrim"] + json[key]["Paired_end"]["R2_rightTrim"] + json[key]["Single_end"]["SE_rightTrim"]

			trimmed_reads += (lefttrimmed_reads + rightrimmed_reads)

			stats_json[key] = {
			 				   "Reads in": json[key]["totalFragmentsInput"],
							   "Reads out": json[key]["totalFragmentsOutput"],
							   "% Discarded" : (discarded_reads / json[key]["totalFragmentsInput"]) * 100,
							   "Notes": json[key]["Notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads,
							  }

		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_reads)
				   }

		return section
