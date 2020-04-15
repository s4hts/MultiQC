from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" NTrimmer submodule for HTStream charts and graphs """

#################################################

class NTrimmer():


	def table(self, json):

		# Table construction. Taken from MultiQC docs.

		headers = OrderedDict()

		headers["Reads in"] = {'namespace': "Reads in", 'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Reads out"] = {'namespace': "Reads out", 'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Avg. BP Trimmed"] = {'namespace': "Avg. BP Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["% Discarded"] = {
						   'namespace': "% Discarded",
						   'description': 'Percentage of Reads (SE and PE) Discarded',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Oranges'
						  }

		headers["Notes"] = {'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, reads):

		# config dict for bar graph
		config = {"title": "HTStream: Trimmed Reads Bargraph"}

		# returns nothing if no reads were trimmed.
		if reads == 0:
			return

		# Bar graph constructor. See MultiQC docs.
		categories  = OrderedDict()

		categories['Left Trimmed Reads'] = {'name': 'Left Trimmed Reads'}
		categories['Right Trimmed Reads'] = {'name': 'Right Trimmed Reads'}

		return bargraph.plot(json, categories, config)


	def execute(self, json):

		stats_json = OrderedDict()

		# accumulator variable. Used to prevent empty bargraphs 
		trimmed_reads = 0

		for key in json.keys():

			# number ofreads discarded
			discarded_reads = json[key]["Single_end"]["discarded"] + json[key]["Paired_end"]["Read1"]["discarded"] + json[key]["Paired_end"]["Read2"]["discarded"]  
			
			# number of trimmed reads by side
			lefttrimmed_reads = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total number of trimmed reads.
			trimmed_reads += (lefttrimmed_reads + rightrimmed_reads)

			# sample entry in stats dictionary
			stats_json[key] = {
			 				   "Reads in": json[key]["Fragment"]["in"],
							   "Reads out": json[key]["Fragment"]["out"],
							   "Avg. BP Trimmed": trimmed_reads / json[key]["Fragment"]["in"],
							   "% Discarded" : (discarded_reads / json[key]["Fragment"]["in"]) * 100,
							   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads,
							  }

		# section and figure function calls
		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_reads)
				   }

		return section
