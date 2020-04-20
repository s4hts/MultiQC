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

		headers["Nt_Reads_in"] = {'title': "Reads in", 'namespace': "Reads in", 'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Nt_Reads_out"] = {'title': "Reads out", 'namespace': "Reads out", 'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Nt_Avg_BP_Trimmed"] = {'title': "Avg. BP Trimmed", 'namespace': "Avg. BP Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["Nt_%_Discarded"] = {'title': "% Discarded",
									 'namespace': "% Discarded",
									 'description': 'Percentage of Reads (SE and PE) Discarded',
									 'suffix': '%',
									 'max': 100,
									 'format': '{:,.2f}',
									 'scale': 'Oranges'
									}

		headers["Nt_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, bps):

		# config dict for bar graph
		config = {"title": "HTStream: Trimmed Basepairs Bargraph",
				  'id': "htstream_ntrimmer_bargraph",
				  'ylab' : "Samples"
				  }

		html = ""

		# returns nothing if no reads were trimmed.
		if bps == 0:
			html = '<div class="alert alert-info"> No basepairs were trimmed from any sample. </div>'	
			return html

		# Bar graph constructor. See MultiQC docs.
		categories  = OrderedDict()

		categories["Nt_Left_Trimmed_bps"] = {'name': 'Left Trimmed Basepairs'}
		categories["Nt_Right_Trimmed_bps"] = {'name': 'Right Trimmed Basepairs'}

		return bargraph.plot(json, categories, config)


	def execute(self, json):

		stats_json = OrderedDict()

		# accumulator variable. Used to prevent empty bargraphs 
		trimmed_bps = 0

		for key in json.keys():

			# number ofreads discarded
			discarded_reads = json[key]["Single_end"]["discarded"] + json[key]["Paired_end"]["Read1"]["discarded"] + json[key]["Paired_end"]["Read2"]["discarded"]  
			
			# number of trimmed reads by side
			lefttrimmed_bps = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_bps = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total number of trimmed reads.
			sample_trimmed_bps = (lefttrimmed_bps + rightrimmed_bps)

			# sample entry in stats dictionary
			stats_json[key] = {
			 				   "Nt_Reads_in": json[key]["Fragment"]["in"],
							   "Nt_Reads_out": json[key]["Fragment"]["out"],
							   "Nt_Avg_BP_Trimmed": sample_trimmed_bps / json[key]["Fragment"]["in"],
							   "Nt_%_Discarded" : (discarded_reads / json[key]["Fragment"]["in"]) * 100,
							   "Nt_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Nt_Left_Trimmed_bps": lefttrimmed_bps,
							   "Nt_Right_Trimmed_bps": rightrimmed_bps,
							  }

			trimmed_bps += sample_trimmed_bps 

		# section and figure function calls
		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_bps)
				   }

		return section
