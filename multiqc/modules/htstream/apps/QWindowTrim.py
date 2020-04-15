from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################

class QWindowTrim():


	def table(self, json):

		# Standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["PE in"] = {'namespace': "PE in", 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'namespace': "PE out", 'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'namespace': "SE in",'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Avg. BP Trimmed"] = {'namespace': "Avg. BP Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["Notes"] = {'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, bp):

		# config dictionary for bar graph
		config = {'title': "HTStream: Trimmed Basepairs Bargraph"}

		# returns nothing if no basepairs were trimmed.
		if bp == 0:
			return

		# standard bar graph construction. See MultiQC docs.
		categories  = OrderedDict()

		categories['Left Trimmed Basepairs'] = {'name': 'Left Trimmed Basepairs'}
		categories['Right Trimmed Basepairs'] = {'name': 'Right Trimmed Basepairs'}

		return bargraph.plot(json, categories, config)



	def execute(self, json):

		stats_json = OrderedDict()

		# accumular variable that prevents empty bar graphs
		total_trimmed_bp = 0

		for key in json.keys():

			# trimmed reads by side
			lefttrimmed_bp = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_bp = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total trimmed reads
			trimmed_bp = (lefttrimmed_bp + rightrimmed_bp)

			# sample dictionary entry
			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "SE in" : json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "Avg. BP Trimmed": trimmed_bp / json[key]["Fragment"]["in"],
							   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Left Trimmed Basepairs": lefttrimmed_bp,
							   "Right Trimmed Basepairs": rightrimmed_bp
						 	  }

			# total basepairs accumlation 
			total_trimmed_bp += trimmed_bp


		# sections and figure function calls
		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Basepairs": self.bargraph(stats_json, total_trimmed_bp)
				   }

		return section
		