from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" Primers submodule for HTStream charts and graphs """

#################################################

class Primers():

	def table(self, json, index):

		# standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["Pr_PE_in" + index] = {'title': "PE in", 'namespace': "PE in", 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Pr_PE_out" + index] = {'title': "PE out", 'namespace': "PE out", 'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Pr_SE_in" + index] = {'title': "SE in", 'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Pr_SE_out" + index] = {'title': "SE out", 'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Pr_Reads_Flipped" + index] = {'title': "Reads Flipped", 'namespace': "Reads Flipped", 'description': 'Number of Flipped Reads', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Pr_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json):

		# bar graph config dict
		config = {'title': "HTStream: Reads with Primers Bargraph",
				  'id': "htstream_primers_bargraph"}

		# bar graph constuctor
		categories  = OrderedDict()

		categories['Primer 1 Only'] = {
									   'name': 'Primer 1 Only',
									   'color': '#4d8de4'
									  }
		categories['Primer 2 Only'] = {
									   'name': 'Primer 2 Only',
									   'color': '#e57433'
									  }
		categories['Both Primers'] = {
									   'name': 'Both Primers',
									   'color': '#33a02c'
									  }

		# data dictionary for bar graph
		data  = OrderedDict()


		'''
		##################
		This tool is currently a work in progress 
		##################
		'''
		for sample in json.keys():

			data[sample] = {}

			for item in json[sample]["Pr_Primer_Counts"]: 

				if item[0] != "None" and item[1] == "None":
					data[sample]["Primer 1 Only"] = item[2]

				elif item[0] == "None" and item[1] != "None":
					data[sample]["Primer 2 Only"] = item[2]

				elif item[0] != "None" and item[1] != "None" :
					data[sample]["Both Primers"] = item[2]

			if data[sample] == {}:
				return ""


		return bargraph.plot(data, categories)


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}


		for key in json.keys():

			overview_dict[key] = {"Output_Bp": json[key]["Fragment"]["basepairs_out"]}

			# dictionary entry for sample
			stats_json[key] = {
			 				   "Pr_PE_in" + index: json[key]["Paired_end"]["in"],
							   "Pr_PE_out" + index: json[key]["Paired_end"]["out"],
							   "Pr_SE_in" + index: json[key]["Single_end"]["in"],
							   "Pr_SE_out" + index: json[key]["Single_end"]["out"],
							   "Pr_Reads_Flipped" + index: json[key]["Fragment"]["flipped"],
							   "Pr_Notes" + index: json[key]["Program_details"]["options"]["notes"],
							   "Pr_Primers": json[key]["Program_details"]["primers"],
							   "Pr_Primer_Counts": json[key]["Fragment"]["primers_counts"]
						 	  }

		# dictionary for sections and figure function calls
		section = {
				   "Table": self.table(stats_json, index),
				   "Reads with Primers": self.bargraph(stats_json),
				   "Overview": overview_dict
				   }

		return section
