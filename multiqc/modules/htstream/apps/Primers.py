from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" Primers submodule for HTStream charts and graphs """

#################################################

class Primers():

	def table(self, json):

		# standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["PE in"] = {'namespace': "PE in", 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'namespace': "PE out", 'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Reads Flipped"] = {'namespace': "Reads Flipped", 'description': 'Number of Flipped Reads', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Notes"] = {'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json):

		# bar graph config dict
		config = {'title': "HTStream: Reads with Primers Bargraph"}

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

			for item in json[sample]["Primer Counts"]: 

				if item[0] != "None" and item[1] == "None":
					data[sample]["Primer 1 Only"] = item[2]

				elif item[0] == "None" and item[1] != "None":
					data[sample]["Primer 2 Only"] = item[2]

				elif item[0] != "None" and item[1] != "None" :
					data[sample]["Both Primers"] = item[2]

			if data[sample] == {}:
				return ""


		return bargraph.plot(data, categories)


	def execute(self, json):

		stats_json = OrderedDict()


		for key in json.keys():

			# dictionary entry for sample
			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "SE in" : json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "Reads Flipped": json[key]["Fragment"]["flipped"],
							   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Primers": json[key]["Program_details"]["primers"],
							   "Primer Counts": json[key]["Fragment"]["primers_counts"]
						 	  }

		# dictionary for sections and figure function calls
		section = {
				   "Table": self.table(stats_json),
				   "Reads with Primers": self.bargraph(stats_json)
				   }

		return section
