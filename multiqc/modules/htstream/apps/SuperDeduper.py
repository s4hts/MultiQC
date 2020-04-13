from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, linegraph

#################################################

""" SuperDeduper submodule for HTStream charts and graphs """

#################################################

class SuperDeduper():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Duplicates"] = {
						   'description': 'Percentage of Duplicate Reads (SE and PE)',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Oranges'
						  }
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def linegraph(self, json):

		config = {'title': "HTStream: Duplicate Saturation",
				  'xlab': "Reads", 'ylab': "Duplicates",
				  'extra_series': []}

		data = {}
		invariant_saturation_list = []

		html = ""

		for key in json.keys():

			if len(json[key]["Saturation"]) == 1:
				info = key.strip() + ": [" + str(json[key]["Saturation"][0][0]) + " Reads, " + str(json[key]["Saturation"][0][1]) + " Dups ]"
				invariant_saturation_list.append(info)

			else:

				data[key] = {}

				for item in json[key]["Saturation"]:

					data[key][item[0]] = item[1] 


		if len(invariant_saturation_list) != 0:
			
			# to include when more is known about handling invariance.
			notice = "<br />".join(invariant_saturation_list)
			html += str("<p>" + "TEMPORARY PLACE HOLDER FOR INVARIANT SATURATION PLOTS <br />" + notice + "</p>")
			

		if data != {}:
			html += linegraph.plot(data, config)

		return html


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			duplicates = (json[key]["Fragment"]["duplicate"] / json[key]["Fragment"]["in"]) * 100

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "SE in": json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "% Duplicates": duplicates,
							   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Saturation": json[key]["Fragment"]["duplicate_saturation"]
						 	  }

		section = {
				   "Table": self.table(stats_json),
				   "Duplicate Saturation": self.linegraph(stats_json)
				   }

		return section