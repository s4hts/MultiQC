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

		config = {'xlab': "Reads", 'ylab': "Duplicates",
				  'extra_series': []}
		data = {}

		for key in json.keys():

			if len(data.keys()) == 0:
				data[key] = {}

				for item in json[key]["Saturation"]:

					data[key][item[0]] = item[1] 

			else:
				series_dict = {
							   'name': key,
        					   'data': []
        					  }

				for item in json[key]["Saturation"]:
					series_dict['data'].append(item)

				config['extra_series'].append(series_dict)

		return linegraph.plot(data, config)


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
				   "Saturation Plot": self.linegraph(stats_json)
				   }

		return section