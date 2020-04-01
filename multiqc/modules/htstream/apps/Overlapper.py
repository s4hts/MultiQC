from collections import OrderedDict
import logging
import numpy as np

from multiqc import config
from multiqc.plots import table, bargraph, linegraph

#################################################

""" Overlapper submodule for HTStream charts and graphs """

#################################################

class Overlapper():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json, inserts):

		if inserts == 0:
			return

		categories  = OrderedDict()

		categories['Sins'] = {
							  'name': 'Short Inserts',
							  'color': '#4d8de4'
							 }
		categories['Mins'] = {
							  'name': 'Medium Inserts',
							  'color': '#e57433'
							 }
		categories['Lins'] = {
							  'name': 'Long Inserts',
							  'color': '#33a02c'
							 }

		return bargraph.plot(json, categories)

	def linegraph(self, json):

		data_list = [] 
		config = {'data_labels': []}

		for key in json.keys():

			data = {}
			data[key] = {}
			config_subdict = {'name': key, 'ylab': 'Frequency', 'xlab': 'Overlapping Read Lengths'}				

			total = sum([ count[1] for count in json[key]["Histogram"]])

			for item in json[key]["Histogram"]:

				data[key][item[0]] = item[1] / total

			data_list.append(data)
			config['data_labels'].append(config_subdict)

								# min             max     
		#x = np.linspace(norm.ppf(0.01), norm.ppf(0.99), 100)

		return linegraph.plot(data_list, config)

	def execute(self, json):

		stats_json = OrderedDict()

		inserts = 0

		for key in json.keys():

			stats_json[key] = {
			 			 	   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
						 	   "Notes": json[key]["Notes"],
							   "Sins": json[key]["sins"],
							   "Mins": json[key]["mins"],
							   "Lins": json[key]["lins"],
							   "Histogram": json[key]["readlength_histogram"]
							  }

			inserts += (json[key]["sins"] + json[key]["mins"] + json[key]["lins"])

		section = {
				   "Table": self.table(stats_json),
				   "Reads with Insertions": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json)
				   }

		return section