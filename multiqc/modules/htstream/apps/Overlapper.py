from collections import OrderedDict
import logging
import os

from multiqc import config
from multiqc.plots import table, bargraph, linegraph

#################################################

""" Overlapper submodule for HTStream charts and graphs """

#################################################

log = logging.getLogger(__name__)

class Overlapper():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Overlapped"] = {'description': 'Percentage of Reads with Overlap.',
								   'suffix': '%',
								   'max': 100,
								   'format': '{:,.2f}',
								   'scale': 'Blues'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json, inserts):

		config = {'title': "HTStream: Overlap Composition Bargraph"}

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

		return bargraph.plot(json, categories, config)

	def linegraph(self, json):

		data_list = [] 
		config = {'title': "HTStream: Overlapped Lengths Density Plots",
				  'data_labels': []}

		for key in json.keys():

			data = {}
			data[key] = {}
			config_subdict = {'name': key, 'ylab': 'Frequency', 'xlab': 'Overlapping Read Lengths'}	

			total = sum([ count[1] for count in json[key]["Histogram"] ])

			sample_list = []

			for item in json[key]["Histogram"]:

				data[key][item[0]] = item[1] / total


			data_list.append(data)
			config['data_labels'].append(config_subdict)


		return linegraph.plot(data_list, config)


	def execute(self, json):

		stats_json = OrderedDict()

		inserts = 0

		for key in json.keys():

			sins = json[key]["Fragment"]["short_inserts"]
			mins = json[key]["Fragment"]["medium_inserts"]
			lins = json[key]["Fragment"]["long_inserts"]

			overlapped_sum = (sins + mins + lins)

			if json[key]["Fragment"]["in"] == 0:

				file_name = os.path.basename(__file__).split(".")[0]
				warning_message = "No Input Reads for HTStream " + file_name + ". Check file for format errors."
				log.warning(warning_message)
				return

			else:
				perc_overlapped = (overlapped_sum / json[key]["Fragment"]["in"]) * 100
				

			stats_json[key] = {
			 			 	   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "SE in" : json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "% Overlapped": perc_overlapped,
						 	   "Notes": json[key]["Program_details"]["options"]["notes"],
							   "Sins": sins,
							   "Mins": mins,
							   "Lins": lins,
							   "Histogram": json[key]["Fragment"]["readlength_histogram"]
							  }

			inserts += overlapped_sum

		section = {
				   "Table": self.table(stats_json),
				   "Overlap Composition": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json)
				   }

		return section