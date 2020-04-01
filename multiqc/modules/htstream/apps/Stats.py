from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import linegraph

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def graph(self, json):

		histograms = ["R1 histogram", "R2 histogram", "SE histogram"]

		config = {'data_labels': [
								  {'name': "R1 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "R2 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "SE histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								 ],
				  'extra_series': [[], [], []]}

		data_list = []

		for i in range(len(histograms)):

			data = {}

			for key in json.keys():

				if len(data.keys()) == 0:
					data[key] = {}

					total = sum([ count[1] for count in json[key][histograms[i]] ])

					for item in json[key][histograms[i]]:

						data[key][item[0]] = item[1] / total

					data_list.append(data)

				else:
					series_dict = {
								   'name': key,
	        					   'data': []
	        					  }

					for item in json[key][histograms[i]]:
						temp = [item[0], (item[1] / total)]
						series_dict['data'].append(temp)

					config['extra_series'][i].append(series_dict)

		return linegraph.plot(data_list, config)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "R1 histogram": json[key]["SE_readlength_histogram"],
			 				   "R2 histogram": json[key]["R1_readlength_histogram"],
			 				   "SE histogram": json[key]["R2_readlength_histogram"]
						 	  }

		section = {
				   "Density Plots": self.graph(stats_json)
				   }

		return section 

