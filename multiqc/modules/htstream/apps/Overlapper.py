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

		config = {'namespace': 'overlapper'}

		# straight forward table construction.
		headers = OrderedDict()

		headers["Ov_PE_in"] = {'title': "PE in", 'namespace': "PE in",'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Ov_PE_out"] = {'title': "PE out", 'namespace': "PE out",'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ov_SE_in"] = {'title': "SE in", 'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Ov_SE_out"] = {'title': "SE out", 'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ov_%_Overlapped"] = {'title': "% Overlapped", 
									  'namespace': "% Overlapped",
									  'description': 'Percentage of Reads with Overlap.',
									  'suffix': '%',
									  'max': 100,
									  'format': '{:,.2f}',
									  'scale': 'Blues'}
		headers["Ov_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, inserts):

		# configuration dictionary for bar graph
		config = {'title': "HTStream: Overlap Composition Bargraph",
				  'id': "htstream_overlapper_bargraph",
				  'ylab' : "Samples"}

		# if no overlaps at all are present, return nothing
		if inserts == 0:
			html = '<div class="alert alert-info"> No overlaps present in samples. </div>'	
			return html

		# bargraph dictionary. Exact use of example in MultiQC docs.
		categories  = OrderedDict()

		categories['Ov_Sins'] = {
							  'name': 'Short Inserts',
							  'color': '#4d8de4'
							 }
		categories['Ov_Mins'] = {
							  'name': 'Medium Inserts',
							  'color': '#e57433'
							 }
		categories['Ov_Lins'] = {
							  'name': 'Long Inserts',
							  'color': '#33a02c'
							 }

		return bargraph.plot(json, categories, config)



	def linegraph(self, json):

		# config dictionary for "density" plots. Its a work in progress. 
		config = {'title': "HTStream: Overlapped Lengths Density Plots",
				  'data_labels': []}

		# initialize data structures
		data_list = [] 

		for key in json.keys():

			# creates empty dictionary to hold data for line graph. 
			data = {}
			data[key] = {}

			# line graph config dictionary
			config_subdict = {'name': key, 'ylab': 'Frequency', 'xlab': 'Overlapping Read Lengths'}	

			# calculate totals for frequency histogran
			total = sum([ count[1] for count in json[key]["Ov_Histogram"] ])

			# iterates over ever value in histogram and adds it to line graph
			for item in json[key]["Ov_Histogram"]:

				data[key][item[0]] = item[1] / total


			# appends data set to data list and config dictionary to data labels section of line graph config
			data_list.append(data)
			config['data_labels'].append(config_subdict)


		return linegraph.plot(data_list, config)


	def execute(self, json):

		stats_json = OrderedDict()

		# accumulator for inserts, used to prevent empty bar graph
		inserts = 0

		for key in json.keys():

			sins = json[key]["Fragment"]["short_inserts"]
			mins = json[key]["Fragment"]["medium_inserts"]
			lins = json[key]["Fragment"]["long_inserts"]

			overlapped_sum = (sins + mins + lins)

			# if input fragements are zero, this will result a warning message and termination of this module
			# 	meant to serve as an example for future reference.
			if json[key]["Fragment"]["in"] == 0:

				file_name = os.path.basename(__file__).split(".")[0]
				warning_message = "No Input Reads for HTStream " + file_name + ". Check file for format errors."
				log.warning(warning_message)
				return ""

			# the INFAMOUS percent overlapped 
			perc_overlapped = ((sins + mins) / json[key]["Paired_end"]["in"]) * 100
				
			# sample instance in dictionary
			stats_json[key] = {
			 			 	   "Ov_PE_in": json[key]["Paired_end"]["in"],
							   "Ov_PE_out": json[key]["Paired_end"]["out"],
							   "Ov_SE_in" : json[key]["Single_end"]["in"],
							   "Ov_SE_out": json[key]["Single_end"]["out"],
							   "Ov_%_Overlapped": perc_overlapped,
						 	   "Ov_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Ov_Sins": sins,
							   "Ov_Mins": mins,
							   "Ov_Lins": lins,
							   "Ov_Histogram": json[key]["Fragment"]["readlength_histogram"]
							  }

			# accumulator accumlating
			inserts += overlapped_sum


		# sections and function calls 
		section = {
				   "Table": self.table(stats_json),
				   "Overlap Composition": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json)
				   }

		return section