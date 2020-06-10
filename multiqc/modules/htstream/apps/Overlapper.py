from collections import OrderedDict
import logging
import os, statistics

from multiqc import config
from multiqc.plots import table, bargraph, linegraph

#################################################

""" Overlapper submodule for HTStream charts and graphs """

#################################################

log = logging.getLogger(__name__)

class Overlapper():

	def table(self, json, index):

		config = {'namespace': 'overlapper'}

		# straight forward table construction.
		headers = OrderedDict()

		headers["Ov_PE_loss" + index] = {'title': "% PE Lost", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.2f}', 
								 'suffix': '%', 'scale': 'Greens' }
		headers["Ov_%_Overlapped" + index] = {'title': "% Overlapped", 
									  'namespace': "% Overlapped",
									  'description': 'Percentage of Reads with Overlap.',
									  'suffix': '%',
									  'format': '{:,.2f}',
									  'scale': 'Blues'}
		headers["Ov_SE_in" + index] = {'title': "SE in", 'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Ov_SE_out" + index] = {'title': "SE out", 'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ov_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, inserts):

		# configuration dictionary for bar graph
		config = {'title': "HTStream: Overlap Composition Bargraph",
				  'id': "htstream_overlapper_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False}

		# if no overlaps at all are present, return nothing
		if inserts == 0:
			html = '<div class="alert alert-info"> No overlaps present in samples. </div>'	
			return html

		if len(json.keys()) > 150:
			html = '<div class="alert alert-info"> Too many samples for bargraph. </div>'	
			return html

		# bargraph dictionary. Exact use of example in MultiQC docs.
		categories  = OrderedDict()

		categories['Ov_Sins'] = {
							  'name': 'Short Inserts',
							  'color': '#779BCC'
							 }
		categories['Ov_Mins'] = {
							  'name': 'Medium Inserts',
							  'color': '#C3C3C3'
							 }
		categories['Ov_Lins'] = {
							  'name': 'Long Inserts',
							  'color': '#D1ADC3'
							 }

		return bargraph.plot(json, categories, config)



	def linegraph(self, json):

		# config dictionary for "density" plots. Its a work in progress. 
		config = {'title': "HTStream: Overlapped Lengths",
				  'ylab': "Counts", "xlab": "Overlap Lengths"}

		# initialize data structures
		multi_line = {}

		for key in json.keys():

			# creates empty dictionary to hold data for line graph. 
			multi_line[key] = {}

			# iterates over ever value in histogram and adds it to line graph
			for item in json[key]["Ov_Histogram"]:

				multi_line[key][item[0]] = item[1]



		return linegraph.plot(multi_line, config)


	def parse_histogram_stats(self, hist):

		hist_stats = {"Max": 0,
					  "Median": 0}

		median_list = []

		for item in hist:
			median_list.append(item[0])

			if hist_stats["Max"] < item[1]:
				hist_stats["Max"] = item[0]


		hist_stats["Median"] = statistics.median(median_list)

		return hist_stats


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		# accumulator for inserts, used to prevent empty bar graph
		inserts = 0

		for key in json.keys():

			sins = json[key]["Fragment"]["inserts"]["short"]
			mins = json[key]["Fragment"]["inserts"]["medium"]
			lins = json[key]["Fragment"]["inserts"]["long"]

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
			perc_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"]) * 100
				
			parsed_hist_stats = self.parse_histogram_stats(json[key]["Fragment"]["overlap_histogram"])

			overview_dict[key] = {
								  "Output_Reads": json[key]["Fragment"]["out"],
								  "Hist_Max": parsed_hist_stats["Max"],
								  "Hist_Med": parsed_hist_stats["Median"],
								  "PE_Lost": json[key]["Paired_end"]["out"] / json[key]["Paired_end"]["in"],
								  "SE_Lost": json[key]["Single_end"]["out"] / json[key]["Paired_end"]["in"],
								  "Bp_Lost":  json[key]["Fragment"]["basepairs_out"] / json[key]["Fragment"]["basepairs_in"],
								  "Sin": sins / json[key]["Fragment"]["in"],
								  "Min": mins / json[key]["Fragment"]["in"],
								  "Lin": lins / json[key]["Fragment"]["in"]
								 }

			# sample instance in dictionary
			stats_json[key] = {
							   "Ov_PE_loss" + index: perc_loss,
							   "Ov_SE_in" + index: json[key]["Single_end"]["in"],
							   "Ov_SE_out" + index: json[key]["Single_end"]["out"],
							   "Ov_%_Overlapped" + index: perc_overlapped,
						 	   "Ov_Notes" + index: json[key]["Program_details"]["options"]["notes"],
							   "Ov_Sins": sins,
							   "Ov_Mins": mins,
							   "Ov_Lins": lins,
							   "Ov_Histogram": json[key]["Fragment"]["overlap_histogram"]
							  }

			# accumulator accumlating
			inserts += overlapped_sum


		# sections and function calls 
		section = {"Table": self.table(stats_json, index),
				   "Overlap Composition": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json),
				   "Overview": overview_dict}
				   
		return section