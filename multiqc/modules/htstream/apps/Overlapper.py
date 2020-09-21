from collections import OrderedDict
import logging
import os, statistics

from multiqc import config
from multiqc.plots import table, bargraph, linegraph

#################################################

""" Overlapper submodule for HTStream charts and graphs """

#################################################

class Overlapper():

	########################
	# Info about App
	def __init__(self):
		self.info = "Attempts to overlap paired end reads to produce the original fragment, trims adapters, and can correct sequencing errors."
		self.type = "read_reducer"


	########################
	# Table Function
	def table(self, json, se_total_gain, index):

		# straight forward table construction.
		headers = OrderedDict()

		headers["Ov_%_Overlapped" + index] = {'title': "% Overlapped", 
									  'namespace': "% Overlapped",
									  'description': 'Percentage of Reads with Overlap.',
									  'suffix': '%',
									  'format': '{:,.2f}',
									  'scale': 'Greens'}

		# If SE were gained, add col
		if se_total_gain != 0:
			headers["Ov_SE_gain" + index] = {'title': "% SE Gained", 'namespace': "% SE Gained",'description': 'Percentage Increase of Single End Reads', 'format': '{:,.2f}', 
											 'suffix': '%', 'scale': 'Blues' }

		headers["Ov_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)


	########################
	# Bargraph Function
	def bargraph(self, json, inserts):

		# configuration dictionary for bar graph
		config = {'title': "HTStream: Overlap Composition Bargraph",
				  'id': "htstream_overlapper_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False}

		# Header 
		html = "<h4> Overlapper: Overlap Composition </h4>\n"
		html += "<p>Plots the quantities of insert types for each sample.</p>"


		# if no overlaps at all are present, return nothing
		if inserts == 0:
			html += '<div class="alert alert-info"> <strong>Notice:</strong> No overlaps present in samples. </div>'	
			return html

		# If too many samples, dont plot
		if len(json.keys()) > 150:
			html += '<div class="alert alert-warning"> <strong>Warning:</strong> Too many samples for bargraph. </div>'	
			return html

		# bargraph dictionary. Exact use of example in MultiQC docs.
		categories  = OrderedDict()

		# Create blocks for bargrapph
		categories['Ov_Sins'] = {'name': 'Short Inserts', 'color': '#779BCC'}
		categories['Ov_Mins'] = {'name': 'Medium Inserts', 'color': '#C3C3C3'}
		categories['Ov_Lins'] = {'name': 'Long Inserts', 'color': '#D1ADC3'}

		# create plot
		html += bargraph.plot(json, categories, config)
		
		return html


	########################
	# Linegraph Function
	def linegraph(self, json, index):

		# config dictionary for "density" plots. Its a work in progress. 
		config = {'id': "htstream_overlapper_linegraph_" + index,
				  'title': "HTStream: Overlapped Lengths",
				  'ylab': "Counts", "xlab": "Overlap Lengths"}

		# initialize data structures
		multi_line = {}

		for key in json.keys():

			# creates empty dictionary to hold data for line graph. 
			multi_line[key] = {}

			# iterates over ever value in histogram and adds it to line graph
			for item in json[key]["Ov_Histogram"]:

				multi_line[key][item[0]] = item[1]

		html = "<h4> Overlapper: Overlapped Lengths </h4>\n"
		html += "<p>Plots the lengths of paired end read overlaps.</p>"
		html += linegraph.plot(multi_line, config)

		return html

	########################
	# Function for parsing histogram in files
	def parse_histogram_stats(self, hist):

		hist_stats = {"Max": 0,
					  "Median": 0}

		median_list = []

		# Find median and max of list in files
		for item in hist:
			median_list.append(item[0])

			if hist_stats["Max"] < item[1]:
				hist_stats["Max"] = item[0]


		hist_stats["Median"] = statistics.median(median_list)

		return hist_stats


	########################
	# Main Function
	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		# accumulator for inserts, used to prevent empty bar graph
		inserts = 0
		se_total_gain = 0

		for key in json.keys():

			sins = json[key]["Fragment"]["inserts"]["short"]
			mins = json[key]["Fragment"]["inserts"]["medium"]
			lins = json[key]["Fragment"]["inserts"]["long"]

			# Total overlap types
			overlapped_sum = (sins + mins + lins)

			# the INFAMOUS percent overlapped 
			perc_overlapped = ((sins + mins) / json[key]["Paired_end"]["in"]) * 100
			perc_pe_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"]) * 100
			
			# if no single end, prevent zero division
			if json[key]["Single_end"]["in"] == 0:
				perc_se_gain = 0

			else:
				perc_se_gain = ((json[key]["Single_end"]["out"] - json[key]["Single_end"]["in"]) / json[key]["Single_end"]["in"]) * 100

			# total SE gain
			se_total_gain += perc_se_gain

			parsed_hist_stats = self.parse_histogram_stats(json[key]["Fragment"]["overlap_histogram"])

			# Overview stats
			overview_dict[key] = {
								  "PE_Output_Reads": json[key]["Paired_end"]["out"],
								  "SE_Output_Reads": json[key]["Single_end"]["out"],
								  "Overlap_Length_Max": parsed_hist_stats["Max"],
								  "Overlap_Length_Med": parsed_hist_stats["Median"],
								  "Sin": sins / json[key]["Fragment"]["in"],
								  "Min": mins / json[key]["Fragment"]["in"],
								  "Lin": lins / json[key]["Fragment"]["in"]
								 }

			# sample instance in dictionary
			stats_json[key] = {
							   "Ov_SE_gain" + index: perc_se_gain,
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
		section = {"Table": self.table(stats_json, se_total_gain, index),
				   "Overlap Composition": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json, index),
				   "Overview": overview_dict}
				   
		return section