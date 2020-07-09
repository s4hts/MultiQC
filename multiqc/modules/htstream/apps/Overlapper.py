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

	def table(self, json, se_total_gain, index):

		config = {'namespace': 'overlapper'}

		# straight forward table construction.
		headers = OrderedDict()

		headers["Ov_%_Overlapped" + index] = {'title': "% Overlapped", 
									  'namespace': "% Overlapped",
									  'description': 'Percentage of Reads with Overlap.',
									  'suffix': '%',
									  'format': '{:,.2f}',
									  'scale': 'Greens'}

		if se_total_gain != 0:
			headers["Ov_SE_gain" + index] = {'title': "% SE Gained", 'namespace': "% SE Gained",'description': 'Percentage Increase of Single End Reads', 'format': '{:,.2f}', 
											 'suffix': '%', 'scale': 'Blues' }

		headers["Ov_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, inserts):

		# configuration dictionary for bar graph
		config = {'title': "HTStream: Overlap Composition Bargraph",
				  'id': "htstream_overlapper_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False}

		html = "<h4> Overlap Composition </h4>\n"

		# if no overlaps at all are present, return nothing
		if inserts == 0:
			html += '<div class="alert alert-info"> <strong>Notice:</strong> No overlaps present in samples. </div>'	
			return html

		if len(json.keys()) > 150:
			html += '<div class="alert alert-warning"> <strong>Warning:</strong> Too many samples for bargraph. </div>'	
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

		html += bargraph.plot(json, categories, config)

		return html



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

		html = "<h4> Overlapped Lengths </h4>\n" +linegraph.plot(multi_line, config)

		return html


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
		se_total_gain = 0

		for key in json.keys():

			sins = json[key]["Fragment"]["inserts"]["short"]
			mins = json[key]["Fragment"]["inserts"]["medium"]
			lins = json[key]["Fragment"]["inserts"]["long"]

			overlapped_sum = (sins + mins + lins)

			# the INFAMOUS percent overlapped 
			perc_overlapped = ((sins + mins) / json[key]["Paired_end"]["in"]) * 100
			perc_pe_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"]) * 100
			
			if json[key]["Single_end"]["in"] == 0:
				perc_se_gain = 0

			else:
				perc_se_gain = ((json[key]["Single_end"]["out"] - json[key]["Single_end"]["in"]) / json[key]["Single_end"]["in"]) * 100

			se_total_gain += perc_se_gain

			parsed_hist_stats = self.parse_histogram_stats(json[key]["Fragment"]["overlap_histogram"])

			overview_dict[key] = {
								  "Output_Reads": json[key]["Fragment"]["out"],
								  "PE_reads_out": (json[key]["Paired_end"]["out"] / json[key]["Fragment"]["out"]) * 100,
								  "SE_reads_out": (json[key]["Single_end"]["out"] / json[key]["Fragment"]["out"]) * 100,
								  "Overlap_Length_Max": parsed_hist_stats["Max"],
								  "Overlap_Length_Med": parsed_hist_stats["Median"],
								  # "Fraction_PE_Lost": (json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Fragment"]["in"],
								  # "Fraction_Bp_Lost": (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) / json[key]["Fragment"]["basepairs_in"],
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