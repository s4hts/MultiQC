from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, linegraph

#################################################

""" SuperDeduper submodule for HTStream charts and graphs """

#################################################

class SuperDeduper():


	def table(self, json):

		# striaght forward table function, right from MultiQC documentation
		headers = OrderedDict()

		headers["Sd_PE_loss"] = {'title': "% PE Lost", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.2f}', 
								 'suffix': '%', 'scale': 'Greens' }
		headers["Sd_SE_in"] = {'title': "SE in", 'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Sd_SE_out"] = {'title': "SE out", 'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Sd_%_Duplicates"] = {'title': "% Duplicates", 
								   'namespace': "% Duplicates", 
								   'description': 'Percentage of Duplicate Reads (SE and PE)',
								   'suffix': '%',
								   'format': '{:,.2f}',
								   'scale': 'Oranges'
								  }
		headers["Sd_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}


		return table.plot(json, headers)



	def linegraph(self, json):

		# plot configurations, list of options in MultiQC docs
		config = {'title': "HTStream: Duplicate Saturation",
				  'xlab': "Total Reads", 'ylab': "Total Reads - Duplicates",
				  'extra_series': []}

		# initialize data structures and variabe;s 
		data = {}
		invariant_saturation_dict = {}
		html = ""

		for key in json.keys():

			# if duplicate saturation histogram has data point, it is added to 'invariant_saturation_dict', where 
			# 	it will be represented as table instead of a hideous graph.
			if len(json[key]["Sd_Saturation"]) == 1:
				invariant_saturation_dict[key] = {"Sd_Total_Reads": json[key]["Sd_Saturation"][0][0], "Sd_Duplicates": json[key]["Sd_Saturation"][0][1]}


			# if more than one data point is identified (low bar, I know), it will be added to the graph's data
			#	dictionary. Data points represented as dictionary: {x: y}.
			else:
				data[key] = {}

				for item in json[key]["Sd_Saturation"]:

					data[key][item[0]] = item[0] - item[1] 


		# checks for any invariant samples and creates an alert div and table to  hold the data.
		if len(invariant_saturation_dict.keys()) != 0:
			# notice
			notice = 'Samples with uniform duplication numbers identified (displayed below). <br />'

			# table
			headers = OrderedDict()
			headers["Sd_Total_Reads"] = {'title': "Total Reads", 'namespace': "Total Reads", 'description': 'Number of Total Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
			headers["Sd_Duplicates"] = {'title': "Total Reads - Duplicates", 'namespace': "Duplicates", 'description': 'Number of Duplicates', 'format': '{:,.0f}', 'scale': 'RdPu'}
			
			# add to output html
			html += '<div class="alert alert-info">{n}</div>'.format(n = notice)	
			html += table.plot(invariant_saturation_dict, headers)
			

		# creates line graph only if samples with more than one data point are presents.
		if data != {}:
			html += linegraph.plot(data, config)

		return html


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			# number of duplicates reletive to input reads 
			perc_duplicates = (json[key]["Fragment"]["duplicate"] / json[key]["Fragment"]["in"]) * 100
			perc_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"])  * 100

			# sample instance in ordered dict
			stats_json[key] = {
			 				   "Sd_PE_loss": perc_loss,
							   "Sd_SE_in": json[key]["Single_end"]["in"],
							   "Sd_SE_out": json[key]["Single_end"]["out"],
							   "Sd_%_Duplicates": perc_duplicates,
							   "Sd_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Sd_Duplicates": json[key]["Fragment"]["duplicate"],
							   "Sd_Saturation": json[key]["Fragment"]["duplicate_saturation"]
						 	  }

		# output dictionary, keys are section, value is function called for figure generation
		section = {
				   "Table": self.table(stats_json),
				   "Duplicate Saturation": self.linegraph(stats_json)
				   }

		return section