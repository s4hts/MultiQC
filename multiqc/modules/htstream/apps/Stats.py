from collections import OrderedDict
import logging, statistics

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, linegraph, heatmap

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def base_by_cycle(self, json, read):

		# config dictionary for line graph
		config = {'title': "HTStream: Base by Cycle",
				  'data_labels': [],
				  'smooth_points_sumcounts': False,
				  'yCeiling': 100,
				  'categories': True,
				  'colors': {
				  			 "A": "#B62612",
				  			 "C": "#82A7E0",
				  			 "G": "#0B8E0B",
				  			 "T": "#DE7D00",
				  			 "N": "black"
				  			},
				  'yPlotBands': [
								{'from': 0, 'to': 40, 'color': '#c3e6c3'},
								{'from': 40, 'to': 60, 'color': '#e6dcc3'},
								{'from': 60, 'to': 100, 'color': '#e6c3c3'},
								]
				  }

		# initalize data structures and important variables
		data_list = []
		color_dict = {}
		html = ""


		for key in json.keys():

			# initializes dat dict. Each key is a line in the graph
			data = {"A": {},
					"C": {},
					"G": {},
					"T": {},
					"N": {}}

			# lists to iterate through
			bases = json[key][read]["data"]
			positions = json[key][read]["col_names"]

			# vairables containing max percentage reached by any nucleotide in the sample
			#	This data is stored so it can be correctly marked in the sample check div.
			sample_color = None
			sample_max = 0

			# iterates through every position
			for i in range(len(positions)):

				# total base calls at that position, for some reason, this is not equal to the number 
				#	of input reads? Potential error. 
				total = bases[0][i] + bases[1][i] + bases[2][i] + bases[3][i] + bases[4][i]

				# list of values for heatmap, just cleaner to put them in a list before hand
				y_value_list = [(bases[0][i] / total) * 100, (bases[1][i] / total) * 100, 
								(bases[2][i] / total) * 100, (bases[3][i] / total) * 100,
								(bases[4][i] / total) * 100]

				# take max for position and compare it to max for entire sample
				sample_max = max([sample_max, max(y_value_list)])

				# add data to dictionary for each base
				data["A"][i] = y_value_list[0]
				data["C"][i] = y_value_list[1]
				data["G"][i] = y_value_list[2]
				data["T"][i] = y_value_list[3]
				data["N"][i] = y_value_list[4]


			# selects color to mark sample if a read has a region of low complextity
			if sample_max >= 60:
				sample_color = '#e6c3c3'
			elif sample_max >= 40:
				sample_color = '#e6dcc3'
			else:
				sample_color = '#c3e6c3'

			# adds color to sample in color dictionary
			color_dict[key] = sample_color

			# this config file is for the individual line of the multiline graph
			config["data_labels"].append({'name': key,'ylab': 'Percentage', 
										  'xlab': 'Cycle', 'yCeiling': 100, 'categories': True, 
										  'smooth_points_sumcounts': False})

			# append base by cycle to data for this to data list
			data_list.append(data)

		# this adds the html output of sample status. This function colors samples
		html += htstream_utils.sample_status(color_dict)

		# add line graphs
		html += linegraph.plot(data_list, config)

		return html


	def quality_by_cycle(self, json, read):

		# Here is the most complicated figure implementation in this whole module.
		#	The issues here are that MultiQC had limited options for displaying 
		#	multiple figures if its a heatmap. Also, it doesnt allow you to switch
		#	back and forth between figure typs. There are workarounds, however, using
		#	javascript and some clever organizations of javascript.


		# config dictionary for mean Q score line graph
		line_config = {
				  'smooth_points_sumcounts': False,
				  'categories': True,
				  'title': "HTStream: Quality by Cycle"
				  	 }

		# config dictionary for heatmaps
		heat_pconfig = {'id' : "",
				   'title': "HTStream: Quality by Cycle", 
				   'yTitle': 'Q Score',
				   'xTitle': 'Cycle',
				   'square' : False,
				   'datalabels': False,
				   'max': 1.0, 
				   'colstops': [
					        [0, '#FFFFFF'],
					        [0.3, '#1DC802'],
					        [0.6, '#F3F943'],
					        [1, '#E70808']
					           ]
    			  }

    	# id of switch buttun, named after read type.
		btn_id = "-".join(read.split(" ")[:2]).lower()

		# In order to be able to switch back and forth between different graph types, we need to add MultiQC's button divs that 
		#	have an onclick attribute for javascript that can switch between graphs.
		wrapper_html = '<div class="btn-group hc_switch_group">\n'
		wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this)" id="htstream_qbc_line_{r}_btn">Linegraph</button>\n'.format(r=btn_id)
		wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this)" id="htstream_qbc_heat_{r}_btn">Heatmaps</button></div>\n'.format(r=btn_id)
		wrapper_html += "<hr>"


		# The heatmaps of this section occur on a per sample basis, meaning we need another subset of buttons to switch between the samples
		html = '<div class="btn-group hc_switch_group">\n'


		# initiate important variables and data structures
		line_data = {}
		first = True
		pid = ""


		for key in json.keys():

			# create dictionary for line graph. Again, format is {x: y}
			line_data[key] = {}

			# creates unique heatmap id that can be queired later by js.
			heat_pconfig["id"] = "htstream_" + key + "_heatmap"

			# creates x and y axis labels for heatmap (categorical)
			x_lab = json[key][read]["col_names"]
			y_lab = json[key][read]["row_names"][::-1] # reverse orientation makes it easier to cycle through
			data = []

			# create variables for range functions in loops. Represents shape of data
			quality_scores = json[key][read]["shape"][0]
			cycles = json[key][read]["shape"][-1]

			# temp total list 
			total = []
			
			# iterates through positions, creates a list of the sum of scores at each position to be used
			#	to calculated frequency for heatmap. Also, calculates avg. Q score for linegraph.
			#	This chunk of code is very ugly, but is a necessary evil. 
			for pos in range(cycles):
				temp = [ score_list[pos] for score_list in json[key][read]["data"] ]
				total.append(sum(temp))
				line_data[key][pos] = statistics.mean(temp)

			# populates data dictionaries for heatmap
			for score in range(quality_scores - 1, -1, -1):

				# create empty list for data. The format is a little strange, each list represents a position 
				#	the value inside of it is the score at that position divided by the total score for that position
				#	giving a frequency.
				data.append([])

				for pos in range(cycles):
					data[-1].append(json[key][read]["data"][score][pos] / total[pos])


			# if this is the first sample process, lucky them, they get to be shown first and marked as active.
			#	This step is necessary otherwise, the plot div is not initialized. The additional calls to the 
			#	heatmap function are simply to add the data to the internal jsons used by MultiQC.
			if first == True:
				active = "active" # button is default active
				hidediv = "" # shows div
				first = False # shuts off first gat
				plot_html = heatmap.plot(data, x_lab, y_lab, heat_pconfig)

			else:
				active = "" # button is default off 
				hidediv = 'style="display: none;"' # div is hidden
				heatmap.plot(data, x_lab, y_lab, heat_pconfig)


			# html div attributes and text
			name = key
			pid = "htstream_" + key + "_btn"

			# add html div for button. Occurs for every sample.
			html += '<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this)" id="{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name)


		# clean up the html and add spacing and heatmap div
		html += '</div>\n\n<br></br>\n\n'
		html += plot_html 

		# this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)

		# line graph div
		wrapper_html += '<div id="htstream_qbc_line_{r}">'.format(r=btn_id)
		wrapper_html += linegraph.plot(line_data, line_config) + "</div>"

		# heatmap div
		wrapper_html += '<div id="htstream_qbc_heat_{r}" style="display:none;">'.format(r=btn_id)
		wrapper_html += html + "</div>"

		html = wrapper_html 

		return html



	def linegraph(self, json, SE_presence):

		# config dict for line graphs
		config = {'title': "HTStream: Read Length Histogram",
				  'data_labels': [
								  {'name': "R1 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "R2 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'}
								  ]}


		# if single end read data is present, this data will be appended to the data set and config dicts.
		if SE_presence == False:
			histograms = ["R1 histogram", "R2 histogram"]
		else:
			histograms = ["R1 histogram", "R2 histogram", "SE histogram"]
			config["data_labels"].append({'name': "SE histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'})


		# initlialize data structures and important variables
		data_list = []
		invariant_dict = {}
		html = ""

		# iterates of all types of read data available (R1, R2, & SE (sometimes))
		for i in range(len(histograms)):

			# creates new data dictionary for every sample. Dictionary contains data points
			#	in following format: {x: y}.
			data = {}

			# iterates over all samples in input dictionary
			for key in json.keys():

				# if read length histogram has one value (ie. all samples have a uniform length),
				#	this data is added to a secondary table as to avoid ugly line graphs.
				if len(json[key][histograms[i]]) == 1:

					# format read name
					read_name = histograms[i].split(" ")[0] + " Length"

					# try appending dictionary, if key doesn't exist, create the instance.
					try:
						invariant_dict[key][read_name] = json[key][histograms[i]][0][0]

					except:
						invariant_dict[key] = {}
						invariant_dict[key][read_name] = json[key][histograms[i]][0][0]

				# executes of more than one data points are found.
				else:

					data[key] = {}

					# sums the total number of reads 
					if histograms[i] == "SE histogram":
						total = json[key]["SE in"]
					else:
						total = json[key]["PE in"]
					
					# populate smaple dictionary with read length and its frequency
					for item in json[key][histograms[i]]:
						data[key][item[0]] = item[1] / total


			# if samples are in read data, append to dataset list for multiple graphs
			if len(data.keys()) != 0:
				data_list.append(data)


		# if samples with uniform read length are present
		if len(invariant_dict.keys()) != 0:

			# notice
			notice = 'Samples with uniform read lengths identified. <br />'

			# table
			headers = OrderedDict()
			table_config = {'table_title': "Length of Uniform Reads"}

			reads = ["R1 Length", "R2 Length", "SE Length"]

			# iterates samples and through possible reads, checks for presenc in dictionary
			for key  in invariant_dict.keys():

				for read_type in reads:

					#	assigns NA value if not present
					if invariant_dict[key].get(read_type, "ERROR") == 'ERROR':
						invariant_dict[key][read_type] = "NA"

			# instantiates table columns
			headers["R1 Length"] = {'namespace': "Read 1 Length",'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Greens' }
			headers["R2 Length"] = {'namespace': "Read 2 Length",'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Oranges' }
			headers["SE Length"] = {'namespace': "Single End Length",'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Blues' }
			
			# add to output html
			html += '<div class="alert alert-info">{n}</div>'.format(n = notice)	
			html += table.plot(invariant_dict, headers, table_config)	
			

		# if data present for line graphs, make graphs
		if len(data_list) != 0:
			html += linegraph.plot(data_list, config)

		return html



	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {}

			# only succeeds if json file contains single end information data in the last instance of hts_Stats,
			#	opens gate for future processing of single end read stats.
			try:
				stats_json[key]["SE histogram"] = json[key][-1]["Single_end"]["readlength_histogram"]
				stats_json[key]["Single End Base by Cycle"] = json[key][-1]["Single_end"]["base_by_cycle"]
				stats_json[key]["Single End Quality by Cycle"] = json[key][-1]["Single_end"]["qualities_by_cycle"]
				stats_json[key]["SE in"] = json[key][-1]["Single_end"]["in"]
							   
				SE_presence = True

			except:
				SE_presence = False

			# sample instance in ordered dict
			stats_json[key]["R1 histogram"] = json[key][-1]["Paired_end"]["Read1"]["readlength_histogram"]
			stats_json[key]["R2 histogram"] = json[key][-1]["Paired_end"]["Read2"]["readlength_histogram"]
			stats_json[key]["Read 1 Base by Cycle"] = json[key][-1]["Paired_end"]["Read1"]["base_by_cycle"]
			stats_json[key]["Read 2 Base by Cycle"] = json[key][-1]["Paired_end"]["Read2"]["base_by_cycle"]
			stats_json[key]["Read 1 Quality by Cycle"] = json[key][-1]["Paired_end"]["Read1"]["qualities_by_cycle"]
			stats_json[key]["Read 2 Quality by Cycle"] =json[key][-1]["Paired_end"]["Read2"]["qualities_by_cycle"]
			stats_json[key]["PE in"] = json[key][-1]["Paired_end"]["in"]


		# output dictionary, keys are section, value is function called for figure generation
		section = {
				   "Density Plots": self.linegraph(stats_json, SE_presence), 
				   "Base by Cycle (Read 1)": self.base_by_cycle(stats_json, "Read 1 Base by Cycle"),
				   "Quality by Cycle (Read 1)": self.quality_by_cycle(stats_json, "Read 1 Quality by Cycle"),
				   "Base by Cycle (Read 2)": self.base_by_cycle(stats_json, "Read 2 Base by Cycle"),
				   "Quality by Cycle (Read 2)": self.quality_by_cycle(stats_json, "Read 2 Quality by Cycle")
				   }

		# only executres if single read data is detected
		if SE_presence == True:
			section["Base by Cycle (Single End)"] = self.base_by_cycle(stats_json, "Single End Base by Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(stats_json, "Single End Quality by Cycle")

	
		return section 

