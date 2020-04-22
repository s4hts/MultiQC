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
		status_dict = {}

		# header read type
		read_header = " ".join(read.split("_")[1:3])

		# section header
		html = '<h4> Base by Cycle: ' + read_header + '</h4>'

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
			sample_status = None
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
				sample_status = 'FAIL'
			elif sample_max >= 40:
				sample_status = 'QUESTIONABLE'
			else:
				sample_status = 'PASS'

			# adds color to sample in color dictionary
			status_dict[key] = sample_status

			# this config file is for the individual line of the multiline graph
			config["data_labels"].append({'name': key,'ylab': 'Percentage', 
										  'xlab': 'Cycle', 'yCeiling': 100, 'categories': True, 
										  'smooth_points_sumcounts': False})

			# append base by cycle to data for this to data list
			data_list.append(data)

		# this adds the html output of sample status. This function colors samples
		html += htstream_utils.sample_status(status_dict)

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

		btn_id = "-".join(read.split("_")[:3]).lower()

		line_data = {}
		status_dict = {}
		first = True
		button_list = []


		for key in json.keys():

			# create dictionary for line graph. Again, format is {x: y}
			line_data[key] = {}

			# creates unique heatmap id that can be queired later by js.
			heat_pconfig["id"] = "htstream_" + btn_id + "_" + key + "_heatmap"

			# creates x and y axis labels for heatmap (categorical)
			x_lab = json[key][read]["col_names"]
			y_lab = json[key][read]["row_names"][::-1] # reverse orientation makes it easier to cycle through

			data = []

			# create variables for range functions in loops. Represents shape of data
			quality_scores = json[key][read]["shape"][0]
			cycles = json[key][read]["shape"][-1]


			if read == "Single End Quality by Cycle":
				input_reads = json[key]["St_SE_in"]
			else:
				input_reads = json[key]["St_PE_in"]

			# temp total list 
			total = []
			
			# iterates through positions, creates a list of the sum of scores at each position to be used
			#	to calculated frequency for heatmap. Also, calculates avg. Q score for linegraph.
			#	This chunk of code is very ugly, but is a necessary evil. 

			num_above_q30 = 0

			for pos in range(cycles):
				temp = [ score_list[pos] for score_list in json[key][read]["data"] ]
				temp_sum = sum(temp)
				total.append(temp_sum)

				# multiples count at poistion by Q Score.
				total_score = sum([(int(p) * int(s)) for p, s in zip(temp, y_lab[::-1])])

				# divides sum of total score by the number of cycles for avg fragments
				line_data[key][pos] = total_score / input_reads 

				if line_data[key][pos] > 30:
					num_above_q30 += 1


			# check to see what percent of bases have a mean Q score of at least 30
			q30_gate = (num_above_q30 / cycles) 

			if q30_gate < 0.6:
				status_dict[key] = "FAIL"

			elif q30_gate < 0.8:
				status_dict[key] = "QUESTIONABLE"

			else:
				status_dict[key] = 'PASS'


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
				heatmap_html = heatmap.plot(data, x_lab, y_lab, heat_pconfig)

			else:
				active = "" # button is default off 
				hidediv = 'style="display: none;"' # div is hidden
				heatmap.plot(data, x_lab, y_lab, heat_pconfig)


			# html div attributes and text
			name = key
			pid = "htstream_" + btn_id + "_" + key + "_btn"

			button_list.append('<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this)" id="{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name))

	
		status_div = htstream_utils.sample_status(status_dict)

		line_plot = linegraph.plot(line_data, line_config)

		html = htstream_utils.qual_by_cycle_html(read, status_div, line_plot, btn_id, button_list, heatmap_html)

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
			histograms = ["St_R1_histogram", "St_R2_histogram"]
			reads = ["St_R1_Length", "St_R2_Length"]
		else:
			histograms = ["St_R1_histogram", "St_R2_histogram", "St_SE_histogram"]
			reads = ["St_R1_Length", "St_R2_Length", "St_SE_Length"]
			config["data_labels"].append({'name': "SE histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'})


		# initlialize data structures and important variables
		data_list = []
		invariant_dict = {}
		html = ""
		reads = []

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

					# format read name for dictionary
					read_name = "St_" + histograms[i].split("_")[1] + "_Length"

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
						total = json[key]["St_SE_in"]
					else:
						total = json[key]["St_PE_in"]
					
					# populate smaple dictionary with read length and its frequency
					for item in json[key][histograms[i]]:
						data[key][item[0]] = item[1] / total


			# if samples are in read data, append to dataset list for multiple graphs
			if len(data.keys()) != 0:
				data_list.append(data)


		# if samples with uniform read length are present
		if len(invariant_dict.keys()) != 0:

			# notice
			notice = 'Samples with uniform read lengths identified (displayed below). <br />'

			# table
			headers = OrderedDict()
			table_config = {'table_title': "Length of Uniform Reads"}

			# iterates samples and through possible reads, checks for presenc in dictionary
			for key  in invariant_dict.keys():

				for read_type in reads:

					#	assigns NA value if not present
					if invariant_dict[key].get(read_type, "ERROR") == 'ERROR':
						invariant_dict[key][read_type] = "NA"

			# instantiates table columns
			headers["St_R1_Length"] = {'title': "Read 1 Length", 'namespace': "Read 1 Length", 'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Greens' }
			headers["St_R2_Length"] = {'title': "Read 2 Length", 'namespace': "Read 2 Length", 'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Oranges' }
			headers["St_SE_Length"] = {'title': "Single End Length", 'namespace': "Single End Length", 'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Blues' }
			
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
				stats_json[key]["St_SE_histogram"] = json[key][-1]["Single_end"]["readlength_histogram"]
				stats_json[key]["St_Single_End_Base_by_Cycle"] = json[key][-1]["Single_end"]["base_by_cycle"]
				stats_json[key]["St_Single_End_Quality_by_Cycle"] = json[key][-1]["Single_end"]["qualities_by_cycle"]
				stats_json[key]["St_SE_in"] = json[key][-1]["Single_end"]["in"]
							   
				SE_presence = True

			except:
				SE_presence = False

			# sample instance in ordered dict
			stats_json[key]["St_R1_histogram"] = json[key][-1]["Paired_end"]["Read1"]["readlength_histogram"]
			stats_json[key]["St_R2_histogram"] = json[key][-1]["Paired_end"]["Read2"]["readlength_histogram"]
			stats_json[key]["St_Read_1_Base_by_Cycle"] = json[key][-1]["Paired_end"]["Read1"]["base_by_cycle"]
			stats_json[key]["St_Read_2_Base_by_Cycle"] = json[key][-1]["Paired_end"]["Read2"]["base_by_cycle"]
			stats_json[key]["St_Read_1_Quality_by_Cycle"] = json[key][-1]["Paired_end"]["Read1"]["qualities_by_cycle"]
			stats_json[key]["St_Read_2_Quality_by_Cycle"] =json[key][-1]["Paired_end"]["Read2"]["qualities_by_cycle"]
			stats_json[key]["St_PE_in"] = json[key][-1]["Paired_end"]["in"]


		# output dictionary, keys are section, value is function called for figure generation
		section = {
				   "Density Plots": self.linegraph(stats_json, SE_presence), 
				   "Base by Cycle (Read 1)": self.base_by_cycle(stats_json, "St_Read_1_Base_by_Cycle"),
				   "Quality by Cycle (Read 1)": self.quality_by_cycle(stats_json, "St_Read_1_Quality_by_Cycle"),
				   "Base by Cycle (Read 2)": self.base_by_cycle(stats_json, "St_Read_2_Base_by_Cycle"),
				   "Quality by Cycle (Read 2)": self.quality_by_cycle(stats_json, "St_Read_2_Quality_by_Cycle")
				   }

		# only executres if single read data is detected
		if SE_presence == True:
			section["Base by Cycle (Single End)"] = self.base_by_cycle(stats_json, "St_Single_End_Base_by_Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(stats_json, "St_Single_End_Quality_by_Cycle")

	
		return section 

