from collections import OrderedDict
import logging, statistics, math
from random import random

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, linegraph, heatmap

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def table(self, json, index):

		# striaght forward table function, right from MultiQC documentation
		headers = OrderedDict()

		headers["St_Fragments_in" + index] = {'title': "Input Reads", 'namespace': "Input Reads", 'description': 'Number of reads', 'format': '{:,.0f}', 'scale': "Greens"}
		headers["St_GC_Content" + index] = {'title': "GC Content", 'namespace': "GC Content", 'description': 'Percentage of bps that are G or C', 
									'format': '{:,.2f}', 'suffix': '%', 'scale': 'RdPu'}
		headers["St_N_Content" + index] = {'title': "N Content", 'namespace': "N Content", 'description': 'Percentage of bps that are N',
								   'format': '{:,.2f}', 'suffix': '%','scale': 'Blues'}
		headers["St_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}


		return table.plot(json, headers)

	def base_by_cycle(self, json, read):

		title_read = " ".join(read.split("_")[1:3])

		# config dictionary for line graph
		config = {'title': "HTStream: Base by Cycle (" + title_read + ")",
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



	def quality_by_cycle(self, json, read, index):

		# Here is the most complicated figure implementation in this module.
		#	The issues here are that MultiQC had limited options for displaying 
		#	multiple figures if its a heatmap. Also, it doesnt allow you to switch
		#	back and forth between figure typs. There are workarounds, however, using
		#	javascript and some clever organizations :).

		title_read = " ".join(read.split("_")[1:3])

		# config dictionary for mean Q score line graph
		line_config = {
				  'smooth_points_sumcounts': False,
				  'categories': True,
				  'title': "HTStream: Mean Quality by Cycle (" + title_read + ")",
				  'xlab': "Cycle",
				  'ylab': "Mean Q Score",
				  }

		# config dictionary for heatmaps
		heat_pconfig = {'id' : "",
				   'title': "HTStream: Quality by Cycle (" + title_read + ")",
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
		unique_id = str(random() % 1000)[2:]
		line_data = {}
		status_dict = {}
		first = True
		button_list = []


		for key in json.keys():

			# create dictionary for line graph. Again, format is {x: y}
			line_data[key] = {}

			# creates unique heatmap id that can be queired later by js.
			heat_pconfig["id"] = "htstream_" + btn_id + "_" + key + "_" + unique_id + "_heatmap_" + index

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

			num_above_q30 = 0

			for pos in range(cycles):
				temp = [ score_list[pos] for score_list in json[key][read]["data"] ]
				temp_sum = sum(temp)
				total.append(temp_sum)

				# multiples count at poistion by Q Score.
				total_score = sum([(int(p) * int(s)) for p, s in zip(temp, y_lab[::-1])])

				# divides sum of total score by the number of cycles for avg fragments
				line_data[key][pos] = total_score / temp_sum # total reads

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
				first = False # shuts off first gat
				heatmap_html = heatmap.plot(data, x_lab, y_lab, heat_pconfig)

			else:
				active = "" # button is default off 
				heatmap.plot(data, x_lab, y_lab, heat_pconfig)



			# html div attributes and text
			name = key
			pid = "htstream_" + btn_id + "_" + key + "_" + unique_id + "_btn"

			button_list.append('<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this, {i})" id="{pid}">{n}</button>\n'.format(a=active, i=index, pid=pid, n=name))

	
		status_div = htstream_utils.sample_status(status_dict)

		line_plot = linegraph.plot(line_data, line_config)

		html = htstream_utils.qual_by_cycle_html(read, status_div, line_plot, unique_id, button_list, heatmap_html, index)

		return html



	def histogram(self, json, read):

		read_keys = {"St_PE_histogram": "PE",
					 "St_SE_histogram": "SE"}


		read_code = read_keys[read]

		data = {}
		invariant_dict = {}
		button_list = []
		first = True
		notice_html = ""
		unique_id = str(random() % 1000)[2:]

		# iterates over all samples in input dictionary
		for key in json.keys():

			# if read length histogram has one value (ie. all samples have a uniform length),
			#	this data is added to a secondary table as to avoid ugly line graphs.
			if len(json[key][read][0]) == 1 and read_keys[read] != "PE":

				# format read name for dictionary
				read_length_col = "St_" + read_code  + "_Length"
				read_count_col = "St_" + read_code  + "_Reads"

					# try appending dictionary, if key doesn't exist, create the instance.
				try:
					invariant_dict[key][read_length_col] = json[key][read][0][0]
					invariant_dict[key][read_count_col] = json[key][read][0][1]

				except:
					invariant_dict[key] = {}
					invariant_dict[key][read_length_col] = json[key][read][0][0]
					invariant_dict[key][read_count_col] = json[key][read][0][1]

			# executes of more than one data points are found.
			else:
				dict_key = key + "_" + unique_id
				data[dict_key] = []

				for x in range(len(json[key][read])):

					max_reads = max([item[0] for item in json[key][read][x]]) + 1

					current = 10
					bins = []
					values = []

					while current < max_reads:
						bins.append(current)
						values.append(1) # pseudo count
						current += 1

					# populate smaple dictionary with read length and its frequency
					for item in json[key][read][x]:

						for x in range(len(bins) - 1, -1, -1):

							if item[0] == bins[x]:
								values[x] += item[1]
								break 

					data[dict_key].append({"bins": bins, "vals": values})

				if read_keys[read] == "SE":
					data[dict_key] = data[dict_key][-1]

				if first == True:
					active = "active" # button is default active
					first = False # shuts off first gat

				else:
					active = "" # button is default off 

				# # html div attributes and text
				pid  = "htstream_stats_" + read + "_" + key + "_" + unique_id + "_btn"
				read_id = read + "_" + unique_id
				button_list.append('<button class="btn btn-default btn-sm hist_btn {a}" onclick="htstream_histogram(\'{r}\', \'{s}_{u}\')" id="{p}">{s}</button>\n'.format(a=active, r=read_id, u=unique_id, s=key, p=pid))


		# if samples with uniform read length are present
		if len(invariant_dict.keys()) != 0:

			# notice
			notice = 'Samples with uniform read lengths identified (displayed below). <br />'

			# table
			headers = OrderedDict()
			table_config = {'table_title': "Length of Uniform Reads"}


			# instantiates table columns
			headers["St_{r}_Length".format(r = read_code)] = {'title': "Read Length ({r})".format(r = read_code), 
									   						  'namespace': "Read Length ({r})".format(r = read_code), 
									  						  'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Greens' }
			headers["St_{r}_Reads".format(r = read_code)] = {'title': "Read Count ({r})".format(r = read_code), 
															 'namespace': "Read Count ({r})".format(r = read_code),
															 'description': 'Length of Read Type', 'format': '{:,.0f}', 'scale': 'Purples' }
			
			# add to output html
			notice_html += '<div class="alert alert-info">{n}</div>'.format(n = notice)	
			notice_html += table.plot(invariant_dict, headers, table_config)	

		html = htstream_utils.stats_histogram_html(read, data, unique_id, button_list, notice_html)

		return html



	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_stats = {} 

		for key in json.keys():

			#
			# STATS FOR TABLE 
			#
			gc_count = (json[key]["Fragment"]["base_composition"]["G"] + json[key]["Fragment"]["base_composition"]["C"])
			gc_content = ( gc_count / json[key]["Fragment"]["basepairs_out"] ) * 100 
			n_content = ( json[key]["Fragment"]["base_composition"]["N"] / json[key]["Fragment"]["basepairs_out"] ) * 100 

			stats_json[key] = {"St_Fragments_in" + index: json[key]["Fragment"]["in"],
							   "St_GC_Content" + index: gc_content,
						       "St_N_Content" + index: n_content,
						       "St_Notes" + index: json[key]["Program_details"]["options"]["notes"]}


			overview_stats[key] = {"Output_Reads": json[key]["Fragment"]["out"],
								   "Output_Bp": json[key]["Fragment"]["basepairs_out"],
								   "gc_content": gc_content,
								   "n_content": n_content,
								   "total_Q30": 0,
								   "Read_Breakdown": {}}
			#
			# SINGLE END STATS
			#
			# only succeeds if json file contains single end information data in the last instance of hts_Stats,
			#	opens gate for future processing of single end read stats.
			try:
				stats_json[key]["St_SE_histogram"] = [json[key]["Single_end"]["readlength_histogram"]]
				stats_json[key]["St_Single_End_Base_by_Cycle"] = json[key]["Single_end"]["base_by_cycle"]
				stats_json[key]["St_Single_End_Quality_by_Cycle"] = json[key]["Single_end"]["qualities_by_cycle"]
				stats_json[key]["St_SE_in"] = json[key]["Single_end"]["in"]
				
				overview_stats[key]["total_Q30"] += json[key]["Single_end"]["total_Q30_basepairs"]
				overview_stats[key]["Read_Breakdown"]["Single_end"] = json[key]["Single_end"]["in"]

				SE_presence = True

			except:
				SE_presence = False


			#
			# PAIRED END STATS
			#
			try:
				# sample instance in ordered dict
				stats_json[key]["St_PE_histogram"] = [json[key]["Paired_end"]["Read1"]["readlength_histogram"],
													  json[key]["Paired_end"]["Read2"]["readlength_histogram"]]
				stats_json[key]["St_Read_1_Base_by_Cycle"] = json[key]["Paired_end"]["Read1"]["base_by_cycle"]
				stats_json[key]["St_Read_2_Base_by_Cycle"] = json[key]["Paired_end"]["Read2"]["base_by_cycle"]
				stats_json[key]["St_Read_1_Quality_by_Cycle"] = json[key]["Paired_end"]["Read1"]["qualities_by_cycle"]
				stats_json[key]["St_Read_2_Quality_by_Cycle"] = json[key]["Paired_end"]["Read2"]["qualities_by_cycle"]
				stats_json[key]["St_PE_in"] = json[key]["Paired_end"]["in"]


				overview_stats[key]["total_Q30"] += json[key]["Paired_end"]["Read1"]["total_Q30_basepairs"] + json[key]["Paired_end"]["Read2"]["total_Q30_basepairs"]
				overview_stats[key]["Read_Breakdown"]["Paired_end"] = json[key]["Paired_end"]["in"]				   

				PE_presence = True

			except:
				PE_presence = False 




		# output dictionary, keys are section, value is function called for figure generation
		section = {"Table": self.table(stats_json, index),
				   "Overview": overview_stats}

		if PE_presence == True:
			section["Read Length Histogram (Paried End)"] = self.histogram(stats_json, "St_PE_histogram")
			section["Base by Cycle (Read 1)"] = self.base_by_cycle(stats_json, "St_Read_1_Base_by_Cycle")
			section["Base by Cycle (Read 2)"] = self.base_by_cycle(stats_json, "St_Read_2_Base_by_Cycle")
			section["Quality by Cycle (Read 1)"] = self.quality_by_cycle(stats_json, "St_Read_1_Quality_by_Cycle", index)
			section["Quality by Cycle (Read 2)"] = self.quality_by_cycle(stats_json, "St_Read_2_Quality_by_Cycle", index)


		#only executres if single read data is detected
		if SE_presence == True:
			section["Read Length Histogram (Single End)"] = self.histogram(stats_json, "St_SE_histogram")
			section["Base by Cycle (Single End)"] = self.base_by_cycle(stats_json, "St_Single_End_Base_by_Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(stats_json, "St_Single_End_Quality_by_Cycle", index)

	
		return section 

