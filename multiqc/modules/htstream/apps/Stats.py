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
		# "St_PE_Fraction" + index
		headers["St_PE_Fraction" + index] = {'title': "% PE", 'namespace': "% PE", 'description': 'percentage of paried end reads', 'format': '{:,.0f}', 
											 'format': '{:,.2f}', 'suffix': '%', 'scale': "Greens"}
		headers["St_SE_Fraction" + index] = {'title': "% SE", 'namespace': "% SE", 'description': 'percentage of paried end reads', 'format': '{:,.0f}', 
											 'format': '{:,.2f}', 'suffix': '%', 'scale': "RdPu"}
		headers["St_R1_Q30" + index] = {'title': "% R1 Q30", 'namespace': "% R1 Q30", 'description': 'percentage of read 1 bps Q30 or greater', 
										'format': '{:,.2f}', 'suffix': '%', 'scale': 'Blues'}
		headers["St_R2_Q30" + index] = {'title': "% R2 Q30", 'namespace': "% R2 Q30", 'description': 'percentage of read 2 bps Q30 or greater', 
										'format': '{:,.2f}', 'suffix': '%','scale': 'Greens'}
		headers["St_SE_Q30" + index] = {'title': "% SE Q30", 'namespace': "% SE Q30", 'description': 'percentage of single end read bps Q30 or greater', 
										'format': '{:,.2f}', 'suffix': '%', 'scale': 'RdPu'}
		headers["St_GC_Content" + index] = {'title': "GC Content", 'namespace': "GC Content", 'description': 'Percentage of bps that are G or C', 
									'format': '{:,.2f}', 'suffix': '%', 'scale': 'Blues'}
		headers["St_N_Content" + index] = {'title': "N Content", 'namespace': "N Content", 'description': 'Percentage of bps that are N',
								   'format': '{:,.4f}', 'suffix': '%','scale': 'Green'}
		headers["St_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}


		return table.plot(json, headers)



	def base_by_cycle(self, json, read, index):

		title_read = " ".join(read.split("_")[1:3])

		# config dictionary for line graph
		config = {'id': "htstream_stats_base_" + read + "_" + index,
				  'title': "HTStream: Base by Cycle (" + title_read + ")",
				  'data_labels': [],
				  'smooth_points_sumcounts': False,
				  'ylab': "Percentage",
				  'xlab': "Cycle",
				  'yCeiling': 100,
				  'categories': True,
				  'tt_decimals': '{:,.2f}',
				  'tt_suffix': "%",     
				  'colors': {
				  			 "Base: A": "#B62612",
				  			 "Base: C": "#82A7E0",
				  			 "Base: G": "#0B8E0B",
				  			 "Base: T": "#DE7D00",
				  			 "Base: N": "black"
				  			},
				  'yPlotBands': [
								{'from': 0, 'to': 40, 'color': '#c3e6c3'},
								{'from': 40, 'to': 60, 'color': '#e6dcc3'},
								{'from': 60, 'to': 100, 'color': '#e6c3c3'},
								]

				  }


		# header read type
		read_header = " ".join(read.split("_")[1:3])

		if read_header == "Paired End":

			midpoint = 0
			uniform_pe = True

			for key in json.keys():

				temp = (json[key][read][0]["shape"][-1] * 2) 

				if midpoint == 0:
					midpoint = temp

				elif midpoint == temp:
					midpoint = midpoint

				else:
					uniform_pe = False
					break

			if uniform_pe == True:
				config['xPlotLines'] = [{'color': '#5D4B87', 
										 "width": 1.5, 
										 "value": (midpoint - 1) / 2, 
										 "dashStyle": 'shortdash',
										 "zIndex": 4}]

				
		# initalize data structures and important variables
		data_list = []
		status_dict = {}

		# section header
		html = '<h4> Base by Cycle: ' + read_header + '</h4>'

		for key in json.keys():

			if read_header == "Paired End":

				temp_data = [ json[key][read][0]["data"][x] + json[key][read][1]["data"][x] for x in range(5) ]
				temp_col_name = [ str(int(x) + int(json[key][read][0]["col_names"][-1]) ) for x in json[key][read][1]["col_names"]]

				json[key][read] = {
					 			   "data": temp_data,
					 			   "col_names": json[key][read][0]["col_names"] + temp_col_name
							 	  }

				del temp_data
				del temp_col_name

			# initializes dat dict. Each key is a line in the graph
			data = {"Base: A": {},
					"Base: C": {},
					"Base: G": {},
					"Base: T": {},
					"Base: N": {}}

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
				data["Base: A"][i] = y_value_list[0]
				data["Base: C"][i] = y_value_list[1]
				data["Base: G"][i] = y_value_list[2]
				data["Base: T"][i] = y_value_list[3]
				data["Base: N"][i] = y_value_list[4]


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
			config["data_labels"].append({'name': key, 'ylab': 'Percentage', 
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

		def temp_split(string):
			return string.split("hts")[0]

		title_read = " ".join(read.split("_")[1:3])

		# config dictionary for mean Q score line graph
		line_config = {
				  'id': "htstream_stats_qbc_" + read + "_" + index,
				  'smooth_points_sumcounts': False,
				  'categories': True,
				  'tt_decimals': '{:,.2f}',
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
					           ],
    			  }

		read_header = " ".join(read.split("_")[1:3])

		# check for uniform pe length
		if read_header == "Paired End":

			midpoint = 0
			uniform_pe = True

			for key in json.keys():

				temp = (json[key][read][0]["shape"][-1] * 2) 

				if midpoint == 0:
					midpoint = temp

				elif midpoint == temp:
					midpoint = midpoint

				else:
					uniform_pe = False
					break

			if uniform_pe == True:
				line_config['xPlotLines'] = [{'color': '#5D4B87', 
											 "width": 1.5, 
											 "value": (midpoint - 1) / 2, 
											 "dashStyle": 'shortdash',
											 "zIndex": 4}]


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

			if read_header == "Paired End":

				length = len(json[key][read][0]["data"])
				temp_data = [ json[key][read][0]["data"][x] + json[key][read][1]["data"][x] for x in range(length) ]
				temp_col_name = [ str(int(x) + int(json[key][read][0]["col_names"][-1]) ) for x in json[key][read][1]["col_names"]]

				json[key][read] = {
					 			   "data": temp_data,
					 			   "col_names": json[key][read][0]["col_names"] + temp_col_name,
					 			   "row_names": json[key][read][0]["row_names"],
					 			   "shape": [json[key][read][0]["shape"][0], json[key][read][0]["shape"][-1] + json[key][read][1]["shape"][-1]]
							 	  }

				del temp_data
				del temp_col_name

			# creates x and y axis labels for heatmap (categorical)
			x_lab = [ str(int(x) - 1) for x in json[key][read]["col_names"]]
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
			suffix = "_heatmap_" + index

			button_list.append('<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this, \'{s}\')" id="{pid}">{n}</button>\n'.format(a=active, s=suffix, pid=pid, n=name))


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
			
			dict_key = key + "_" + unique_id
			data[dict_key] = []

			for x in range(len(json[key][read])):

				max_reads = max([item[0] for item in json[key][read][x]]) + 1

				current = 10
				bins = []
				values = []

				while current < max_reads:
					values.append(0)
					bins.append(current)
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

		html = htstream_utils.stats_histogram_html(read, data, unique_id, button_list, notice_html)

		return html



	def execute(self, json, index):

		stats_json = OrderedDict()
		SE_json = OrderedDict()
		PE_json = OrderedDict()
		overview_stats = {} 

		for key in json.keys():

			#
			# STATS FOR TABLE 
			#
			total_frags = json[key]["Fragment"]["out"]
			total_bps = json[key]["Fragment"]["basepairs_out"]
			gc_count = (json[key]["Fragment"]["base_composition"]["G"] + json[key]["Fragment"]["base_composition"]["C"])
			gc_content = gc_count / total_bps
			n_content = json[key]["Fragment"]["base_composition"]["N"] / total_bps

			stats_json[key] = {"St_Fragments_in" + index: total_frags,
							   "St_GC_Content" + index: gc_content * 100 ,
						       "St_N_Content" + index: n_content * 100 ,
						       "St_Notes" + index: json[key]["Program_details"]["options"]["notes"]}


			overview_stats[key] = {
								   "GC_Content": gc_content,
								   "N_Content": n_content,
								   "Q30_Fraction": 0
								   }
			#
			# SINGLE END STATS
			#
			# only succeeds if json file contains single end information data in the last instance of hts_Stats,
			#	opens gate for future processing of single end read stats.
			try:
				SE_json[key] = {}
				SE_json[key]["St_SE_histogram"] = [json[key]["Single_end"]["readlength_histogram"]]
				SE_json[key]["St_Single_End_Base_by_Cycle"] = json[key]["Single_end"]["base_by_cycle"]
				SE_json[key]["St_Single_End_Quality_by_Cycle"] = json[key]["Single_end"]["qualities_by_cycle"]
				SE_json[key]["St_SE_in"] = json[key]["Single_end"]["in"]

				stats_json[key]["St_SE_Fraction" + index] = (json[key]["Single_end"]["out"] / total_frags) * 100
				stats_json[key]["St_SE_Q30" + index] = ( json[key]["Single_end"]["total_Q30_basepairs"] / json[key]["Single_end"]["basepairs_in"] ) * 100

				overview_stats[key]["Q30_Fraction"] += (json[key]["Single_end"]["total_Q30_basepairs"] / total_bps)
				overview_stats[key]["SE_Fraction"] = json[key]["Single_end"]["out"] / total_frags	
				overview_stats[key]["SE_Output_Reads"] = json[key]["Single_end"]["out"]
				overview_stats[key]["SE_Output_Bps"] = json[key]["Single_end"]["basepairs_out"] 
				

			except:
				overview_stats[key]["Q30_Fraction"] += 0
				overview_stats[key]["SE_Fraction"] = 0
				overview_stats[key]["SE_Output_Reads"] = 0 
				overview_stats[key]["SE_Output_Bps"] = 0


				del SE_json[key]



			#
			# PAIRED END STATS
			#
			try:
				PE_json[key] = {}	
				PE_json[key]["St_PE_histogram"] = [json[key]["Paired_end"]["Read1"]["readlength_histogram"],
													  json[key]["Paired_end"]["Read2"]["readlength_histogram"]]
				PE_json[key]["St_Paired_End_Base_by_Cycle"] = [json[key]["Paired_end"]["Read1"]["base_by_cycle"],
																  json[key]["Paired_end"]["Read2"]["base_by_cycle"]]
				PE_json[key]["St_Paired_End_Quality_by_Cycle"] = [json[key]["Paired_end"]["Read1"]["qualities_by_cycle"],
																 json[key]["Paired_end"]["Read2"]["qualities_by_cycle"]]
				PE_json[key]["St_PE_in"] = json[key]["Paired_end"]["in"]

				stats_json[key]["St_PE_Fraction" + index] = (json[key]["Paired_end"]["out"] / total_frags) * 100
				stats_json[key]["St_R1_Q30" + index] = ( json[key]["Paired_end"]["Read1"]["total_Q30_basepairs"] / json[key]["Paired_end"]["Read1"]["basepairs_in"] ) * 100 
				stats_json[key]["St_R2_Q30" + index] = ( json[key]["Paired_end"]["Read2"]["total_Q30_basepairs"] / json[key]["Paired_end"]["Read2"]["basepairs_in"] ) * 100 

				overview_stats[key]["Q30_Fraction"] += ((json[key]["Paired_end"]["Read1"]["total_Q30_basepairs"] + json[key]["Paired_end"]["Read2"]["total_Q30_basepairs"]) / total_bps)
				overview_stats[key]["PE_Fraction"] = json[key]["Paired_end"]["out"] / total_frags
				overview_stats[key]["PE_Output_Reads"] = json[key]["Paired_end"]["out"]
				overview_stats[key]["PE_Output_Bps"] = json[key]["Paired_end"]["Read1"]["basepairs_out"] + json[key]["Paired_end"]["Read2"]["basepairs_out"]
				
			except:
				overview_stats[key]["Q30_Fraction"] += 0
				overview_stats[key]["PE_Fraction"] = 0	
				overview_stats[key]["PE_Output_Reads"] = 0
				overview_stats[key]["PE_Output_Bps"] = 0 		   

				del PE_json[key]


		# output dictionary, keys are section, value is function called for figure generation
		section = {"Table": self.table(stats_json, index),
				   "Overview": overview_stats}

		stats_json.clear()

		if len(PE_json.keys()) != 0:
			section["Read Length Histogram (Paried End)"] = self.histogram(PE_json, "St_PE_histogram")
			section["Base by Cycle (Paired End)"] = self.base_by_cycle(PE_json, "St_Paired_End_Base_by_Cycle", index)
			section["Quality by Cycle (Paired End)"] = self.quality_by_cycle(PE_json, "St_Paired_End_Quality_by_Cycle", index)

		PE_json.clear()

		#only executres if single read data is detected
		if len(SE_json.keys()) != 0:
			section["Read Length Histogram (Single End)"] = self.histogram(SE_json, "St_SE_histogram")
			section["Base by Cycle (Single End)"] = self.base_by_cycle(SE_json, "St_Single_End_Base_by_Cycle", index)
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(SE_json, "St_Single_End_Quality_by_Cycle", index)
		
	
		return section 

