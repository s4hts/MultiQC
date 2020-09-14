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

	def __init__(self):
		self.info = "Generates a JSON formatted file containing a set of statistical measures about the input read data."
		self.type = "both"

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



	def base_by_cycle(self, json, read):

		read_keys = {"St_PE_Base_by_Cycle": "PE",
					 "St_SE_Base_by_Cycle": "SE",
					 "PE": "Paired End",
				     "SE": "Single End"}


		read_code = read_keys[read]
		unique_id = str(random() % 1000)[5:]


		config = {'id': "htstream_stats_entropy_" + read + "_" + unique_id,
				  'title': "HTStream: Base by Cycle (" + read_code + ")",
				  'smooth_points_sumcounts': False,
				  'ylab': "Avg. Difference from 25%",
				  'xlab': "Cycle",
				  'categories': True,
				  'tt_decimals': '{:,.2f}',
				  'yPlotBands': [
								{'from': 0, 'to': 8, 'color': '#c3e6c3'},
								{'from': 8, 'to': 35, 'color': '#e6dcc3'},
								{'from': 35, 'to': 100, 'color': '#e6c3c3'},
								]

				  }

		samp_config = {'id': "htstream_stats_base_line_" + read + "_" + unique_id,
					  'title': "HTStream: Base by Cycle (" + read_code + ")",
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


		if read_code == "PE":

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


		line_data = {}
		data_list = []

		for samp in json.keys():

			line_data[samp] = {}
			samp_data = {"Base: A": {},
						 "Base: C": {},
						 "Base: G": {},
						 "Base: T": {},
						 "Base: N": {}}


			if read_code == "PE":
				data = [ json[samp][read][0]["data"][x] + json[samp][read][1]["data"][x] for x in range(5) ]


			else:
				data = json[samp][read]["data"]


			bases = len(data[0])

			for i in range(bases):

				total = sum([base[i] for base in data])

				samp_data["Base: A"][i + 1] = (data[0][i] / total) * 100
				samp_data["Base: C"][i + 1] = (data[1][i] / total) * 100
				samp_data["Base: G"][i + 1] = (data[2][i] / total) * 100
				samp_data["Base: T"][i + 1] = (data[3][i] / total) * 100
				samp_data["Base: N"][i + 1] = (data[4][i] / total) * 100
								
				avg = sum([ abs(((data[x][i] / total) * 100) - 25) for x in range(4) ]) / 4

				line_data[samp][i + 1] = avg



			# this config file is for the individual line of the multiline graph
			samp_config["data_labels"].append({'name': samp, 'ylab': 'Percentage', 
											   'xlab': 'Cycle', 'yCeiling': 100, 'categories': True, 
											   'smooth_points_sumcounts': False})

			# append base by cycle to data for this to data list
			data_list.append(samp_data)


		header_html = '<h4> Base by Cycle: ' + read_keys[read_code] + '</h4>'
		header_html += '''<p> Provides a measure of the uniformity of a distribution. The higher the average is at a certain position,
							the more unequal the base pair composition. N's are excluded from this calculation. </p>'''

		btn_label_1 = "Avg. Diff. from  25%"
		btn_label_2 = "Base by Cycle"


		line_1_id = "htstream_stats_entropy_{r}_{b}".format(r=read_code, b=unique_id)
		line_2_id = "htstream_stats_base_line_{r}_{b}".format(r=read_code, b=unique_id)


		line_1 = linegraph.plot(line_data, config)
		line_2 = linegraph.plot(data_list, samp_config)


		html = htstream_utils.multi_plot_html(header_html,
											  btn_label_1, btn_label_2,
											  line_1_id, line_2_id,
											  line_1, line_2)


		return html



	def quality_by_cycle(self, json, read):

	
		read_keys = {"St_PE_Quality_by_Cycle": "PE",
					 "St_SE_Quality_by_Cycle": "SE",
					 "PE": "Paired End",
				     "SE": "Single End"}


		read_code = read_keys[read]
		unique_id = str(random() % 1000)[5:]

		# config dictionary for mean Q score line graph
		line_config = {
				  'id': "htstream_stats_qbc_line_" + read_code + "_" + unique_id,
				  'smooth_points_sumcounts': False,
				  'categories': True,
				  'tt_decimals': '{:,.2f}',
				  'title': "HTStream: Mean Quality by Cycle (" + read_code + ")",
				  'xlab': "Cycle",
				  'ylab': "Mean Q Score",
				  'colors': {}
				  }


		# config dictionary for heatmaps
		heat_pconfig = {
				   'title': "HTStream: Quality by Cycle (" + read_code + ")",
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


		# check for uniform pe length, if so, add midpoint line to denote paired end reads.
		if read_code == "PE":

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
		


		line_data = {}
		first = True
		button_list = []


		for key in json.keys():

			# create dictionary for line graph. Again, format is {x: y}
			line_data[key] = {}

			# creates unique heatmap id that can be queired later by js.
			heat_pconfig["id"] = "htstream_stats_qbc_heat_" + read_code + "_" + key + "_" + unique_id

			if read_code == "PE":

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
			x_lab = [ str(int(x)) for x in json[key][read]["col_names"]]
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
				line_data[key][pos + 1] = total_score / temp_sum # total reads

				if line_data[key][pos + 1] > 30:
					num_above_q30 += 1


			# check to see what percent of bases have a mean Q score of at least 30
			q30_gate = (num_above_q30 / cycles) 

			if q30_gate < 0.6:
				line_config['colors'][key] = "#E16B6B"

			elif q30_gate < 0.8:
				line_config['colors'][key] = "#E8A243"

			else:
				line_config['colors'][key] = "#78D578"


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
			pid = heat_pconfig["id"] + "_btn"

			button_list.append('<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this)" id="{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name))


		# section heade
		header_html = '<h4> Quality by Cycle: ' + read_keys[read_code] + '</h4>'
		header_html += '''<p> Mean quality score for each position along the read. 
							  Sample is colored red if less than 60% of bps have mean score of at least Q30, 
							  orange if between 60% and 80%, and green otherwise.</p>'''

		btn_label_1 = "Mean Quality"
		btn_label_2 = "Quality by Cycle"


		line_1_id = "htstream_qbc_line_{r}_{u}".format(r=read_code, u=unique_id)
		line_2_id = "htstream_qbc_heat_{r}_{u}".format(r= read_code, u=unique_id)


		line_plot = linegraph.plot(line_data, line_config)
		heatmap_plot = htstream_utils.multi_heatmap_html(button_list, heatmap_html)

		html = htstream_utils.multi_plot_html(header_html,
											  btn_label_1, btn_label_2,
											  line_1_id, line_2_id,
											  line_plot, heatmap_plot,
											  exempt=False)

		return html




	def read_length(self, json, read):

		read_keys = {"St_PE_Read_Lengths": "PE",
					 "St_SE_Read_Lengths": "SE"}

		read_code = read_keys[read]
		unique_id = str(random() % 1000)[5:]


		# config dictionary for heatmaps
		heat_pconfig = {'id' : "htstream_stats_read_lengths_" + read_code + "_" + unique_id,
				   'title': "HTStream:  Read Length Heatmap (" +  read_code + ")",
				   'yTitle': 'Sample',
				   'xTitle': 'Length',
				   'square' : False,
				   'datalabels': False,
				   'max': 1.0, 
				   'colstops': [
					        [0, '#FFFFFF'],
					        [0.3, '#F8D527'],
					        [0.6, '#F8A627'],
					        [1, '#E70808']
					           ],
    			  }


		readlength_data = []
		lengths = []
		samples = []
		max_length = 0 
		min_length = 0

		for samp in json.keys():

			# paired end reads require the histograms be concatenated
			if read_code == "SE":
				read_lengths = json[samp][read][0]
				total = json[samp]["St_SE_in"]
				
			else:
				read_lengths = json[samp][read][0]
				total = json[samp]["St_PE_in"]
				r2 = [ [read_lengths[-1][0] + x[0], x[1]] for x in json[samp][read][1]] 

				# concat reads
				read_lengths += r2 

	
			data = []


			# populate data 
			for length in read_lengths:

				data.append(length[1] / total)


			# check if max read length is longest, if not update length of lists 
			if max_length < read_lengths[-1][0]:

				max_length = read_lengths[-1][0]

				for i in range(len(readlength_data)):

					temp = readlength_data[i]

					if len(temp) < (max_length - min_length): 

						readlength_data[i] = temp + [0] * (len(temp) - (max_length - min_length))  


			# check if min read length is longest
			if min_length > read_lengths[0][0] or min_length == 0:

				min_length = read_lengths[0][0]

				for i in range(len(readlength_data)):

					temp = readlength_data[i]

					if len(temp) < (max_length - min_length): 

						readlength_data[i] = [0] * (len(temp) - (max_length - min_length)) + temp  


			
			readlength_data.append(data)
			samples.append(samp)


		lengths = [i for i in range(min_length, max_length + 1)]


		if read_code == "PE":
			read = "Paired End"

		else:
			read = "Single End"


		html = '<h4> Read Lengths: ' + read + ' </h4>'
		html += '''<p> Distribution of read lengths for each sample. </p>'''


		html += heatmap.plot(readlength_data, lengths, samples, heat_pconfig)

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
				SE_json[key]["St_SE_Read_Lengths"] = [json[key]["Single_end"]["readlength_histogram"]]
				SE_json[key]["St_SE_Base_by_Cycle"] = json[key]["Single_end"]["base_by_cycle"]
				SE_json[key]["St_SE_Quality_by_Cycle"] = json[key]["Single_end"]["qualities_by_cycle"]
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
				PE_json[key]["St_PE_Read_Lengths"] = [json[key]["Paired_end"]["Read1"]["readlength_histogram"],
													  json[key]["Paired_end"]["Read2"]["readlength_histogram"]]
				PE_json[key]["St_PE_Base_by_Cycle"] = [json[key]["Paired_end"]["Read1"]["base_by_cycle"],
																  json[key]["Paired_end"]["Read2"]["base_by_cycle"]]
				PE_json[key]["St_PE_Quality_by_Cycle"] = [json[key]["Paired_end"]["Read1"]["qualities_by_cycle"],
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


		if len(PE_json.keys()) != 0:
			section["Read Length Histogram (Paried End)"] = self.read_length(PE_json, "St_PE_Read_Lengths")
			section["Base by Cycle (Paired End)"] = self.base_by_cycle(PE_json, "St_PE_Base_by_Cycle")
			section["Quality by Cycle (Paired End)"] = self.quality_by_cycle(PE_json, "St_PE_Quality_by_Cycle")

		
		#only executres if single read data is detected
		if len(SE_json.keys()) != 0:
			section["Read Length Histogram (Single End)"] = self.read_length(SE_json, "St_SE_Read_Lengths")
			section["Base by Cycle (Single End)"] = self.base_by_cycle(SE_json, "St_SE_Base_by_Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(SE_json, "St_SE_Quality_by_Cycle")
		

		return section 

