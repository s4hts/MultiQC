from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import linegraph, heatmap

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def base_by_cycle_R1(self, json, read):

		config = {'data_labels': [],
				  'extra_series': []}

		data_list = []

		for key in json.keys():

			config["data_labels"].append({'name': key, 'ylab': 'Percentage', 'xlab': 'Cycle'})
			config["extra_series"].append([])

			series_list = []

			C_series_dict = {'name': "C", 'data': []}
			G_series_dict = {'name': "G", 'data': []}
			T_series_dict = {'name': "T", 'data': []}
			N_series_dict = {'name': "N", 'data': []}

			series_list = [C_series_dict, G_series_dict, T_series_dict, N_series_dict]

			data = {"A": {}}
			total = []


			for pos in range(json[key][read]["shape"][-1]):
				total.append(sum(json[key][read]["data"][sublist][pos] for sublist in range(5)) / 100)

			bases = json[key][read]["data"]
			positions = json[key][read]["col_names"]

			for i in range(len(positions)):

				data["A"][positions[i]] = bases[0][i] / total[i]
				C_series_dict['data'].append([positions[i], bases[1][i] / total[i] ])
				G_series_dict['data'].append([positions[i], bases[2][i] / total[i] ])
				T_series_dict['data'].append([positions[i], bases[3][i] / total[i] ])
				N_series_dict['data'].append([positions[i], bases[4][i] / total[i] ])


			data_list.append(data)
			config['extra_series'][-1] = series_list


		return linegraph.plot(data_list, config)

	def quality_by_cycle(self, json, read):

		pconfig = {'yTitle': 'Q Score',
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

		for key in json.keys():

			x_lab = json[key][read]["col_names"]
			y_lab = json[key][read]["row_names"][::-1]
			data = []

			quality_scores = json[key][read]["shape"][0]
			cycles = json[key][read]["shape"][-1]

			for score in range(quality_scores - 1, -1, -1):
				data.append([])

				for pos in range(cycles):
					data[-1].append(json[key][read]["data"][score][pos] / cycles)


			break


		return heatmap.plot(data, x_lab, y_lab, pconfig)


	def graph(self, json):

		histograms = ["R1 histogram", "R2 histogram", "SE histogram"]

		config = {'data_labels': [
								  {'name': "R1 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "R2 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "SE histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								 ],
				  'extra_series': [[], [], []]}

		data_list = []

		for i in range(len(histograms)):

			data = {}

			for key in json.keys():

				if len(data.keys()) == 0:
					data[key] = {}

					total = sum([ count[1] for count in json[key][histograms[i]] ])

					for item in json[key][histograms[i]]:

						data[key][item[0]] = item[1] / total

					data_list.append(data)

				else:
					series_dict = {
								   'name': key,
	        					   'data': []
	        					  }

					for item in json[key][histograms[i]]:
						temp = [item[0], (item[1] / total)]
						series_dict['data'].append(temp)

					config['extra_series'][i].append(series_dict)

		return linegraph.plot(data_list, config)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "R1 histogram": json[key][1]["Single_end"]["readlength_histogram"],
			 				   "R2 histogram": json[key][1]["Paired_end"]["Read1"]["readlength_histogram"],
			 				   "SE histogram": json[key][1]["Paired_end"]["Read2"]["readlength_histogram"],
			 				   "Read 1 Base by Cycle": json[key][1]["Paired_end"]["Read1"]["base_by_cycle"],
			 				   "Read 2 Base by Cycle": json[key][1]["Paired_end"]["Read2"]["base_by_cycle"],
			 				   "Single Base by Cycle": json[key][1]["Single_end"]["base_by_cycle"],
			 				   "Read 1 Quality by Cycle": json[key][1]["Paired_end"]["Read1"]["qualities_by_cycle"],
			 				   "Read 2 Quality by Cycle": json[key][1]["Paired_end"]["Read2"]["qualities_by_cycle"],
			 				   "Single Quality by Cycle": json[key][1]["Single_end"]["qualities_by_cycle"]
						 	  }


		section = {
				   "Base by Cycle (Read 1)": self.base_by_cycle_R1(stats_json, "Read 1 Base by Cycle"),
				   "Base by Cycle (Read 2)": self.base_by_cycle_R1(stats_json, "Read 2 Base by Cycle"),
				   "Base by Cycle (Single End)": self.base_by_cycle_R1(stats_json, "Single Base by Cycle"),
				   "Quality by Cycle" : self.quality_by_cycle(stats_json, "Read 1 Quality by Cycle"),
				   "Density Plots": self.graph(stats_json)
				   }

		return section 

