from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import linegraph, heatmap

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def base_by_cycle(self, json, read):

		config = {'data_labels': [],
				  'smooth_points_sumcounts': False,
				  'yCeiling': 100,
				  'categories': True,
				  }

		data_list = []

		for key in json.keys():

			data = {"A": {},
					"C": {},
					"G": {},
					"T": {},
					"N": {}}

			bases = json[key][read]["data"]
			positions = json[key][read]["col_names"]

			for i in range(len(positions)):

				total = bases[0][i] + bases[1][i] + bases[2][i] + bases[3][i] + bases[4][i]

				data["A"][i] = (bases[0][i] / total) * 100
				data["C"][i] = (bases[1][i] / total) * 100
				data["G"][i] = (bases[2][i] / total) * 100
				data["T"][i] = (bases[3][i] / total) * 100
				data["N"][i] = (bases[4][i] / total) * 100

			config["data_labels"].append({'name': key, 'ylab': 'Percentage', 'xlab': 'Cycle', 'yCeiling': 100, 'categories': True, 'smooth_points_sumcounts': False})
			data_list.append(data)

		return linegraph.plot(data_list, config)


	def quality_by_cycle(self, json, read):

		plot_list = []

		pconfig = {'id' : "",
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

		html = '<div class="btn-group hc_switch_group">\n'

		first = True
		pid = ""

		for key in json.keys():

			pconfig["id"] = "htstream_" + key + "_heatmap"

			x_lab = json[key][read]["col_names"]
			y_lab = json[key][read]["row_names"][::-1]
			data = []

			quality_scores = json[key][read]["shape"][0]
			cycles = json[key][read]["shape"][-1]

			total = []

			for pos in range(cycles):
				total.append(sum([ score_list[pos] for score_list in json[key][read]["data"] ]))

			for score in range(quality_scores - 1, -1, -1):
				data.append([])

				for pos in range(cycles):
					data[-1].append(json[key][read]["data"][score][pos] / total[pos])

			if first == True:
				active = "active"
				hidediv = ""
				first = False
				plot_html = heatmap.plot(data, x_lab, y_lab, pconfig)
			else:
				active = ""
				hidediv = 'style="display: none;"'
				heatmap.plot(data, x_lab, y_lab, pconfig)

			name = key
			pid = "htstream_" + key + "_btn"

			plot_list.append(plot_html)

			html += '<button class="btn btn-default btn-sm {a}" onclick="htstream_switch(this)" id="{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name)

		html += '</div>\n\n<br></br>\n\n'
		html += plot_html 

		return html



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

				try: 
					json[key][histograms[i]]

				except:
					break

				if len(data.keys()) == 0:
					data[key] = {}

					total = sum([ count[1] for count in json[key][histograms[i]]])

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

			stats_json[key] = {}

			try:

				stats_json[key]["SE histogram"] = json[key][-1]["Single_end"]["readlength_histogram"]
				stats_json[key]["Single Base by Cycle"] = json[key][-1]["Single_end"]["base_by_cycle"]
				stats_json[key]["Single Quality by Cycle"] = json[key][-1]["Single_end"]["qualities_by_cycle"]
							   
				SE_presence = True

			except:

				SE_presence = False
				pass

			stats_json[key]["R1 histogram"] = json[key][-1]["Paired_end"]["Read1"]["readlength_histogram"]
			stats_json[key]["R2 histogram"] = json[key][-1]["Paired_end"]["Read2"]["readlength_histogram"]
			stats_json[key]["Read 1 Base by Cycle"] = json[key][-1]["Paired_end"]["Read1"]["base_by_cycle"]
			stats_json[key]["Read 2 Base by Cycle"] = json[key][-1]["Paired_end"]["Read2"]["base_by_cycle"]
			stats_json[key]["Read 1 Quality by Cycle"] = json[key][-1]["Paired_end"]["Read1"]["qualities_by_cycle"]
			stats_json[key]["Read 2 Quality by Cycle"] =json[key][-1]["Paired_end"]["Read2"]["qualities_by_cycle"]


		section = {
				   "Density Plots": self.graph(stats_json), 
				   "Base by Cycle (Read 1)": self.base_by_cycle(stats_json, "Read 1 Base by Cycle"),
				   "Quality by Cycle (Read 1)": self.quality_by_cycle(stats_json, "Read 1 Quality by Cycle"),
				   "Base by Cycle (Read 2)": self.base_by_cycle(stats_json, "Read 2 Base by Cycle"),
				   "Quality by Cycle (Read 2)": self.quality_by_cycle(stats_json, "Read 2 Quality by Cycle")
				   }

		if SE_presence == True:

			section["Base by Cycle (Single End)"] = self.base_by_cycle(stats_json, "Single Base by Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(stats_json, "Single Quality by Cycle")

	
		return section 

