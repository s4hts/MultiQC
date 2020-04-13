from collections import OrderedDict
import logging, statistics

from multiqc import config
from multiqc.plots import linegraph, heatmap

#################################################

""" Stats submodule for HTStream charts and graphs """

#################################################

class Stats():

	def base_by_cycle(self, json, read):

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

		data_list = []

		for key in json.keys():

			data = {"A": {},
					"C": {},
					"G": {},
					"T": {},
					"N": {}}

			bases = json[key][read]["data"]
			positions = json[key][read]["col_names"]

			sample_color = None
			sample_max = 0

			for i in range(len(positions)):

				total = bases[0][i] + bases[1][i] + bases[2][i] + bases[3][i] + bases[4][i]

				y_value_list = [(bases[0][i] / total) * 100, (bases[1][i] / total) * 100, 
								(bases[2][i] / total) * 100, (bases[3][i] / total) * 100,
								(bases[4][i] / total) * 100]

				sample_max = max(y_value_list)

				data["A"][i] = y_value_list[0]
				data["C"][i] = y_value_list[1]
				data["G"][i] = y_value_list[2]
				data["T"][i] = y_value_list[3]
				data["N"][i] = y_value_list[4]

			if sample_max >= 0.6:
				sample_color = '#e6c3c3'
			elif sample_max >= 0.4:
				sample_color = '#e6dcc3'
			else:
				sample_color = '#c3e6c3'


			config["data_labels"].append({'name': key,'ylab': 'Percentage', 
										  'xlab': 'Cycle', 'yCeiling': 100, 'categories': True, 
										  'smooth_points_sumcounts': False})
			data_list.append(data)


		return linegraph.plot(data_list, config)


	def quality_by_cycle(self, json, read):

		plot_list = []

		line_config = {
				  'smooth_points_sumcounts': False,
				  'categories': True,
				  'title': "HTStream: Quality by Cycle"
				  	 }

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

		btn_id = "-".join(read.split(" ")[:2]).lower()

		wrapper_html = '<div class="btn-group hc_switch_group">\n'
		wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this)" id="htstream_qbc_line_{r}_btn">Linegraph</button>\n'.format(r=btn_id)
		wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this)" id="htstream_qbc_heat_{r}_btn">Heatmaps</button></div>\n'.format(r=btn_id)
		wrapper_html += "<hr>"

		html = '<div class="btn-group hc_switch_group">\n'

		first = True
		pid = ""

		line_data = {}

		for key in json.keys():

			line_data[key] = {}

			heat_pconfig["id"] = "htstream_" + key + "_heatmap"

			x_lab = json[key][read]["col_names"]
			y_lab = json[key][read]["row_names"][::-1]
			data = []

			quality_scores = json[key][read]["shape"][0]
			cycles = json[key][read]["shape"][-1]

			total = []
			
			for pos in range(cycles):
				total.append(sum([ score_list[pos] for score_list in json[key][read]["data"] ]))
				line_data[key][pos] = statistics.mean([ score_list[pos] for score_list in json[key][read]["data"] ])

			for score in range(quality_scores - 1, -1, -1):
				data.append([])

				for pos in range(cycles):
					data[-1].append(json[key][read]["data"][score][pos] / total[pos])

			if first == True:
				active = "active"
				hidediv = ""
				first = False
				plot_html = heatmap.plot(data, x_lab, y_lab, heat_pconfig)
			else:
				active = ""
				hidediv = 'style="display: none;"'
				heatmap.plot(data, x_lab, y_lab, heat_pconfig)

			name = key
			pid = "htstream_" + key + "_btn"

			plot_list.append(plot_html)

			html += '<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this)" id="{pid}">{n}</button>\n'.format(a=active, pid=pid, n=name)



		html += '</div>\n\n<br></br>\n\n'
		html += plot_html 


		wrapper_html += '<div id="htstream_qbc_line_{r}">'.format(r=btn_id)
		wrapper_html += linegraph.plot(line_data, line_config) + "</div>"
		wrapper_html += '<div id="htstream_qbc_heat_{r}" style="display:none;">'.format(r=btn_id)
		wrapper_html += html + "</div>"

		html = wrapper_html 

		return html



	def linegraph(self, json, SE_presence):

		config = {'title': "HTStream: Read Length Histogram",
				  'data_labels': [
								  {'name': "R1 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'},
								  {'name': "R2 histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'}
								  ]}

		if SE_presence == False:
			histograms = ["R1 histogram", "R2 histogram"]
		else:
			histograms = ["R1 histogram", "R2 histogram", "SE histogram"]
			config["data_labels"].append({'name': "SE histogram", 'ylab': 'Frequency', 'xlab': 'Read Lengths'})


		invariant_dict = {}
		html = ""
		data_list = []

		for i in range(len(histograms)):

			data = {}

			for key in json.keys():

				if len(json[key][histograms[i]]) == 1:

					read_name = histograms[i].split(" ")[0]
					info = "[ " + str(json[key][histograms[i]][0][0]) + " Read Length at 1.0 Frequency ]"

					try:
						invariant_dict[key] += ", " + read_name + ": " + info

					except:
						invariant_dict[key] = read_name + ": " + info

				else:

					data[key] = {}

					total = sum([ count[1] for count in json[key][histograms[i]]])

					for item in json[key][histograms[i]]:

						data[key][item[0]] = item[1] / total

			if len(data.keys()) != 0:
				data_list.append(data)

		if len(list(invariant_dict.keys())) != 0:

			for sample, info in invariant_dict.items():
				notice = sample + ": " + info + "<br />"
			
			# to include when more is known about handling invariance.
			html += str("<p>" + "TEMPORARY PLACE HOLDER FOR READ HISTOGRAM PLOTS <br />" + notice + "</p>")
			
		if len(data_list) != 0:
			html += linegraph.plot(data_list, config)

		return html




	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {}

			try:

				stats_json[key]["SE histogram"] = json[key][-1]["Single_end"]["readlength_histogram"]
				stats_json[key]["Single End Base by Cycle"] = json[key][-1]["Single_end"]["base_by_cycle"]
				stats_json[key]["Single End Quality by Cycle"] = json[key][-1]["Single_end"]["qualities_by_cycle"]
							   
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
				   "Density Plots": self.linegraph(stats_json, SE_presence), 
				   "Base by Cycle (Read 1)": self.base_by_cycle(stats_json, "Read 1 Base by Cycle"),
				   "Quality by Cycle (Read 1)": self.quality_by_cycle(stats_json, "Read 1 Quality by Cycle"),
				   "Base by Cycle (Read 2)": self.base_by_cycle(stats_json, "Read 2 Base by Cycle"),
				   "Quality by Cycle (Read 2)": self.quality_by_cycle(stats_json, "Read 2 Quality by Cycle")
				   }

		if SE_presence == True:

			section["Base by Cycle (Single End)"] = self.base_by_cycle(stats_json, "Single End Base by Cycle")
			section["Quality by Cycle (Single End)"] = self.quality_by_cycle(stats_json, "Single End Quality by Cycle")

	
		return section 

