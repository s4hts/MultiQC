
from collections import OrderedDict
import logging
import numpy as np
import math

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, linegraph, scatter


class OverviewStats():

	def composition_and_reduction(self, json, app_list, data_type):

		line_config = {
					  'smooth_points_sumcounts': False,
					  'categories': True,
					  'yCeiling': 100,
					  'xlab': "Tool",
					  'ylab': "Percentage",
					  'colors': {
				  			 "SE": "#1EC2D0",
				  			 "PE": "#E8961B",
				  			 	},
					  'data_labels': []
					  }

		if data_type == "read":

			config = {'table_title': 'Fragment Reduction', 'id': "htstream_overview_read_reduction"}
			line_config['title'] = "HTStream: Fragment Composition"
			reducers = ["hts_SeqScreener", "hts_SuperDeduper", 
						"hts_Overlapper", "hts_LengthFilter", "hts_Stats"]
			table_suffix = " (read)"
			line_index = "_reads_out"
			index = "Reads"
			description = "Number of Output Fragments for "
			html_title = " Fragment Reduction "
			notice = "No Read Reducing Apps were found."

		else:

			config = {'table_title': 'Basepair Reduction', 'id': "htstream_overview_bp_reduction"}
			line_config['title'] = "HTStream: Basepair Composition"
			reducers = ["hts_AdapterTrimmer", "hts_CutTrim", 
						"hts_NTrimmer", "hts_QWindowTrim", "hts_PolyATTrim", "hts_Stats"]
			table_suffix = " (bps)"
			line_index = "_bps_out"
			index = "Bp"
			description = "Number of Output Bps for "
			html_title = " Basepair Reduction "
			notice = "No Read Reducing Apps were found."

		table_data = {}
		line_data_list = []

		# Table constructor. Just like the MultiQC docs.
		headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']

		first_app = list(json.keys())[0]
		samples = list(json[first_app].keys())

		app_list = ["Pipeline Input"] + app_list
		app_subsset = ["Pipeline Input"]

		for samp in samples:

			temp = {}
			line_data_temp = {"SE": {}, "PE": {}}

			for app in app_list:

				include = False

				if app[:-2] in reducers:
					total = json[app][samp]["Output_" + index]
					temp[app + table_suffix] = total
					app_subsset.append(app)
					include = True

				elif app == "Pipeline Input":
					temp[app + table_suffix] = json[app][samp]["Input_" + index]
					include = True


				if include == True:

					try:
						line_data_temp["SE"][app] = json[app][samp]["SE" + line_index] 
					except:
						line_data_temp["SE"][app] = 0

					try:
						line_data_temp["PE"][app] = json[app][samp]["PE" + line_index]
					except:
						line_data_temp["PE"][app] = 0


			table_data[samp] = temp
			line_data_list.append(line_data_temp)
			line_config['data_labels'].append({"name": samp, 'ylab': 'Percentage', 'xlab': 'Tool'})


		for i in range(len(app_subsset)):

			app = app_subsset[i]
			color = color_rotations[i % 4]
			headers[app + table_suffix] = {'title': app, 'namespace': app, 'description': description + app,
										   'format': '{:,.0f}', 'scale': color}


		title = '<h4> {t} </h4>'.format(t=html_title)

		if len(headers.keys()) < 2:
			html = title + "\n<br>"
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	
			return html

		else:	
			table_html = table.plot(table_data, headers, config)
			line_html = linegraph.plot(line_data_list, line_config)


		html = htstream_utils.composition_html(title, table_html, line_html, data_type) 

		return 	html


	def hts_pca(self, json):

		pca_plot = {}

		keys = list(json.keys())
		samples_list = list(json[keys[0]].keys())
		row_length = len(samples_list)

		exclude_list = ["Output_Reads", "Output_Bp", 
						"PE_reads_out", "PE_bps_out",
						"SE_reads_out", "SE_bps_out"]

		special_list = ["Overlap_Length_Max", "Overlap_Length_Med"]

		data = [[] for x in range(row_length)]

		stats_order = []
		stats_bool = True

		for x in range(row_length):

			sample = samples_list[x]

			for key in keys:

				if key != "Pipeline Input":
					sample_json = json[key][sample]
					temp = []

					for k, v in sample_json.items():
						if k not in exclude_list:
							temp.append(v)

							if stats_bool == True:
								stats_order.append(str(key + ": " + k))

					data[x] += temp

			stats_bool = False

		# prepe matrix
		data = np.asarray(data).T
		
		# normalize 
		data, stats_order, raw_data = htstream_utils.normalize(data, samples_list, stats_order, special_list)

		# pca 
		data, loadings, pc_perc = htstream_utils.pca(data, stats_order)			


		for x in range(len(data[0])):

			pca_plot[samples_list[x]] = {"x": data[0, x],
										 "y": data[1, x]}


		config = {'title': "HTStream: PCA Plot",
				  'xlab': "PC1" + " ({:.2f}%)".format(pc_perc[0]),
				  'ylab': "PC2" + " ({:.2f}%)".format(pc_perc[1]),
				  'data_labels': [
									{'name': 'Samples', 'xlab': "PC1" + " ({:.2f}%)".format(pc_perc[0]),
														'ylab': "PC2" + " ({:.2f}%)".format(pc_perc[1])},
									{'name': 'Loadings', 'xlab': "PC1" + " ({:.2f}%)".format(pc_perc[0]),
														 'ylab': "PC2" + " ({:.2f}%)".format(pc_perc[1])}
								  ]}

		data = [pca_plot, loadings]

		html = "<hr><h4> PCA Plot </h4>\n"
		html += scatter.plot(data, config) +  "\n<br>"

		
		return html, raw_data


	def execute(self, json, app_list):


		read_html = self.composition_and_reduction(json, app_list, "read") +  "\n<br>"
		bps_html = self.composition_and_reduction(json, app_list, "bp")
		scatter_html, pca_data = self.hts_pca(json)

		html = scatter_html + read_html + bps_html 

		return html, pca_data

