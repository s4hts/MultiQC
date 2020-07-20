
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
					  'id': "htstream_overview_" + data_type + "composition",
					  'smooth_points_sumcounts': False,
					  'logswitch': True,
					  'categories': True,
					  'xlab': "Tool",
					  'ylab': "Counts",
					  'tt_decimals': "{point.y:.0f}'",
					  'colors': {
				  			 "SE Reads": "#1EC2D0",
				  			 "PE Reads": "#E8961B",
				  			 	},
					  'data_labels': []
					  }

		if data_type == "read":

			config = {'table_title': 'Fragment Reduction', 'id': "htstream_overview_read_reduction"}
			line_config['title'] = "HTStream: Fragment Composition"
			reducers = ["hts_SeqScreener", "hts_SuperDeduper", 
						"hts_Overlapper", "hts_LengthFilter", "hts_Stats"]
			table_suffix = " (read)"
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
			index = "Bps"
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
		app_subset = ["Pipeline Input"]

		for samp in samples:

			temp = {}
			line_data_temp = {"SE "  + index: {}, 
							  "PE "  + index: {}}

			for app in app_list:

				include = False

				if app[:-2] in reducers:
					total = json[app][samp]["PE_Output_" + index] + json[app][samp]["SE_Output_" + index]
					app_subset.append(app)
					prefix = "Output_"
					include = True

				elif app == "Pipeline Input":
					total = json[app][samp]["PE_Input_" + index] + json[app][samp]["SE_Input_" + index]
					prefix = "Input_"
					include = True


				if include == True:

					temp[app + table_suffix] = total

					try:
						line_data_temp["SE " + index][app] = json[app][samp]["SE_" + prefix + index]
					except:
						line_data_temp["SE " + index][app] = 0

					try:
						line_data_temp["PE " + index][app] = json[app][samp]["PE_" + prefix + index]
					except:
						line_data_temp["PE "  + index][app] = 0


			table_data[samp] = temp
			line_data_list.append(line_data_temp)
			line_config['data_labels'].append({"name": samp, 'ylab': 'Counts', 'xlab': 'Tool'})


		for i in range(len(app_subset)):

			app = app_subset[i]
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


	def hts_line(self, json):

		line_config = {
					  'id': "htstream_overview_linechart",
					  'title': "HTStream: Preprocessing Statistics",
					  'smooth_points_sumcounts': False,
					  'categories': True,
					  'yCeiling': 1,
					  'xlab': "Tool",
					  'ylab': "Value",
					  'tt_decimals': "{point.y:.3f}'",  
					  'data_labels': []
					  }


		# pca_plot = {}

		keys = list(json.keys())
		samples_list = list(json[keys[0]].keys())
		row_length = len(samples_list)

		exclude_list = ["PE_Output_Reads", "PE_Output_Bps",
						"SE_Output_Reads", "SE_Output_Bps"]

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
		data = np.array(data).T

		# normalize 
		data, stats_order, raw_data = htstream_utils.normalize(data, samples_list, stats_order)

		data_dict = {}

		for x in range(len(samples_list)):

			samp = samples_list[x]
			data_dict[samp] = {}
			values = data[:,x]

			for y in range(len(values)):
				data_dict[samp][stats_order[y]] = values[y]


			line_config['data_labels'].append({"name": samp, 'ylab': 'Value', 'xlab': 'Tool'})


		html = "<hr><h4> Preprocessing Statistics </h4>\n"
		html += linegraph.plot(data_dict, line_config) +  "\n<br>"

		
		return html, raw_data


	def execute(self, json, app_list):


		read_html = self.composition_and_reduction(json, app_list, "read") +  "\n<br>"
		bps_html = self.composition_and_reduction(json, app_list, "bp")
		line_html, line_data = self.hts_line(json)

		html = line_html + read_html + bps_html 

		return html, line_data

