
from collections import OrderedDict
import logging
import numpy as np
import math

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, scatter

class OverviewStats():

	def table(self, json, app_list):


		read_data = {}
		bp_data = {}

		read_config = {'table_title': 'Input Bp Reduction'}
		bp_config = {'table_title': 'Input Fragment Reduction'}

		read_reducers = ["Pipeline Input", "hts_SeqScreener", "hts_SuperDeduper", 
						 "hts_Overlapper", "hts_LengthFilter", "hts_Stats"]
		bp_reducers = ["Pipeline Input", "hts_AdapterTrimmer", "hts_CutTrim", 
						 "hts_NTrimmer", "hts_QWindowTrim", "hts_PolyATTrim", "hts_Stats"]

		# Table constructor. Just like the MultiQC docs.
		read_headers = OrderedDict()
		bp_headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']


		first_app = list(json.keys())[0]
		samples = list(json[first_app].keys())

		for samp in samples:

			read_temp = {}
			bp_temp = {}

			for app in app_list:

				if app in read_reducers:
					read_temp[app] = json[app][samp]["Input_Reads"]

				if app in bp_reducers:
					bp_temp[app] = json[app][samp]["Input_Bp"]


			read_data[samp] = read_temp
			bp_data[samp] = bp_temp


		for i in range(len(app_list)):

			app = app_list[i]
			read_description = "Number of Input Fragments for " + app
			bp_description = "Number of Input Bps for " + app
			color = color_rotations[i % 4]
			read_headers[app] = {'title': app, 'namespace': app, 'description': read_description, 'format': '{:,.0f}', 'scale': color}
			bp_headers[app] = {'title': app, 'namespace': app, 'description': bp_description, 'format': '{:,.0f}', 'scale': color}

		html = '<h4>  Input Fragment Reduction </h4>'
		if len(read_headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	
		else:	
			html += table.plot(read_data, read_headers, read_config)

		html += '<br>\n<h4>  Input Basepair Reduction </h4>'
		if len(bp_headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	
		else:	
			html += table.plot(bp_data, bp_headers, bp_config)

		return 	html


	def hts_pca(self, json):

		mds_plot = {}

		config = {'title': "HTStream: PCA Plot"}

		keys = list(json.keys())
		samples_list = list(json[keys[0]].keys())
		row_length = len(samples_list)


		data = [[] for x in range(row_length)]

		stats_order = []
		stats_bool = True

		for x in range(row_length):

			sample = samples_list[x]

			for key in keys:

				sample_json = json[key][sample]

				if key == "hts_Stats" or key == "Pipeline Input":

					total_frags = sample_json["Fragment_Section"]["in"]
					total_bp = sample_json["Fragment_Section"]["basepairs_in"]

					try:
						fraction_se = sample_json["Read_Breakdown"]["Single_end"] / total_frags # fraction SE
					except:
						fraction_se = 0

					gc_content = (sample_json["Fragment_Section"]["base_composition"]["G"] / total_bp) + (sample_json["Fragment_Section"]["base_composition"]["C"] / total_bp)
					n_content = sample_json["Fragment_Section"]["base_composition"]["N"] / total_bp

					temp = [
							total_frags, # fragments in 
							sample_json["total_Q30"] / total_bp, # fraction Q30
							sample_json["Read_Breakdown"]["Paired_end"] / total_frags, # fraction PE
							fraction_se, # fraction SE
							gc_content, # GC Content
							n_content # N Content
							]

					data[x] += temp

					if stats_bool == True:
						stats_order += ["Total Frags", "Q30 Faction", "Fraction PE", "Fraction SE", "GC", "N"]

				else:

					temp = []

					for k, v in sample_json.items():
						if k != "Input_Reads":
							temp.append(v)

							if stats_bool == True:
								stats_order.append(k)
					
					data[x] += temp		

			stats_bool = False

		data = np.asarray(data).T
			
		data = htstream_utils.pca(data)
			
		for x in range(row_length):

			mds_plot[samples_list[x]] = {"x": data[0, x],
										 "y": data[1, x]}


		html = "<hr><h4>  Sample PCA Plot </h4>\n"
		html += scatter.plot(mds_plot, config)
		
		return html


	def execute(self, json, app_list):

		html = "" 
		html = self.table(json, app_list)
		html += self.hts_pca(json)
			
		return html

