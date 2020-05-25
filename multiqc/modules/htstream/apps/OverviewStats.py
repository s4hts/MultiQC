
from collections import OrderedDict
import logging
import numpy as np
import math

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, scatter

class OverviewStats():

	def table(self, json, app_list):

		print(json)
		return ""

		config = {'table_title': 'Input Fragment Reduction'}
		read_reducers = {"Pipeline Input", "hts_SeqScreener", "hts_SuperDeduper", 
						 "hts_Overlapper", "hts_LengthFilter", "hts_Stats"}

		# Table constructor. Just like the MultiQC docs.
		headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']

		html = '<h4>  Input Fragment Reduction </h4>'

		for key in json.keys():

			temp = {}

			for app in app_list:

				temp[app] = json[key][app]["InputFragments"]

			json[key] = temp

		table_data = False
		for i in range(len(app_list)):

			app = app_list[i]

			if app in read_reducers:

				description = "Number of Input Fragments for " + app
				color = color_rotations[i % 4]

				headers[app] = {'title': app, 'namespace': app, 'description': description, 'format': '{:,.0f}', 'scale': color}


		if len(headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	

		else:	
			html += table.plot(json, headers, config)

		return 	html


	def hts_mds(self, json):

		mds_plot = {}

		config = {'title': "HTStream: MDS Plot"}

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

					temp = [v for k,v in sample_json.items()]
					data[x] += temp

					if stats_bool == True:
						stats_order += [k for k,v in sample_json.items()]

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
		#html = self.table(json, app_list)
		html += self.hts_mds(json)
			
		return html

