
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

		read_config = {'table_title': 'Fragment Reduction', 'id': "htstream_overview_read_reduction"}
		bp_config = {'table_title': 'Basepair Reduction', 'id': "htstream_overview_bp_reduction"}

		read_reducers = ["hts_SeqScreener", "hts_SuperDeduper", 
						 "hts_Overlapper", "hts_LengthFilter", "hts_Stats"]
		bp_reducers = ["hts_AdapterTrimmer", "hts_CutTrim", 
						 "hts_NTrimmer", "hts_QWindowTrim", "hts_PolyATTrim", "hts_Stats"]

		# Table constructor. Just like the MultiQC docs.
		read_headers = OrderedDict()
		bp_headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']

		first_app = list(json.keys())[0]
		samples = list(json[first_app].keys())

		app_list = ["Pipeline Input"] + app_list

		for samp in samples:

			read_temp = {}
			bp_temp = {}


			for app in app_list:

				if app[:-2] in read_reducers:
					read_temp[app + " (read)"] = json[app][samp]["Output_Reads"]

				if app[:-2] in bp_reducers:
					bp_temp[app + " (bps)"] = json[app][samp]["Output_Bp"]

				if app == "Pipeline Input":
					read_temp[app + " (read)"] = json[app][samp]["Input_Reads"]
					bp_temp[app + " (bps)"] = json[app][samp]["Input_Bp"]



			read_data[samp] = read_temp
			bp_data[samp] = bp_temp


		for i in range(len(app_list)):

			app = app_list[i]
			read_description = "Number of Output Fragments for " + app
			bp_description = "Number of Output Bps for " + app
			color = color_rotations[i % 4]
			read_headers[app + " (read)"] = {'title': app, 'namespace': app, 'description': read_description, 'format': '{:,.0f}', 'scale': color}
			bp_headers[app + " (bps)"] = {'title': app, 'namespace': app, 'description': bp_description, 'format': '{:,.0f}', 'scale': color}


		html = '<h4> Fragment Reduction </h4>'
		if len(read_headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	
		else:	
			html += table.plot(read_data, read_headers, read_config)


		html += '<br>\n<h4> Basepair Reduction </h4>'
		if len(bp_headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	
		else:	
			html += table.plot(bp_data, bp_headers, bp_config)

		return 	html


	def hts_pca(self, json):

		mds_plot = {}

		keys = list(json.keys())
		samples_list = list(json[keys[0]].keys())
		row_length = len(samples_list)

		exclude_list = ["Output_Reads", "Output_Bp"]
		special_list = ["Hist_Med", "Hist_Max"]

		data = [[] for x in range(row_length)]
		data_out = {}

		stats_order = []
		stats_bool = True

		for x in range(row_length):

			sample = samples_list[x]

			for key in keys:

				sample_json = json[key][sample]

				if "hts_Stats" in key:

					total_frags = sample_json["Output_Reads"]
					total_bp = sample_json["Output_Bp"]
					gc_content = sample_json["gc_content"] 
					n_content = sample_json["n_content"] 

					try:
						fraction_se = sample_json["Read_Breakdown"]["Single_end"] / total_frags # fraction SE
					except:
						fraction_se = 0

					try:
						fraction_pe = sample_json["Read_Breakdown"]["Paired_end"] / total_frags # fraction PE
					except:
						fraction_pe = 0

					temp = [
							sample_json["total_Q30"] / total_bp, # fraction Q30
							fraction_pe, 
							fraction_se, 
							gc_content, # GC Content
							n_content # N Content
							]
				
					if stats_bool == True:
						stats_order += [key + ": Q30 Fraction",
										key + ": Fraction PE",
										key + ": Fraction SE",
										key + ": GC Content",
										key + ": N Content"]

					data[x] += temp	

				elif key != "Pipeline Input":

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
		n, m = data.shape # rows, col
		to_delete = []

		# format dictionary for output pca stats (raw data)
		for x in range(len(samples_list)):
			data_out[samples_list[x]] = dict(zip(stats_order, data[:,x]))


		# normalize 
		for x in range(n):

			row = data[x,:]
			mean = np.mean(row)
			std = np.std(row)

			# remove rows with no variation, also, mean center and normalize variance
			if np.all(row == row[0]):
				to_delete.append(x)

			elif any(s in stats_order[x] for s in special_list):
				data[x,:] = (row - mean) / mean

			elif any(i > 1 for i in row):
				data[x,:] = (row - mean) / std 

			else:
				data[x,:] = (row - np.mean(row))



		# remove indeterminant columns
		to_delete = sorted(to_delete, reverse=True)
		for x in to_delete:	
			data = np.delete(data, x, 0)
			stats_order.remove(stats_order[x])

		
		# pca function
		data, loadings, pc_perc = htstream_utils.pca(data, stats_order)			


		for x in range(len(data[0])):

			mds_plot[samples_list[x]] = {"x": data[0, x],
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

		data = [mds_plot, loadings]

		html = "<hr><h4> PCA Plot </h4>\n"
		html += scatter.plot(data, config)

		
		return html, data_out


	def execute(self, json, app_list):

		table_html = self.table(json, app_list)
		scatter_html, pca_data = self.hts_pca(json)

		html = table_html + scatter_html

		return html, pca_data

