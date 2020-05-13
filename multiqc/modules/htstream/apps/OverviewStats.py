
from collections import OrderedDict
import logging
from sklearn import manifold
import numpy as np

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, scatter

class OverviewStats():

	def table(self, json, app_list):

		config = {'table_title': 'Input Fragment Reduction'}

		# Table constructor. Just like the MultiQC docs.
		headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']

		html = '<h4>  Input Fragment Reduction </h4>'

		for key in json.keys():

			temp = {}

			for app in app_list:

				temp[app] = json[key][app]["InputFragments"]

			json[key] = temp


		for i in range(len(app_list)):

			app = app_list[i]

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
		length_keys = len(json.keys())

		data = np.zeros((length_keys, 5))


		for x in range(length_keys):

			sample = json[keys[x]]

			data[x,:] = [sample["total_Q30"]["Read1"],
						 sample["total_Q30"]["Read2"],
						  sample["total_Q30"]["Single_end"],
						 sample["Read_Breakdown"]["Paired_end"],
						 sample["Read_Breakdown"]["Single_end"]]
			

		data = manifold.MDS(n_components=2).fit(data).embedding_

		for x in range(length_keys):

			mds_plot[keys[x]] = {"x": data[x,0] // 100000000,
								 "y": data[x,1] // 100000000}


		html = "<hr><h4>  Random Data MDS Plot </h4>\n"
		html += scatter.plot(mds_plot, config)
		
		return html


	def execute(self, json, app_list):

			html = self.table(json, app_list)
			#html += self.hts_mds(json)
		
			return html

