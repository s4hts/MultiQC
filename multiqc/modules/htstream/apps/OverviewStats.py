from collections import OrderedDict
import logging, statistics

from multiqc import config
from multiqc.plots import table, linegraph

class OverviewStats():

	def table(self, json, app_list):

		config = {'table_title': 'Input Fragment Reduction'}

		# Table constructor. Just like the MultiQC docs.
		headers = OrderedDict()

		color_rotations = ['Greens', 'RdPu', 'Blues', 'Oranges']

		html = '<h4>  Input Fragment Reduction </h4>'

		for i in range(len(app_list)):

			app = app_list[i]
			header_title = app + "_InputFragments"
			description = "Number of Input Fragments for " + app
			color = color_rotations[i % 4]

			headers[header_title] = {'title': app, 'namespace': app, 'description': description, 'format': '{:,.0f}', 'scale': color}


		if len(headers.keys()) < 2:
			notice = "No Read Reducing Apps were found."
			html = '<div class="alert alert-info">{n}</div>'.format(n = notice)	

		else:	
			html += table.plot(json, headers, config)

		return 	html


	def execute(self, json, app_list):

			html = self.table(json, app_list)
		
			return html

