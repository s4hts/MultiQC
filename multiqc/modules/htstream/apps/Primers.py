from collections import OrderedDict
import logging
from random import random

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, heatmap

#################################################

""" Primers submodule for HTStream charts and graphs """

#################################################


class Primers():


	def __init__(self):
		self.info = "Identifies primer sequences located on the 5' ends of R1 and R2, or 5' and 3' end of SE reads."
		self.type = "bp_reducer"		


	def table(self, json, index, total_flipped):

		# standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["Pr_%_BP_Lost" + index] = {'title': "% Bp Lost",
										   'namespace': "% Bp Lost",
										   'description': 'Percentage of bps lost.',
										   'suffix': '%',
										   'format': '{:,.4f}',
										   'scale': 'Greens'}

		if total_flipped != 0:
			headers["Pr_Reads_Flipped" + index] = {'title': "Reads Flipped", 'namespace': "Reads Flipped", 'description': 'Number of Flipped Reads', 'format': '{:,.0f}', 'scale': 'Blues'}
		
		headers["Pr_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)


	def heatmap(self, json, index):

		# config dictionary for heatmaps
		heat_pconfig = {'id' : "htstream_primers_bargraph_" + index,
					   'title': "HTStream: Primers Heatmap",
					   'square' : False,
					   'datalabels': False,
					   'colstops': [
						        [0, '#FFFFFF'],
						        [1, '#1DC802']
						           ]
    			  }


		btn_id = "primers_" + index 
		unique_id = str(random() % 1000)[2:]
		first = True
		button_list = []
	
		for key in json.keys():

			# creates unique heatmap id that can be queired later by js.
			heat_pconfig["id"] = "htstream_" + btn_id + "_" + key + "_" + unique_id + "_heatmap_" + index

			data = []
			labs = []
			counts_list = json[key]["Pr_Primer_Counts" + index]

			for x in range(len(counts_list)):
				temp = counts_list[x]
				labs += temp[:-1]

			labs = list(set(labs))

			data = [ [0] * len(labs) for i in range(len(labs)) ] 

			for x in range(len(counts_list)):
				x_pos = labs.index(counts_list[x][0])
				y_pos = labs.index(counts_list[x][1])
				data[x_pos][y_pos] = counts_list[x][-1]
				data[y_pos][x_pos] = counts_list[x][-1]
			

			# if this is the first sample process, lucky them, they get to be shown first and marked as active.
			#	This step is necessary otherwise, the plot div is not initialized. The additional calls to the 
			#	heatmap function are simply to add the data to the internal jsons used by MultiQC.
			if first == True:
				active = "active" # button is default active
				first = False # shuts off first gat
				heatmap_html = heatmap.plot(data, labs, labs, heat_pconfig)

			else:
				active = "" # button is default off 
				heatmap.plot(data, labs, labs, heat_pconfig)


			# html div attributes and text
			name = key
			pid = "htstream_" + btn_id + "_" + key + "_" + unique_id + "_btn"

			button_list.append('<button class="btn btn-default btn-sm {a}" onclick="htstream_div_switch(this, {i})" id="{pid}">{n}</button>\n'.format(a=active, i=index, pid=pid, n=name))


		heatmap_plot = htstream_utils.multi_heatmap_html(button_list, heatmap_html)


		wrapper_html = '<h4> Primers: Primer Counts </h4>'
		wrapper_html += '''<p>Heatmap indicating abundance of primer combinations.</p>'''

		# Heatmaps
		wrapper_html  += '''<div class="mqc_hcplot_plotgroup">'''
		wrapper_html += '<div id="htstream_heat_primers_{u}" class="htstream_fadein">'.format(u=unique_id)
		wrapper_html += heatmap_plot + "</div></div>"

		final_html = wrapper_html 

		return wrapper_html


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		total_flipped = 0

		for key in json.keys():

			reads_lost = ( json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"] ) / json[key]["Fragment"]["in"]
			bp_lost = ( json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"] )
			perc_bp_lost = bp_lost / json[key]["Fragment"]["basepairs_in"]

			total_flipped += json[key]["Fragment"]["flipped"]

			overview_dict[key] = {
								  "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"] + json[key]["Paired_end"]["Read2"]["basepairs_out"],
								  "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
								  "Fraction_Bp_Lost": bp_lost / json[key]["Fragment"]["basepairs_in"]
								  }

			stats_json[key] = {
							   "Pr_%_BP_Lost" + index: perc_bp_lost * 100,
							   "Pr_Primer_Counts" + index: json[key]["Fragment"]['primers_counts'],
							   "Pr_Reads_Flipped" + index: json[key]["Fragment"]["flipped"],
							   "Pr_Notes" + index: json[key]["Program_details"]["options"]["notes"],
							  }


		# dictionary for sections and figure function calls
		section = {
				   "Table": self.table(stats_json, index, total_flipped),
				   "Primer Counts": self.heatmap(stats_json, index),
				   "Overview": overview_dict
				   }

		return section
