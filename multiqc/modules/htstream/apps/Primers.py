from collections import OrderedDict
import logging
from random import random

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, heatmap

#################################################

""" Primers submodule for HTStream charts and graphs """

#################################################

'''
##################
This tool is currently a work in progress 
##################
'''

class Primers():

	def table(self, json, index, total_flipped):

		# standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["Pr_%_BP_Lost" + index] = {'title': "% Bp Lost",
										   'namespace': "% Bp Lost",
										   'description': 'Percentage of bps lost.',
										   'suffix': '%',
										   'format': '{:,.2f}',
										   'scale': 'Greens'}
		headers["Pr_BP_Lost" + index] = {'title': "Bp Lost", 'namespace': "Bp Lost", 'description': 'Number of basepairs lost', 'format': '{:,.0f}', 'scale': 'RdPu'}

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

		html = htstream_utils.primers_heatmap_html(unique_id, button_list, heatmap_html)

		return html


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		total_flipped = 0

		for key in json.keys():

			bp_lost = ( json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"] )
			perc_bp_lost = bp_lost / json[key]["Fragment"]["basepairs_in"]

			total_flipped += json[key]["Fragment"]["flipped"]

			overview_dict[key] = {
								  "Output_Bp": json[key]["Fragment"]["basepairs_out"],
								  "Bp_Lost": bp_lost / json[key]["Fragment"]["basepairs_in"]
								  }

			stats_json[key] = {
							   "Pr_%_BP_Lost" + index: perc_bp_lost * 100,
							   "Pr_BP_Lost" + index: bp_lost,
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
