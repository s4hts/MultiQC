from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" LengthFilter submodule for HTStream charts and graphs """

#################################################



# NOT YET IMPLEMENTED



class LengthFilter():

	def table(self, json, PE_presence, total, index):

	# Basic table constructor. See MultiQC docs.
		headers = OrderedDict()

		if total == 0:
			html = '<div class="alert alert-info"> No reads discarded in any sample. </div>'	
			return html

		if PE_presence == True:
			headers["Lf_PE_loss" + index] = {'title': "% PE Lost", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.2f}', 
								 'suffix': "%", 'scale': 'Greens' }

		headers["Lf_SE_in" + index] = {'title': "SE in", 'namespace': 'SE in', 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Lf_SE_out" + index] = {'title': "SE out", 'namespace': 'SE out','description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Lf_Notes" + index] = {'title': "Notes", 'namespace': 'Notes', 'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		total_loss = 0 

		for key in json.keys():

			try:
				perc_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"])  * 100
				frag_out = json[key]["Fragment"]["out"]
				reads_lost = json[key]["Fragment"]["out"] / json[key]["Fragment"]["in"]
				PE_presence = True


			except:
				perc_loss = 0
				frag_out = 0
				reads_lost = 0
				PE_presence = False 

			total_loss += perc_loss

			overview_dict[key] = {
								  "Output_Reads": frag_out,
								  "Reads_Lost": reads_lost
								 }

			# sample entry for stats dictionary
			stats_json[key] = {
			 				   "Lf_PE_loss" + index: perc_loss,
							   "Lf_SE_in" + index: json[key]["Single_end"]["in"],
							   "Lf_SE_out" + index: json[key]["Single_end"]["out"],
							   "Lf_Notes" + index: json[key]["Program_details"]["options"]["notes"],
						 	  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json, PE_presence, total_loss, index),
				   "Overview": overview_dict}


		return section 