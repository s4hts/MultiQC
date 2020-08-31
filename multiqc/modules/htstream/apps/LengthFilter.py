from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" LengthFilter submodule for HTStream charts and graphs """

#################################################


class LengthFilter():


	def __init__(self):
		self.info = "Discards reads below a minimum length threshold."
		self.type = "read_reducer"
	

	def table(self, json, pe_total_loss, se_total_loss, pe_orphaned_total, index):

	# Basic table constructor. See MultiQC docs.
		headers = OrderedDict()

		if pe_total_loss == 0 and se_total_loss == 0 and pe_orphaned_total == 0:
			html = '<div class="alert alert-info"> No reads discarded in any sample. </div>'	
			return html

		if pe_total_loss != 0:
			headers["Lf_PE_loss" + index] = {'title': "% PE Lost", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.2f}', 
								 'suffix': "%", 'scale': 'Greens' }

		if pe_orphaned_total != 0:
			headers["Lf_PE_Orphaned" + index] = {'title': "% PE Orphaned", 'namespace': "% PE Orphaned", 'description': 'Percentage of Paired End Reads Orphaned (Now Single End)', 
											 'format': '{:,.2f}', 'suffix': "%", 'scale': 'RdPu' }

		if se_total_loss != 0:
			headers["Lf_SE_loss" + index] = {'title': "% SE Lost", 'namespace': "% SE Lost", 'description': 'Percentage of Single End Reads Lost', 'format': '{:,.2f}', 
											   'suffix': "%", 'scale': 'Blues' }


		headers["Lf_Notes" + index] = {'title': "Notes", 'namespace': 'Notes', 'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		pe_total_loss = 0
		se_total_loss = 0
		pe_orphaned_total = 0

		for key in json.keys():

			try:
				reads_lost = (json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"]) / json[key]["Fragment"]["in"]
			
			except:
				reads_lost = 0

			try:
				pe_perc_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"])  * 100
				pe_orphaned = ((json[key]["Paired_end"]["Read1"]["discarded"] + json[key]["Paired_end"]["Read2"]["discarded"]) / json[key]["Paired_end"]["in"])  * 100

			except:
				pe_perc_loss = 0
				pe_orphaned = 0

			if json[key]["Single_end"]["in"] == 0:
				se_perc_loss = 0

			else:
				se_perc_loss = (json[key]["Single_end"]["discarded"] / json[key]["Single_end"]["in"])  * 100


			pe_total_loss += pe_perc_loss
			se_total_loss += se_perc_loss
			pe_orphaned_total += pe_orphaned

			overview_dict[key] = {
								  "PE_Output_Reads": json[key]["Paired_end"]["out"],
								  "SE_Output_Reads": json[key]["Single_end"]["out"],
								  "Fraction_Reads_Lost": reads_lost
								 }

			# sample entry for stats dictionary
			stats_json[key] = {
			 				   "Lf_PE_loss" + index: pe_perc_loss,
							   "Lf_PE_Orphaned" + index: pe_orphaned,
							   "Lf_SE_loss" + index: se_perc_loss,
							   "Lf_Notes" + index: json[key]["Program_details"]["options"]["notes"],
						 	  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json, pe_total_loss, se_total_loss, pe_orphaned_total, index),
				   "Overview": overview_dict}


		return section 