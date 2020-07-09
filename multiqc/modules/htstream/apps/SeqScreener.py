from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################

class SeqScreener():

	def table(self, json, pe_total, se_total, index):

		# Basic table constructor. See MultiQC docs.
		headers = OrderedDict()

		if (pe_total + se_total) == 0:
			html = '<div class="alert alert-info"> <strong>Notice:</strong> No hits in any sample. </div>'	
			return html

		if pe_total != 0:
			headers["Ss_PE_hits" + index] = {'title': "PE hits", 'namespace': 'PE hits','description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
			headers["Ss_PE_%_hits" + index] = {'title': "% PE Hits", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.4f}', 
											   'suffix': "%", 'scale': 'Greens' }

		if se_total != 0:
			headers["Ss_SE_hits" + index] = {'title': "SE hits", 'namespace': 'SE hits','description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Greens'}
			headers["Ss_SE_%_hits" + index] = {'title': "% SE Hits", 'namespace': "% SE Lost",'description': 'Percentage of Single End Reads Lost', 'format': '{:,.4f}', 
											   'suffix': "%", 'scale': 'RdPu' }

		headers["Ss_Notes" + index] = {'title': "Notes", 'namespace': 'Notes', 'description': 'Notes'}

		return table.plot(json, headers)



	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		pe_total_hits = 0
		se_total_hits = 0 


		for key in json.keys():

			try:		
				pe_hits = json[key]["Paired_end"]["hits"]
				perc_pe_hits = (pe_hits / json[key]["Paired_end"]["in"])  * 100


			except:
				pe_hits = 0
				perc_pe_hits = 0 
				


			try:
				se_hits = json[key]["Single_end"]["hits"]
				perc_se_hits = (se_hits / json[key]["Single_end"]["in"]) * 100

			except:
				se_hits = 0
				perc_se_hits = 0

			pe_total_hits += pe_hits
			se_total_hits += se_hits


			overview_dict[key] = {
								  "Output_Reads": json[key]["Fragment"]["out"],
								  "PE_reads_out": (json[key]["Paired_end"]["out"] / json[key]["Fragment"]["out"] * 100),
								  "SE_reads_out": (json[key]["Single_end"]["out"] / json[key]["Fragment"]["out"]) * 100,
								  "Fraction_Reads_Lost": (json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"]) / json[key]["Fragment"]["in"],
								  "Percent_Hits": (pe_hits + se_hits) / json[key]["Fragment"]["in"]
								 }

			# sample entry for stats dictionary
			stats_json[key] = {
			 				   "Ss_PE_%_hits" + index: perc_pe_hits,
							   "Ss_PE_hits" + index: pe_hits,
							   "Ss_SE_%_hits"  + index: perc_se_hits,
							   "Ss_SE_hits" + index: se_hits,
							   "Ss_Notes" + index: json[key]["Program_details"]["options"]["notes"],
						 	  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json, pe_total_hits, se_total_hits, index),
				   "Overview": overview_dict}

		return section

