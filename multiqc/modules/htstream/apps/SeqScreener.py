from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################

class SeqScreener():

	def table(self, json, PE_presence, total):

		# Basic table constructor. See MultiQC docs.
		headers = OrderedDict()

		if total == 0:
			html = '<div class="alert alert-info"> No hits in any sample. </div>'	
			return html

		if PE_presence == True:
			headers["Ss_PE_loss"] = {'title': "% PE Lost", 'namespace': "% PE Lost",'description': 'Percentage of Paired End Reads Lost', 'format': '{:,.2f}', 
								 'suffix': "%", 'scale': 'Greens' }
			headers["Ss_PE_hits"] = {'title': "PE hits", 'namespace': 'PE hits','description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}

		headers["Ss_SE_in"] = {'title': "SE in", 'namespace': 'SE in', 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Ss_SE_out"] = {'title': "SE out", 'namespace': 'SE out','description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ss_SE_hits"] = {'title': "SE hits", 'namespace': 'SE hits', 'description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Ss_Notes"] = {'title': "Notes", 'namespace': 'Notes', 'description': 'Notes'}

		return table.plot(json, headers)



	def execute(self, json):

		stats_json = OrderedDict()
		overview_dict = {}

		total_hits = 0 

		for key in json.keys():

			try:
				perc_loss = ((json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"])  * 100
				PE_presence = True
				pe_hits = json[key]["Paired_end"]["hits"]


			except:
				perc_loss = 0
				pe_hits = 0
				PE_presence = False 


			try:
				se_hits = json[key]["Single_end"]["hits"]

			except:
				se_hits = 0

			total_hits += pe_hits + se_hits

			overview_dict[key] = {
								  "Reads_Lost": json[key]["Fragment"]["out"] / json[key]["Fragment"]["in"]
								 }

			# sample entry for stats dictionary
			stats_json[key] = {
			 				   "Ss_PE_loss": perc_loss,
							   "Ss_PE_hits": pe_hits,
							   "Ss_SE_in" : json[key]["Single_end"]["in"],
							   "Ss_SE_out": json[key]["Single_end"]["out"],
							   "Ss_SE_hits": se_hits,
							   "Ss_Notes": json[key]["Program_details"]["options"]["notes"],
						 	  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json, PE_presence, total_hits),
				   "Overview": overview_dict}

		return section

