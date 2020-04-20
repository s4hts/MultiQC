from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################

class SeqScreener():

	def table(self, json):

		# Basic table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["Ss_PE_in"] = {'title': "PE in", 'namespace': 'PE in', 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Ss_PE_out"] = {'title': "PE out", 'namespace': 'PE out','description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ss_PE_hits"] = {'title': "PE hits", 'namespace': 'PE hits','description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Ss_SE_in"] = {'title': "SE in", 'namespace': 'SE in', 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Ss_SE_out"] = {'title': "SE out", 'namespace': 'SE out','description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Ss_SE_hits"] = {'title': "SE hits", 'namespace': 'SE hits', 'description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Ss_Notes"] = {'title': "Notes", 'namespace': 'Notes', 'description': 'Notes'}

		return table.plot(json, headers)



	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			# sample entry for stats dictionary
			stats_json[key] = {
			 				   "Ss_PE_in": json[key]["Paired_end"]["in"],
							   "Ss_PE_out": json[key]["Paired_end"]["out"],
							   "Ss_PE_hits": json[key]["Paired_end"]["hits"],
							   "Ss_SE_in" : json[key]["Single_end"]["in"],
							   "Ss_SE_out": json[key]["Single_end"]["out"],
							   "Ss_SE_hits": json[key]["Single_end"]["hits"],
							   "Ss_Notes": json[key]["Program_details"]["options"]["notes"],
						 	  }

		# sections and figure function calls
		section = {
				   "Table": self.table(stats_json)
				   }

		return section

