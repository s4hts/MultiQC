from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table

#################################################

""" SeqScreener submodule for HTStream charts and graphs """

#################################################

class SeqScreener():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["PE hits"] = {'description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE hits"] = {'description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "PE hits": json[key]["Paired_end"]["PE_hits"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
							   "SE hits": json[key]["Single_end"]["SE_hits"],
							   "Notes": json[key]["Notes"],
						 	  }

		section = {
				   "Table": self.table(stats_json)
				   }

		return section

