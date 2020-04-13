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

		headers["PE in"] = {'namespace': "PE in", 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'namespace': "PE out", 'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["PE hits"] = {'namespace': "PE hits", 'description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["SE in"] = {'namespace': "SE in", 'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE hits"] = {'namespace': "SE hits", 'description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Blues'}
		headers["Notes"] = {'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)

	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["in"],
							   "PE out": json[key]["Paired_end"]["out"],
							   "PE hits": json[key]["Paired_end"]["hits"],
							   "SE in" : json[key]["Single_end"]["in"],
							   "SE out": json[key]["Single_end"]["out"],
							   "SE hits": json[key]["Single_end"]["hits"],
							   "Notes": json[key]["Program_details"]["options"]["notes"],
						 	  }

		section = {
				   "Table": self.table(stats_json)
				   }

		return section

