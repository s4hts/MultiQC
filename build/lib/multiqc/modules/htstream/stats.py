#!/usr/bin/env python

""" MultiQC submodule for HTStream charts and graphs """

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph, linegraph, table


#################################################
# HTStream App Classes

#################################################
class AdapterTrimmer():


	def table(self, json):

		headers = OrderedDict()

		headers["Reads in"] = {'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Reads out"] = {'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Adapters"] = {
						   'description': 'Percentage of Reads (SE and PE) with an Adapter',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Blues'
						  }

		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():
			
			adapter_reads = json[key]["Single_end"]["SE_adapterTrim"] + json[key]["Paired_end"]["PE_adapterTrim"]  

			stats_json[key] = {
							   "Reads in": json[key]["totalFragmentsInput"],
							   "Reads out": json[key]["totalFragmentsOutput"],
							   "% Adapters" : (adapter_reads / json[key]["totalFragmentsInput"]) * 100,
							   "Notes": json[key]["Notes"]
							  }

		section = {"Table": self.table(stats_json)}

		return section

#################################################
class CutTrim():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)

		section = {
				   "Table": plot,
				   }

		return section 

#################################################
class Overlapper():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json, inserts):

		if inserts == 0:
			return

		categories  = OrderedDict()

		categories['Sins'] = {
							  'name': 'Short Inserts',
							  'color': '#4d8de4'
							 }
		categories['Mins'] = {
							  'name': 'Medium Inserts',
							  'color': '#e57433'
							 }
		categories['Lins'] = {
							  'name': 'Long Inserts',
							  'color': '#33a02c'
							 }

		return bargraph.plot(json, categories)

	def linegraph(self, json):

		data_list = [] 
		config = {'data_labels': []}

		for key in json.keys():

			data = {}
			data[key] = {}
			config_subdict = {'name': key, 'ylab': 'Frequency', 'xlab': 'Overlapping Read Lengths'}				

			total = sum([ count[1] for count in json[key]["Histogram"]])

			for item in json[key]["Histogram"]:

				data[key][item[0]] = item[1] / total

			data_list.append(data)
			config['data_labels'].append(config_subdict)


		return linegraph.plot(data_list, config)

	def execute(self, json):

		stats_json = OrderedDict()

		inserts = 0

		for key in json.keys():

			stats_json[key] = {
			 			 	   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
						 	   "Notes": json[key]["Notes"],
							   "Sins": json[key]["sins"],
							   "Mins": json[key]["mins"],
							   "Lins": json[key]["lins"],
							   "Histogram": json[key]["readlength_histogram"]
							  }

			inserts += (json[key]["sins"] + json[key]["mins"] + json[key]["lins"])

		section = {
				   "Table": self.table(stats_json),
				   "Reads with Insertions": self.bargraph(stats_json, inserts),
				   "Overlapped Lengths Density Plots": self.linegraph(stats_json)
				   }

		return section

#################################################
class QWindowTrim():


	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def bargraph(self, json, reads):

		if reads == 0:
			return

		categories  = OrderedDict()

		categories['Left Trimmed Reads'] = {'name': 'Left Trimmed Reads'}
		categories['Right Trimmed Reads'] = {'name': 'Right Trimmed Reads'}

		return bargraph.plot(json, categories)

	def execute(self, json):

		stats_json = OrderedDict()

		trimmed_reads = 0

		for key in json.keys():

			lefttrimmed_reads = json[key]["Paired_end"]["R1_leftTrim"] + json[key]["Paired_end"]["R2_leftTrim"] + json[key]["Single_end"]["SE_leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["R1_rightTrim"] + json[key]["Paired_end"]["R2_rightTrim"] + json[key]["Single_end"]["SE_rightTrim"]

			trimmed_reads += (lefttrimmed_reads + rightrimmed_reads)

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
							   "Notes": json[key]["Notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads
						 	  }

		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_reads)
				   }

		return section

#################################################
class NTrimmer():


	def table(self, json):

		headers = OrderedDict()

		headers["Reads in"] = {'description': 'Number of Input Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Reads out"] = {'description': 'Number of Output Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Discarded"] = {
						   'description': 'Percentage of Reads (SE and PE) Discarded',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Blues'
						  }

		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)


	def bargraph(self, json, reads):

		if reads == 0:
			return

		categories  = OrderedDict()

		categories['Left Trimmed Reads'] = {'name': 'Left Trimmed Reads'}
		categories['Right Trimmed Reads'] = {'name': 'Right Trimmed Reads'}

		return bargraph.plot(json, categories)


	def execute(self, json):

		stats_json = OrderedDict()

		trimmed_reads = 0

		for key in json.keys():

			discarded_reads = json[key]["Single_end"]["SE_discarded"] + json[key]["Paired_end"]["R1_discarded"] + json[key]["Paired_end"]["R2_discarded"]  

			lefttrimmed_reads = json[key]["Paired_end"]["R1_leftTrim"] + json[key]["Paired_end"]["R2_leftTrim"] + json[key]["Single_end"]["SE_leftTrim"]
			rightrimmed_reads = json[key]["Paired_end"]["R1_rightTrim"] + json[key]["Paired_end"]["R2_rightTrim"] + json[key]["Single_end"]["SE_rightTrim"]

			trimmed_reads += (lefttrimmed_reads + rightrimmed_reads)

			stats_json[key] = {
			 				   "Reads in": json[key]["totalFragmentsInput"],
							   "Reads out": json[key]["totalFragmentsOutput"],
							   "% Discarded" : (discarded_reads / json[key]["totalFragmentsInput"]) * 100,
							   "Notes": json[key]["Notes"],
							   "Left Trimmed Reads": lefttrimmed_reads,
							   "Right Trimmed Reads": rightrimmed_reads,
							  }

		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_reads)
				   }

		return section

#################################################
class PolyATTrim():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)

		section = {
				   "Table": plot,
				   }

		return #section 

#################################################
class SeqScreener():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["PE hits"] = {'description': 'Number of Paired End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Oranges'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE hits"] = {'description': 'Number of Single End Reads with Sequence', 'format': '{:,.0f}', 'scale': 'Oranges'}
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

#################################################
class SuperDeduper():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["% Duplicates"] = {
						   'description': 'Percentage of Duplicate Reads (SE and PE)',
						   'suffix': '%',
						   'max': 100,
						   'format': '{:,.2f}',
						   'scale': 'Oranges'
						  }
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def linegraph(self, json):

		config = {'xlab': "Duplication Size", 'ylab': "Counts",
				  'extra_series': []}
		data = {}

		for key in json.keys():

			if len(data.keys()) == 0:
				data[key] = {}

				for item in json[key]["Saturation"]:

					data[key][item[0]] = item[1] 

			else:
				series_dict = {
							   'name': key,
        					   'data': []
        					  }

				for item in json[key]["Saturation"]:
					series_dict['data'].append(item)

				config['extra_series'].append(series_dict)

		return linegraph.plot(data, config)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			duplicates = (json[key]["duplicate"] / json[key]["totalFragmentsInput"]) * 100

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in": json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
							   "% Duplicates": duplicates,
							   "Notes": json[key]["Notes"],
							   "Saturation": json[key]["duplicate_saturation"]
						 	  }

		section = {
				   "Table": self.table(stats_json),
				   "Saturation Plot": self.linegraph(stats_json)
				   }

		return section

#################################################
class Primers():

	def table(self, json):

		headers = OrderedDict()

		headers["PE in"] = {'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["PE out"] = {'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["SE in"] = {'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["SE out"] = {'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Reads Flipped"] = {'description': 'Number of Flipped Reads', 'format': '{:,.0f}', 'scale': 'Oranges'}
		headers["Notes"] = {'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json):

		categories  = OrderedDict()

		categories['Primer 1 Only'] = {
									   'name': 'Primer 1 Only',
									   'color': '#4d8de4'
									  }
		categories['Primer 2 Only'] = {
									   'name': 'Primer 2 Only',
									   'color': '#e57433'
									  }
		categories['Both Primers'] = {
									   'name': 'Both Primers',
									   'color': '#33a02c'
									  }
		

		data  = OrderedDict()

		for sample in json.keys():

			data[sample] = {}

			for item in json[sample]["Primer Counts"]: 

				if item[0] != "None" and item[1] == "None":
					data[sample]["Primer 1 Only"] = item[2]

				elif item[0] == "None" and item[1] != "None":
					data[sample]["Primer 2 Only"] = item[2]

				elif item[0] != "None" and item[1] != "None" :
					data[sample]["Both Primers"] = item[2]


		return bargraph.plot(data, categories)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "PE in": json[key]["Paired_end"]["PE_in"],
							   "PE out": json[key]["Paired_end"]["PE_out"],
							   "SE in" : json[key]["Single_end"]["SE_in"],
							   "SE out": json[key]["Single_end"]["SE_out"],
							   "Reads Flipped": json[key]["Flipped"],
							   "Notes": json[key]["Notes"],
							   "Primers": json[key]["Primers"],
							   "Primer Counts": json[key]["Primers_counts"]
						 	  }

		section = {
				   "Table": self.table(stats_json),
				   "Reads with Primers": self.bargraph(stats_json)
				   }

		return section

#################################################
class Stats():

	def graph(self, json):

		histograms = ["R1 histogram", "R2 histogram", "SE histogram"]

		config = {'data_labels': [
								  {'name': "R1 histogram", 'ylab': 'Density', 'xlab': 'Read Lengths'},
								  {'name': "R2 histogram", 'ylab': 'Density', 'xlab': 'Read Lengths'},
								  {'name': "SE histogram", 'ylab': 'Density', 'xlab': 'Read Lengths'},
								 ],
				  'extra_series': [[], [], []]}

		data_list = []

		for i in range(len(histograms)):

			data = {}

			for key in json.keys():

				if len(data.keys()) == 0:
					data[key] = {}

					total = sum([ count[1] for count in json[key][histograms[i]] ])

					for item in json[key][histograms[i]]:

						data[key][item[0]] = item[1] / total

					data_list.append(data)

				else:
					series_dict = {
								   'name': key,
	        					   'data': []
	        					  }

					for item in json[key][histograms[i]]:
						temp = [item[0], (item[1] / total)]
						series_dict['data'].append(temp)

					config['extra_series'][i].append(series_dict)

		return linegraph.plot(data_list, config)


	def execute(self, json):

		stats_json = OrderedDict()

		for key in json.keys():

			stats_json[key] = {
			 				   "R1 histogram": json[key]["SE_readlength_histogram"],
			 				   "R2 histogram": json[key]["R1_readlength_histogram"],
			 				   "SE histogram": json[key]["R2_readlength_histogram"]
						 	  }

		section = {
				   "Density Plots": self.graph(stats_json)
				   }

		return section 

