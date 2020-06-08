from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" AdapterTrimmer submodule for HTStream charts and graphs """

#################################################

class AdapterTrimmer():

	def table(self, json, total, zeroes):

		# Table constructor. Just like the MultiQC docs.

		if total == 0:
			return ""
		
		headers = OrderedDict()

		if zeroes == False:

			headers["At_%_BP_Lost"] = {'title': "% Bp Lost",
										'namespace': "% Bp Lost",
										'description': 'Percentage of Input bps (SE and PE) trimmed.',
										'suffix': '%',
										'format': '{:,.2f}',
										'scale': 'RdPu'
										}
			headers["At_%_Adapters"] = {'title': "% Adapters",
										'namespace': "% Adapters",
										'description': 'Percentage of Reads (SE and PE) with an Adapter',
										'suffix': '%',
										'format': '{:,.2f}',
										'scale': 'Blues'
										}

		else:

			headers["At_BP_Lost"] = {'title': "Bp Lost", 'namespace': "Bp Lost", 'description': 'Input bps (SE and PE) trimmed.', 'scale': 'RdPu', 'format': '{:,.0f}'}
			headers["At_Adapters"] = {'title': "Adapters", 'namespace': "Adapters", 'description': 'Reads (SE and PE) with an Adapter', 'scale': 'Blues', 'format': '{:,.0f}'}


		headers["At_Avg_BP_Trimmed"] = {'title': "Avg. Bps Trimmed", 'namespace': "Avg. Bps Trimmed", 'description': 'Average Number of basepairs trimmed from reads', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["At_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json, avg_bp_trimmed):

		# configuration dictionary for bar graph
		config = {'title': "HTStream: AdapterTrimmer Bp Composition Bargraph",
				  'id': "htstream_adaptertrimmer_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False}

		# if no overlaps at all are present, return nothing
		if avg_bp_trimmed == 0:
			html = '<div class="alert alert-info"> No adapters were trimmed from samples. </div>'	
			return html

		if len(json.keys()) > 150:
			html = '<div class="alert alert-info"> Too many samples for bargraph. </div>'	
			return html

		# bargraph dictionary. Exact use of example in MultiQC docs.
		categories  = OrderedDict()

		categories['At_R1'] = {
							  'name': 'Read 1',
							  'color': '#779BCC'
							 }
		categories['At_R2'] = {
							  'name': 'Read 2',
							  'color': '#C3C3C3'
							 }
		categories['At_SE'] = {
							  'name': 'Single End',
							  'color': '#D1ADC3'
							 }

		return bargraph.plot(json, categories, config)


	def execute(self, json):

		stats_json = OrderedDict()
		overview_dict = {}

		total = 0
		zeroes = False

		for key in json.keys():

			frag_in = json[key]["Fragment"]["in"]
			
			# calculations for reads with adapters and bps trimmed
			adapter_reads = json[key]["Single_end"]["adapterTrim"] + json[key]["Paired_end"]["Read1"]["adapterTrim"] + json[key]["Paired_end"]["Read2"]["adapterTrim"] # total reads trimmed
			bp_trimmed = json[key]["Single_end"]["adapterBpTrim"] + json[key]["Paired_end"]["Read1"]["adapterBpTrim"] + json[key]["Paired_end"]["Read2"]["adapterBpTrim"] # total basepairs trimmed
			perc_bp_lost = ( (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) / json[key]["Fragment"]["basepairs_in"] ) * 100

			# if adapter trim is zero, so is the percentage and the avg basepair trimmed. This prevents division by zero error
			if adapter_reads == 0:
				perc_adapters = 0
				avg_bp_trimmed = 0

			else:
				perc_adapters = (adapter_reads / frag_in) * 100
				avg_bp_trimmed = (bp_trimmed / adapter_reads)

			total += avg_bp_trimmed

			if perc_bp_lost < 0.01 and zeroes == False:
				zeroes = True


			overview_dict[key] = {
								  "Output_Bp": json[key]["Fragment"]["basepairs_out"],
								  "Bp_Lost": json[key]["Fragment"]["basepairs_out"] / json[key]["Fragment"]["basepairs_in"],
								  "R1_Bp_Trim": json[key]["Paired_end"]["Read1"]["adapterBpTrim"] / frag_in,
								  "R1_Read_Trim": json[key]["Paired_end"]["Read1"]["adapterTrim"] / frag_in,
								  "R2_Bp_Trim": json[key]["Paired_end"]["Read2"]["adapterBpTrim"] / frag_in,
								  "R2_Read_Trim": json[key]["Paired_end"]["Read2"]["adapterTrim"] / frag_in,
								  "SE_Bp_Trim": json[key]["Single_end"]["adapterBpTrim"] / frag_in,
								  "SE_Read_Trim": json[key]["Single_end"]["adapterTrim"] / frag_in
								  }

			# sample dictionary entry
			stats_json[key] = {
							   "At_%_BP_Lost": perc_bp_lost,
							   "At_%_Adapters": perc_adapters,
							   "At_BP_Lost": bp_trimmed,
							   "At_Adapters": adapter_reads,
							   "At_Avg_BP_Trimmed": avg_bp_trimmed,
							   "At_Notes": json[key]["Program_details"]["options"]["notes"],
							   "At_R1": json[key]["Paired_end"]["Read1"]["adapterBpTrim"],
							   "At_R2": json[key]["Paired_end"]["Read2"]["adapterBpTrim"],
							   "At_SE": json[key]["Single_end"]["adapterBpTrim"]
							  }


		# sections and figure function calls
		section = {"Table": self.table(stats_json, total, zeroes),
				   "Bp Composition Bargraph": self.bargraph(stats_json, total),
				   "Overview": overview_dict}

		return section