from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################

class QWindowTrim():

	########################
	# Info about App
	def __init__(self):
		self.info = "Uses a sliding window approach to remove the low quality ends of reads."
		self.type = "bp_reducer"		
	

	########################
	# Table Function
	def table(self, json, pe_bps, se_bps, index):

		# Table construction. Taken from MultiQC docs.

		# If no SE and BP lost, dont add table
		if (pe_bps + se_bps) == 0:
			return ""

		headers = OrderedDict()

		headers["Qt_%_BP_Lost" + index] = {'title': "% Bp Lost", 'namespace': "% Bp Lost", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
								   'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}
		
		# If PE data, add cols
		if pe_bps != 0:
			headers["Qt_%_R1_BP_Lost" + index] = {'title': "% R1 of Bp Lost", 'namespace': "% Bp Lost from R1", 'description': 'Percentage of total trimmed bps.',
										  'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
			headers["Qt_%_R2_BP_Lost" + index] = {'title': "% R2 of Bp Lost", 'namespace': "% Bp Lost from R2", 'description': 'Percentage of total trimmed bps.',
										  'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}

		# If SE data, add cols
		if se_bps != 0:
			headers["Qt_%_SE_BP_Lost" + index] = {'title': "% SE of Bp Lost", 'namespace': "% Bp Lost from SE", 'description': 'Percentage of total trimmed bps.',
										  'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}

		headers["Qt_Avg_BP_Trimmed" + index] = {'title': "Avg. Bps Trimmed", 'namespace': "Avg. Bpss Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 
										'format': '{:,.2f}', 'scale': 'Blues'}
		headers["Qt_%_Discarded" + index] = {'title': "% Discarded",
									 'namespace': "% Discarded",
									 'description': 'Percentage of Reads (SE and PE) Discarded',
									 'suffix': '%',
									 'format': '{:,.2f}',
									 'scale': 'Oranges'
									}

		headers["Qt_Notes" + index] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)


	########################
	# Bargraph Function
	def bargraph(self, json, bps, index):

		# config dict for bar graph
		config = {
				  "title": "HTStream: QWindowTrim Trimmed Basepairs Bargraph",
				  'id': "htstream_qwindowtrimmer_bargraph_" + index,
				  'ylab' : "Basepairs",
				  'cpswitch_c_active': False,
				  'data_labels': [{'name': "Read 1"},
       							 {'name': "Read 2"},
       							 {'name': "Single End"}]
				  }

		# Header 
		html = "<h4> QWindowTrim: Trimmed Basepairs Composition </h4>\n" 
		html += "<p>Plots the number of low quality basepairs trimmed from ends of paired end and single end reads.</p>"


		# If too many samples, don't add bargraph
		if len(json.keys()) > 150:
			html += '<div class="alert alert-warning"> <strong>Warning:</strong> Too many samples for bargraph. </div>'	
			return html
			
		r1_data = {}
		r2_data = {}
		se_data = {}

		# Construct data for multidataset bargraph
		for key in json:

			r1_data[key] = {"LT_R1": json[key]["Qt_Left_Trimmed_R1"],
						    "RT_R1": json[key]["Qt_Right_Trimmed_R1"]}

			r2_data[key] = {"LT_R2": json[key]["Qt_Left_Trimmed_R2"],
						    "RT_R2": json[key]["Qt_Right_Trimmed_R2"]}

			se_data[key] = {"LT_SE": json[key]["Qt_Left_Trimmed_SE"],
						    "RT_SE": json[key]["Qt_Right_Trimmed_SE"]}

		# returns nothing if no reads were trimmed.
		if bps == 0:
			html += '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'	
			return html

		# Create categories for multidataset bargraph
		cats = [OrderedDict(), OrderedDict(), OrderedDict()]
		cats[0]["LT_R1"] =   {'name': 'Left Trimmmed'}
		cats[0]["RT_R1"] =  {'name': 'Right Trimmmed'}
		cats[1]["LT_R2"] =   {'name': 'Left Trimmmed'}
		cats[1]["RT_R2"] =  {'name': 'Right Trimmmed'}
		cats[2]["LT_SE"] =   {'name': 'Left Trimmmed'}
		cats[2]["RT_SE"] =  {'name': 'Right Trimmmed'}

		# Create bargraph
		html += bargraph.plot([r1_data, r2_data, se_data], cats, config)

		return html


	########################
	# Main Function
	def execute(self, json, index):

		stats_json = OrderedDict()
		overview_dict = {}

		# accumular variable that prevents empty bar graphs
		overall_se = 0
		overall_pe = 0

		for key in json.keys():

			total_bp_lost = (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) 

			# If no bps lost, prevent zero division
			if total_bp_lost == 0:
				perc_bp_lost = 0
				total_r1 = 0 
				total_r2 = 0
				total_se = 0 
				total_pe = 0

			else:
				perc_bp_lost = (total_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100

				# Will fail if no PE data
				try:
					total_r1 = ( (json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]) / total_bp_lost ) * 100
					total_r2 = ( (json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]) / total_bp_lost) * 100
					left_pe_trimmed = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"]
					right_pe_trimmed = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"]
					total_pe = left_pe_trimmed + right_pe_trimmed

				except:
					total_r1 = 0
					total_r2 = 0
					total_pe = 0
					
				# Will fail if no SE data
				try:
					total_se = ( (json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]) / total_bp_lost ) * 100
					left_se_trimmed = json[key]["Single_end"]["leftTrim"]
					right_se_trimmed = json[key]["Single_end"]["rightTrim"]

				except:
					total_se = 0


			bp_in = json[key]["Fragment"]["basepairs_in"]

			# overview data
			overview_dict[key] = {
								  "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"] + json[key]["Paired_end"]["Read2"]["basepairs_out"],
								  "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
								  "Fraction_R1_Bp_Trimmed_Left": json[key]["Paired_end"]["Read1"]["leftTrim"] / bp_in, 
								  "Fraction_R1_Bp_Trimmed_Right": json[key]["Paired_end"]["Read1"]["rightTrim"] / bp_in, 
								  "Fraction_R2_Bp_Trimmed_Left": json[key]["Paired_end"]["Read2"]["leftTrim"] / bp_in, 
								  "Fraction_R2_Bp_Trimmed_Right": json[key]["Paired_end"]["Read2"]["rightTrim"] / bp_in, 
								  "Fraction_SE_Bp_Trimmed_Left": json[key]["Single_end"]["leftTrim"] / bp_in, 
								  "Fraction_SE_Bp_Trimmed_Right": json[key]["Single_end"]["rightTrim"] / bp_in, 
								  }

			# sample dictionary entry
			stats_json[key] = {
							   "Qt_%_BP_Lost" + index: perc_bp_lost,
							   "Qt_%_R1_BP_Lost" + index: total_r1,
							   "Qt_%_R2_BP_Lost" + index: total_r2,
							   "Qt_%_SE_BP_Lost" + index: total_se,
							   "Qt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
							   "Qt_Notes" + index: json[key]["Program_details"]["options"]["notes"],
							   "Qt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
							   "Qt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
							   "Qt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
							   "Qt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
							   "Qt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
							   "Qt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"]
						 	  }

			# total basepairs accumlation 
			overall_pe += total_pe
			overall_se += total_se


		# sections and figure function calls
		section = {"Table": self.table(stats_json, overall_pe, overall_se, index),
				   "Trimmed Basepairs": self.bargraph(stats_json, (overall_pe + overall_se), index),
				   "Overview": overview_dict}

		return section
	