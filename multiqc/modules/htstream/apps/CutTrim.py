from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" CutTrim submodule for HTStream charts and graphs """

#################################################

class CutTrim():

	def table(self, json, bps):

		# returns nothing if no reads were trimmed.
		if bps == 0:
			return ""


		# Table constructor. Just like the MultiQC docs.
		
		headers = OrderedDict()

		headers["Ct_%_BP_Lost"] = {'title': "% Bp Lost", 'namespace': "% Bp Lost", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
								   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Ct_%_R1_BP_Lost"] = {'title': "% Bp Lost from R1", 'namespace': "% Bp Lost from R1", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
								   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Ct_%_R2_BP_Lost"] = {'title': "% Bp Lost from R2", 'namespace': "% Bp Lost from R2", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
								   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Ct_%_SE_BP_Lost"] = {'title': "% Bp Lost from SE", 'namespace': "% Bp Lost from SE", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
								   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Ct_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)

	def bargraph(self, json, bps):

		# config dict for bar graph
		config = {
				  "title": "HTStream: CutTrim Trimmed Basepairs Bargraph",
				  'id': "htstream_qwindowtrimmer_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False,
				  'data_labels': [{'name': "Read 1"},
       							 {'name': "Read 2"},
       							 {'name': "Single End"}]
				  }
				  
		# returns nothing if no reads were trimmed.
		if bps == 0:
			html = '<div class="alert alert-info"> No basepairs were trimmed from any sample. </div>'	
			return html

		if len(json.keys()) > 150:
			html = '<div class="alert alert-info"> Too many samples for bargraph. </div>'	
			return html

		html = ""

		r1_data = {}
		r2_data = {}
		se_data = {}

		for key in json:

			r1_data[key] = {"LT_R1": json[key]["Ct_Left_Trimmed_R1"],
						    "RT_R1": json[key]["Ct_Right_Trimmed_R1"]}

			r2_data[key] = {"LT_R2": json[key]["Ct_Left_Trimmed_R2"],
						    "RT_R2": json[key]["Ct_Right_Trimmed_R2"]}

			se_data[key] = {"LT_SE": json[key]["Ct_Left_Trimmed_SE"],
						    "RT_SE": json[key]["Ct_Right_Trimmed_SE"]}


		cats = [OrderedDict(), OrderedDict(), OrderedDict()]
		cats[0]["LT_R1"] =   {'name': 'Left Trimmmed'}
		cats[0]["RT_R1"] =  {'name': 'Right Trimmmed'}
		cats[1]["LT_R2"] =   {'name': 'Left Trimmmed'}
		cats[1]["RT_R2"] =  {'name': 'Right Trimmmed'}
		cats[2]["LT_SE"] =   {'name': 'Left Trimmmed'}
		cats[2]["RT_SE"] =  {'name': 'Right Trimmmed'}


		return bargraph.plot([r1_data, r2_data, se_data], cats, config)



	def execute(self, json):

		stats_json = OrderedDict()
		overview_dict = {}

		total = 0

		for key in json.keys():
			
			total_bp_lost = (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) 

			if total_bp_lost == 0:
				perc_bp_lost  = 0
				total_r1 = 0 
				total_r2 = 0
				total_se = 0 

			else:
				perc_bp_lost = ( (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) / json[key]["Fragment"]["basepairs_in"] ) * 100

				total_r1 = ( (json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]) / total_bp_lost ) * 100
				total_r2 = ( (json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]) / total_bp_lost) * 100
				total_se = ( (json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]) / total_bp_lost ) * 100
				
			total += total_bp_lost

			overview_dict[key] = {
								  "Bp_Lost": json[key]["Fragment"]["basepairs_out"] / json[key]["Fragment"]["basepairs_in"] 
								 }

			# sample dictionary entry
			stats_json[key] = {
							   "Ct_%_BP_Lost": perc_bp_lost,
							   "Ct_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Ct_%_R1_BP_Lost": total_r1,
							   "Ct_%_R2_BP_Lost": total_r2,
							   "Ct_%_SE_BP_Lost": total_se,
							   "Ct_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"], 
							   "Ct_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
							   "Ct_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"], 
							   "Ct_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
							   "Ct_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
							   "Ct_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"]
							  }

		# sections and figure function calls
		section = {"Table": self.table(stats_json, total),
				   "Trimmed Bp Composition Bargraph": self.bargraph(stats_json, total),
				   "Overview": overview_dict}


		return section
