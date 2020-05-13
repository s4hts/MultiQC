from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" NTrimmer submodule for HTStream charts and graphs """

#################################################

class NTrimmer():


	def table(self, json, bps, zeroes):

		# Table construction. Taken from MultiQC docs.

		if bps == 0:
			return ""

		headers = OrderedDict()

		if zeroes == False:
			headers["Nt_%_BP_Lost"] = {'title': "% Bp Lost", 'namespace': "% Bp Lost", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}
		else:
			headers["Nt_BP_Lost"] = {'title': "Total Bp Lost", 'namespace': "Total Bp Lost", 'description': 'Total input bps (SE and PE) trimmed.',
									 'format': '{:,.0f}', 'scale': 'Greens'}

		headers["Nt_%_R1_BP_Lost"] = {'title': "% Bp Lost from R1", 'namespace': "% Bp Lost from R1", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Nt_%_R2_BP_Lost"] = {'title': "% Bp Lost from R2", 'namespace': "% Bp Lost from R2", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}
		headers["Nt_%_SE_BP_Lost"] = {'title': "% Bp Lost from SE", 'namespace': "% Bp Lost from SE", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}


		if zeroes == False:
			headers["Nt_Avg_BP_Trimmed"] = {'title': "Avg. Bps Trimmed", 'namespace': "Avg. Bps Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Blues'}
			

		headers["Nt_%_Discarded"] = {'title': "% Discarded",
									 'namespace': "% Discarded",
									 'description': 'Percentage of Reads (SE and PE) Discarded',
									 'suffix': '%',
									 'max': 100,
									 'format': '{:,.2f}',
									 'scale': 'Oranges'
									}

		headers["Nt_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, bps):

		# returns nothing if no reads were trimmed.
		if bps == 0:
			html = '<div class="alert alert-info"> No basepairs were trimmed from any sample. </div>'	
			return html

		# config dict for bar graph
		config = {
				  "title": "HTStream: NTrimmer Trimmed Basepairs Bargraph",
				  'id': "htstream_ntrimmer_bargraph",
				  'ylab' : "Samples",
				  'cpswitch_c_active': False,
				  'data_labels': [{'name': "Read 1"},
       							 {'name': "Read 2"},
       							 {'name': "Single End"}]
				  }

		html = ""

		r1_data = {}
		r2_data = {}
		se_data = {}

		for key in json:

			r1_data[key] = {"LT_R1": json[key]["Nt_Left_Trimmed_R1"],
						    "RT_R1": json[key]["Nt_Right_Trimmed_R1"]}

			r2_data[key] = {"LT_R2": json[key]["Nt_Left_Trimmed_R2"],
						    "RT_R2": json[key]["Nt_Right_Trimmed_R2"]}

			se_data[key] = {"LT_SE": json[key]["Nt_Left_Trimmed_SE"],
						    "RT_SE": json[key]["Nt_Right_Trimmed_SE"]}



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

		# accumulator variable. Used to prevent empty bargraphs 
		trimmed_bps = 0
		zeroes = False

		for key in json.keys():

			total_bp_lost = (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) 

			if total_bp_lost == 0:
				perc_bp_lost = 0
				total_r1 = 0 
				total_r2 = 0
				total_se = 0 

			else:
				perc_bp_lost = ( total_bp_lost / json[key]["Fragment"]["basepairs_in"] ) * 100

				total_r1 = ( (json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]) / total_bp_lost ) * 100
				total_r2 = ( (json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]) / total_bp_lost) * 100
				total_se = ( (json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]) / total_bp_lost ) * 100
				

			# number ofreads discarded
			discarded_reads = json[key]["Single_end"]["discarded"] + json[key]["Paired_end"]["discarded"] 
			
			# number of trimmed reads by side
			lefttrimmed_bps = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_bps = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total number of trimmed reads.
			sample_trimmed_bps = (lefttrimmed_bps + rightrimmed_bps)

			if perc_bp_lost < 0.01 and zeroes == False:
				zeroes = True

			# sample entry in stats dictionary
			stats_json[key] = {
							   "Nt_%_BP_Lost": perc_bp_lost,
							   "Nt_BP_Lost": total_bp_lost,
							   "Nt_%_R1_BP_Lost": total_r1,
							   "Nt_%_R2_BP_Lost": total_r2,
							   "Nt_%_SE_BP_Lost": total_se,
							   "Nt_Avg_BP_Trimmed": total_bp_lost / json[key]["Fragment"]["in"],
							   "Nt_%_Discarded" : (discarded_reads  / json[key]["Fragment"]["in"]) * 100,
							   "Nt_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Nt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
							   "Nt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
							   "Nt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
							   "Nt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
							   "Nt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
							   "Nt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"]
							  }

			trimmed_bps += sample_trimmed_bps 

		# section and figure function calls
		section = {
				   "Table": self.table(stats_json, trimmed_bps, zeroes),
				   "Trimmed Reads": self.bargraph(stats_json, trimmed_bps)
				   }

		return section
