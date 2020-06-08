from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" PolyATTrim submodule for HTStream charts and graphs """

#################################################

class PolyATTrim():

	def table(self, json, bps, zeroes):

		# Table construction. Taken from MultiQC docs.

		if bps == 0:
			return ""

		headers = OrderedDict()

		if zeroes == False:
			headers["Pt_%_BP_Lost"] = {'title': "% Bp Lost", 'namespace': "% Bp Lost", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}
		else:
			headers["Pt_BP_Lost"] = {'title': "Total Bp Lost", 'namespace': "Total Bp Lost", 'description': 'Total input bps (SE and PE) trimmed.',
									 'format': '{:,.0f}', 'scale': 'Greens'}

		headers["Pt_%_R1_BP_Lost"] = {'title': "% Bp Lost from R1", 'namespace': "% Bp Lost from R1", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}
		headers["Pt_%_R2_BP_Lost"] = {'title': "% Bp Lost from R2", 'namespace': "% Bp Lost from R2", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'Greens'}
		headers["Pt_%_SE_BP_Lost"] = {'title': "% Bp Lost from SE", 'namespace': "% Bp Lost from SE", 'description': 'Percentage of Input bps (SE and PE) trimmed.',
									   'suffix': '%', 'format': '{:,.2f}', 'scale': 'RdPu'}


		if zeroes == False:
			headers["Pt_Avg_BP_Trimmed"] = {'title': "Avg. Bps Trimmed", 'namespace': "Avg. Bps Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Blues'}
			


		headers["Pt_Notes"] = {'title': "Notes", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)


	def execute(self, json):

		stats_json = OrderedDict()
		overview_dict = {}

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
				
			
			# number of trimmed reads by side
			lefttrimmed_bps = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_bps = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total number of trimmed reads.
			sample_trimmed_bps = (lefttrimmed_bps + rightrimmed_bps)

			if perc_bp_lost < 0.01 and zeroes == False:
				zeroes = True

			overview_dict[key] = {
								  "Output_Bp": json[key]["Fragment"]["basepairs_out"],
								  "Bp_Lost": json[key]["Fragment"]["basepairs_out"] / json[key]["Fragment"]["basepairs_in"]
								  }

			# sample entry in stats dictionary
			stats_json[key] = {
							   "Pt_%_BP_Lost": perc_bp_lost,
							   "Pt_BP_Lost": total_bp_lost,
							   "Pt_%_R1_BP_Lost": total_r1,
							   "Pt_%_R2_BP_Lost": total_r2,
							   "Pt_%_SE_BP_Lost": total_se,
							   "Pt_Avg_BP_Trimmed": total_bp_lost / json[key]["Fragment"]["in"],
							   "Pt_Notes": json[key]["Program_details"]["options"]["notes"],
							  }

			trimmed_bps += sample_trimmed_bps 

		# section and figure function calls
		section = {"Table": self.table(stats_json, trimmed_bps, zeroes),
				   "Overview": overview_dict}


		return section 