from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################

class QWindowTrim():


	def table(self, json):

		# Standard table constructor. See MultiQC docs.
		headers = OrderedDict()

		headers["Qt_PE_in"] = {'title': "PE in", 'namespace': "PE in", 'description': 'Number of Input Paired End Reads', 'format': '{:,.0f}', 'scale': 'Greens' }
		headers["Qt_PE_out"] = {'title': "PE out", 'namespace': "PE out", 'description': 'Number of Output Paired End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Qt_SE_in"] = {'title': "SE in", 'namespace': "SE in",'description': 'Number of Input Single End Reads', 'format': '{:,.0f}', 'scale': 'Greens'}
		headers["Qt_SE_out"] = {'title': "SE out", 'namespace': "SE out", 'description': 'Number of Output Single End Reads', 'format': '{:,.0f}', 'scale': 'RdPu'}
		headers["Qt_Avg_BP_Trimmed"] = {'title': "Avg. BP Trimmed", 'namespace': "Avg. BP Trimmed", 'description': 'Average Number of Basepairs Trimmed per Read', 'format': '{:,.2f}', 'scale': 'Oranges'}
		headers["Qt_Notes"] = {'title': "Avg. BP Trimmed", 'namespace': "Notes", 'description': 'Notes'}

		return table.plot(json, headers)



	def bargraph(self, json, bps):

		# config dict for bar graph
		config = {
				  "title": "HTStream: QWindowTrim Trimmed Basepairs Bargraph",
				  'id': "htstream_qwindowtrimmer_bargraph",
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

			r1_data[key] = {"LT_R1": json[key]["Qt_Left_Trimmed_R1"],
						    "RT_R1": json[key]["Qt_Right_Trimmed_R1"]}

			r2_data[key] = {"LT_R2": json[key]["Qt_Left_Trimmed_R2"],
						    "RT_R2": json[key]["Qt_Right_Trimmed_R2"]}

			se_data[key] = {"LT_SE": json[key]["Qt_Left_Trimmed_SE"],
						    "RT_SE": json[key]["Qt_Right_Trimmed_SE"]}

		# returns nothing if no reads were trimmed.
		if bps == 0:
			html = '<div class="alert alert-info"> No basepairs were trimmed from any sample. </div>'	
			return html


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

		# accumular variable that prevents empty bar graphs
		total_trimmed_bp = 0

		for key in json.keys():

			# trimmed reads by side
			lefttrimmed_bp = json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"] + json[key]["Single_end"]["leftTrim"]
			rightrimmed_bp = json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"] + json[key]["Single_end"]["rightTrim"]

			# total trimmed reads
			trimmed_bp = (lefttrimmed_bp + rightrimmed_bp)

			# sample dictionary entry
			stats_json[key] = {
			 				   "Qt_PE_in": json[key]["Paired_end"]["in"],
							   "Qt_PE_out": json[key]["Paired_end"]["out"],
							   "Qt_SE_in" : json[key]["Single_end"]["in"],
							   "Qt_SE_out": json[key]["Single_end"]["out"],
							   "Qt_Avg_BP_Trimmed": trimmed_bp / json[key]["Fragment"]["in"],
							   "Qt_Notes": json[key]["Program_details"]["options"]["notes"],
							   "Qt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
							   "Qt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
							   "Qt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
							   "Qt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
							   "Qt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
							   "Qt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"]
						 	  }

			# total basepairs accumlation 
			total_trimmed_bp += trimmed_bp


		# sections and figure function calls
		section = {
				   "Table": self.table(stats_json),
				   "Trimmed Basepairs": self.bargraph(stats_json, total_trimmed_bp)
				   }

		return section
	