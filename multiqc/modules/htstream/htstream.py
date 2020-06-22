#!/usr/bin/env python

""" MultiQC module to parse output from HTStream """

from __future__ import print_function
from collections import OrderedDict
import logging
import re, json, os, operator

# HTStream Apps
from .apps import AdapterTrimmer, CutTrim, LengthFilter, Overlapper, QWindowTrim, NTrimmer
from .apps import PolyATTrim, SeqScreener, SuperDeduper, Primers, Stats, OverviewStats, htstream_utils

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule


#################################################

# Logger Initialization
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

	def __init__(self):

		self.sample_statistics = {}

		# Initialise the parent object
		super(MultiqcModule, self).__init__(name='HTStream',
		anchor='htstream', href='https://s4hts.github.io/HTStream/',
		info=" quality control and processing pipeline for High Throughput Sequencing data ")


		# Initialize ordered dictionary (key: samples, values: their respective json files)
		self.data = OrderedDict()

		# Import js and css functions.
		self.js = { 'assets/js/htstream.js' : os.path.join(os.path.dirname(__file__), 'assets', 'js', 'htstream.js') }
		self.css = { 'assets/css/htstream.css' : os.path.join(os.path.dirname(__file__), 'assets', 'css', 'htstream.css') }

		 # iterates through files found by "find_log_files" (located in base_module.py, re patterns found in search_patterns.yml)
		for file in self.find_log_files('htstream'):

			self.add_data_source(file) # write file to MultiQC souce file 

			self.s_name = file['s_name'] # sample name
			self.file_data = self.parse_json(file['s_name'], file['f']) # parse stats file. Should return json directory of apps and their stats 

			self.data[self.s_name] = self.file_data # add sample and stats to OrderedDict


		# make sure samples are being processed 
		if len(self.data) == 0:
			raise UserWarning

		# remove excluded samples 
		self.data = self.ignore_samples(self.data)

		# parse json containing stats on each sample
		self.generate_reports(self.data) 

		# general stats table, can't upload dictionary of dictionaries :/
		#self.general_stats_addcols(self.data)


	#################################################
	# Json and stats parsing functions

	def parse_json(self, name, f):

		app_dict = {}
		apps = json.loads(f)

		try:

			for a in apps:
				i = 1
				app_name = a["Program_details"]["program"] + "_" + str(i)

				while app_name in app_dict.keys():
					i += 1
					app_name = a["Program_details"]["program"] + "_" + str(i)

				app_dict[app_name] = a	

		except:

			app_dict = json.loads(f, object_pairs_hook=htstream_utils.resolve)
			log.warning("Sample " + name + " uses old json format. Please update to a newer version of HTStream.")


		return app_dict


	def generate_reports(self, json):

		# preserves order of apps in report
		self.programs = {
						'AdapterTrimmer': {"app": AdapterTrimmer.AdapterTrimmer(),
										   "description": "Trims adapters which are sequenced when the fragment insert length is shorter than the read length."},

						'CutTrim': {"app": CutTrim.CutTrim(),
								    "description": "Trims a fixed number of bases from the 5' and/or 3' end of each read."},

						'LengthFilter': {"app": LengthFilter.LengthFilter(),
								    "description": "Discards reads below a minimum length threshold."},

						'NTrimmer': {"app": NTrimmer.NTrimmer(),
									 "description": "Trims reads to the longest subsequence that contains no Ns."},
						
						'Overlapper': {"app": Overlapper.Overlapper(),
									   "description": "Attempts to overlap paired end reads to produce the original fragment, trims adapters, and can correct sequencing errors."},

						'PolyATTrim': {"app": PolyATTrim.PolyATTrim(),
									   "description": "Attempts to trim poly-A and poly-T sequences from the end of reads."},

						'Primers': {"app": Primers.Primers(),
									"description": "Identifies primer sequences located on the 5' ends of R1 and R2, or 5' and 3' end of SE reads."},

						'QWindowTrim': {"app": QWindowTrim.QWindowTrim(),
										"description": "Uses a sliding window approach to remove the low quality ends of reads."},

						'SeqScreener': {"app": SeqScreener.SeqScreener(),
										"description": "A simple sequence screening tool which uses a kmer lookup approach to identify reads from an unwanted source."},

						'SuperDeduper': {"app": SuperDeduper.SuperDeduper(),
										 "description": "A reference free duplicate read removal tool."},

						'Stats': {"app": Stats.Stats(),
								  "description": "Generates a JSON formatted file containing a set of statistical measures about the input read data."}
						}


		self.report_sections = {}

		# checks that order is consistent within stats files 
		app_order = []
		stats_section = ""

		for key in json.keys():

			temp = list(json[key].keys())
			
			if app_order == []:
				app_order = temp

			elif app_order == temp:
				continue

			else:
				log.error("Inconsistent order of HTStream applications.")



		# scold people that don't read the documentation
		if "hts_Stats_1" not in app_order:
			log.warning("hts_Stats not found. It is recommended you run this app before and after pipeline.")
		
		self.overview_stats = {"Pipeline Input": {}}

		# sort list of samples
		sample_keys = list(sorted(json.keys()))
		excludes = []
		stats_wrapper = False
		pipeline_input = True


		############################
		# GENERATE REPORT SECTIONS 

		for i in range(len(app_order)):

			app = app_order[i]
			program = app.split("hts_")[-1].split("_")[0]

			if program not in self.programs.keys():
				log.warning("hts_" + program + " is currently not supported by MultiQC: HTStrean. Apps currently supported: " + htstream_utils.key_print(self.programs))
				continue

			# creat app specific dictionary, each entry will be a sample
			stats_dict = OrderedDict()


			for key in sample_keys:

				stats_dict[key] = json[key][app]

				if pipeline_input == True:

					self.overview_stats["Pipeline Input"][key] = {
																 "Input_Reads": json[key][app]["Fragment"]["in"],
																 "Input_Bp": json[key][app]["Fragment"]["basepairs_in"]
																 }
					

			pipeline_input = False					

			# if data exists for app, execute app specific stats processing
			if len(stats_dict.keys()) != 0:

				app_name = app 
				app = program
				index = app_name.split("_")[-1]

				# dictionary of subsections
				section_dict = self.programs[app]["app"].execute(stats_dict, index)

				# if dictionary is not empty
				if len(section_dict.keys()) != 0:

					self.overview_stats[app_name] = section_dict["Overview"] 

					# construct html for section
					html = ""
					for title, section in section_dict.items():

						if section != "" and title != "Overview":
							html += section + '<br>\n'

					# remove trailing space
					html = html[:-5]


					# add description for app 
					description = self.programs[app]["description"]

					self.report_sections[app_name] = {'description': description,
													  'html': html}



		############################
		# ADD SECTIONS TO FILE

		# add pipeline overview section if appropriate
		if self.overview_stats != {}:

			try:
				app = OverviewStats.OverviewStats()

				description = "General statistics from the HTStream pipeline."
				html, pca_data = app.execute(self.overview_stats, app_order)
					
				self.write_data_file(pca_data, 'htstream_raw_pca_data_matrix')
				
				self.add_section(name = "Processing Overview",
								 description = description,
								 content = html) 

			except:
			 	log.warning("Report Section for Processing Overview Failed.")




		# add app sections
		for section, content in self.report_sections.items():

				temp_list = section.split("_")
				name = temp_list[0] + "_" + temp_list[1] + " " + temp_list[-1]

				try:
					self.add_section(name = name,
									 description = content["description"],
									 content = content["html"])


				except:
					msg = "Report Section for " + section + " Failed."
					log.warning(msg)




