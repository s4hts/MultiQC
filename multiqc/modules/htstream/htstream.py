#!/usr/bin/env python

""" MultiQC module to parse output from HTStream """

from __future__ import print_function
from collections import OrderedDict
import logging
import re, json, os, operator

# HTStream Apps
from .apps import AdapterTrimmer, CutTrim, Overlapper, QWindowTrim, NTrimmer, PolyATTrim
from .apps import  SeqScreener, SuperDeduper, Primers, Stats, OverviewStats, htstream_utils

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
		anchor='htstream', href='https://ibest.github.io/HTStream/',
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
			self.file_data = self.parse_json(file['f']) # parse stats file. Should return json directory of apps and their stats 

			self.data[self.s_name] = self.file_data # add sample and stats to OrderedDict


		# make sure samples are being processed 
		if len(self.data) == 0:
			raise UserWarning

		# remove excluded samples 
		self.data = self.ignore_samples(self.data)

		# parse json containing stats on each sample
		self.parse_stats(self.data) 

		# general stats table, can't upload dictionary of dictionaries :/
		#self.general_stats_addcols(self.data)


	#################################################
	# Json and stats parsing functions

	def parse_json(self, f):

		return json.loads(f, object_pairs_hook=htstream_utils.resolve)


	def parse_stats(self, json):

		# preserves order of apps in report
		self.programs = {
						'AdapterTrimmer': {"app": AdapterTrimmer.AdapterTrimmer(),
										   "description": "Trims adapters which are sequenced when the fragment insert length is shorter than the read length."},

						'CutTrim': {"app": CutTrim.CutTrim(),
								    "description": "Trims a fixed number of bases from the 5' and/or 3' end of each read."},

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
		self.summary_stats = {}


		# checks that order is consistent within stats files 
		app_order = []	
		stats_section = ""			
		for key in json.keys():

			sample = key
			if app_order == []:
				app_order = list(json[key].keys())

			elif app_order == list(json[key].keys()):
				app_order = list(json[key].keys())

			else:
				log.error("Inconsistent order of HTStream applications.")


		# scold people that don't read the documentation
		if "hts_Stats" not in app_order:
			log.warning("hts_Stats not found. It is recommended you run this app before and after pipeline.")
			self.overview_stats = {}

		else:
			if len(json[key]["hts_Stats"]) == 2:
				self.overview_stats = {"Pipeline Input": {},
											"hts_Stats": {}}
			else:
				self.overview_stats = {"hts_Stats": {}}


		# sort list of samples
		sample_keys = list(sorted(json.keys()))
		excludes = []
		stats_wrapper = False


		for i in range(len(app_order)):

			app = app_order[i]
			program = app.split("hts_")[-1]

			if program not in self.programs.keys():
				log.warning(app + " is currently not supported by MultiQC: HTStrean.")
				excludes.append(app)
				continue

			
			# creat app specific dictionary, each entry will be a sample
			stats_dict = OrderedDict()

			for key in sample_keys:

				stats_dict[key] = json[key][app]

				if app == "hts_Stats":

					if len(json[key][app]) == 2:

						stats_wrapper = True

						self.overview_stats["Pipeline Input"][key] = {
																	 "Fragment_Section": json[key][app][0]["Fragment"],
																	 "total_Q30": json[key][app][0]["Paired_end"]["Read1"]["total_Q30_basepairs"] + json[key][app][0]["Paired_end"]["Read2"]["total_Q30_basepairs"],
																	 "Read_Breakdown":{
																	 				   "Paired_end": json[key][app][0]["Paired_end"]["in"]
																	 				   },
																	 "Input_Reads": json[key][app][0]["Fragment"]["in"],
																	 "Input_Bp": json[key][app][0]["Fragment"]["basepairs_in"]
																	}



						try:
							self.overview_stats["Pipeline Input"][key]["total_Q30"] += json[key][app][0]["Single_end"]["total_Q30_basepairs"]
							self.overview_stats["Pipeline Input"][key]["Read_Breakdown"]["Single_end"] = json[key][app][0]["Single_end"]["in"]

						except:
							pass

					self.overview_stats[app][key] = {
													 "Fragment_Section": json[key][app][-1]["Fragment"],
													 "total_Q30": json[key][app][-1]["Paired_end"]["Read1"]["total_Q30_basepairs"] + json[key][app][-1]["Paired_end"]["Read2"]["total_Q30_basepairs"],
													 "Read_Breakdown": {
													 					"Paired_end": json[key][app][-1]["Paired_end"]["in"]
													 					},
													 "Input_Reads": json[key][app][-1]["Fragment"]["in"],
													 "Input_Bp": json[key][app][-1]["Fragment"]["basepairs_in"]
													}


					try:
						self.overview_stats[app][key]["total_Q30"] += json[key][app][-1]["Single_end"]["total_Q30_basepairs"]
						self.overview_stats[app][key]["Read_Breakdown"]["Single_end"] = json[key][app][-1]["Single_end"]["in"]

					except:
							pass	
					

			# if data exists for app, execute app specific stats processing
			if len(stats_dict.keys()) != 0:

				app = app.split("hts_")[-1]

				# dictionary of subsections
				section_dict = self.programs[app]["app"].execute(stats_dict)

				# if dictionary is not empty
				if len(section_dict.keys()) != 0:

					if app != "Stats":
						self.overview_stats["hts_" + app] = section_dict["Overview"] 

					html = ""

					for title, section in section_dict.items():

						if section != "" and title != "Overview":
							html += section + '<br>\n'

					# remove trailing space
					html = html[:-5]

					description = self.programs[app]["description"]

					if app == "Stats":

						self.summary_stats = {'description': description,
											  'html': html}

					else:

						self.report_sections[app] = {'description': description,
													 'html': html}


		# create logical order for general overview table
		overview_order = []
		
		if stats_wrapper == True:
			overview_order.append("Pipeline Input")

			for a in app_order:
				if a != "hts_Stats" and a not in excludes:
					overview_order.append(a)

			overview_order.append("hts_Stats")

		else:
			for a in app_order:
				if a not in excludes:
					overview_order.append(a)


		# add pipeline overview 

		if self.overview_stats != {}:

			app = OverviewStats.OverviewStats()
			description = "General statistics from the HTStream pipeline."
			html = app.execute(self.overview_stats, overview_order)
			
			self.add_section(name = "Processing Overview",
							 description = description,
							 content = html) 


		# add apps

		for section, content in self.report_sections.items():

				try:
					self.add_section(name = section,
									 description = content["description"],
									 content = content["html"])


				except:
					msg = "Report Section for " + section + " Failed."
					log.warning(msg)


		# add summary stats

		if self.summary_stats != {}:

			try:
				self.add_section(name = "Summary Stats",
								 description = self.summary_stats["description"],
								 content = self.summary_stats["html"])


			except:

				log.warning("Report Section for hts_Stats Failed.")


