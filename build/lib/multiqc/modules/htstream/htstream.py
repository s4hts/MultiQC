#!/usr/bin/env python

""" MultiQC module to parse output from HTStream """

from __future__ import print_function
from collections import OrderedDict
import logging
import re, json

from .apps import AdapterTrimmer, CutTrim, Overlapper, QWindowTrim, NTrimmer
from .apps import PolyATTrim, SeqScreener, SuperDeduper, Primers, Stats
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

		 # iterates through files found by "find_log_files" (located in base_module.py, re patterns found in search_patterns.yml)
		for file in self.find_log_files('htstream'):

			self.add_data_source(file) # write file to MultiQC souce file 

			self.s_name = file['s_name'] # sample name
			self.file_data = self.parse_json(file['f']) # parse stats file. Should return json directory of apps and their stats 

			self.data[self.s_name] = self.file_data # add sample and stats to OrderedDict

			
		# remove excluded samples 
		self.data = self.ignore_samples(self.data)

		# make sure samples are being processed 
		if len(self.data) == 0:
			raise UserWarning

		# parse json containing stats on each sample
		self.parse_stats(self.data) 

		# general stats table, can't upload dictionary of dictionaries :/
		#self.general_stats_addcols(self.data)


	#################################################
	# Json and stats parsing functions

	def parse_json(self, f):

		return json.loads(f)


	def parse_stats(self, json):

		self.apps = {
            'AdapterTrimmer': AdapterTrimmer.AdapterTrimmer(),
            'CutTrim': CutTrim.CutTrim(),
            'Overlapper': Overlapper.Overlapper(),
            'QWindowTrim': QWindowTrim.QWindowTrim(),
            'NTrimmer': NTrimmer.NTrimmer(),
            'PolyATTrim': PolyATTrim.PolyATTrim(),
            'SeqScreener': SeqScreener.SeqScreener(),
            'SuperDeduper': SuperDeduper.SuperDeduper(),
            'Primers': Primers.Primers(),
            'Stats': Stats.Stats()
            }

		for app in self.apps.keys():

			stats_dict = OrderedDict()

			for key in json.keys():

				stats_boolean = False

				for subkey in json[key].keys():	

					if app == 'Stats' and app in subkey:
							if stats_boolean == False:
								stats_boolean = True
							else:
								stats_dict[key] = json[key][subkey]
								break

					elif app in subkey:
						stats_dict[key] = json[key][subkey]
						break

			if len(stats_dict.keys()) != 0:

				section_dict = self.apps[app].execute(stats_dict)

				if section_dict != None:

					section = "hts_" + app

					for key, value in section_dict.items():

						try:

							self.add_section(name = str(section + ": " + key),
											 plot = section_dict[key])

						except:
							pass


				# Possibly will be of use when more is known about what to include 

				# self.add_section(name = 'HTStream',
				# 				 anchor = section,
				# 				 description = 'This plot shows some really nice data.',
				# 				 helptext = 'This longer string (can be **markdown**) helps explain how to interpret the plot',
				# 				 plot = plot
				# 				 )



