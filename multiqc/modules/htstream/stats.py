#!/usr/bin/env python

""" MultiQC submodule for HTStream charts and graphs """

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph, linegraph


#################################################
# HTStream App Classes

class  AdapterTrimmer():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot  

class CutTrim():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 


class Overlapper():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 


class QWindowTrim():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 

class NTrimmer():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 


class PolyATTrim():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 


class SeqScreener():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 


class SuperDeduper():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 

class Primers():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 

class Stats():

	def execute(self, json):

		categories = ["totalFragmentsInput", "totalFragmentsOutput"]
		plot = bargraph.plot(json, categories)
		return plot 

