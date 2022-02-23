from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import linegraph

#################################################

""" SuperDeduper submodule for HTStream charts and graphs """

#################################################


class SuperDeduper:

    ########################
    # Info about App
    def __init__(self):
        self.info = "A reference free duplicate read removal tool."
        self.type = "read_reducer"

    ########################
    # Linegraph Function
    def linegraph(self, json, index):

        # plot configurations, list of options in MultiQC docs
        config = {
            "id": "htstream_superdedup_" + index,
            "title": "HTStream: Duplicates",
            "xlab": "Percentage of Reads",
            "ylab": "Percent Duplicates",
            "data_labels": [{"name": "Percent Duplicates", 
                                "xlab": "Percentage of Reads",
                                "ylab": "Percent Duplicates"},
                            {"name": "Duplicate Saturation", 
                                "xlab": "Total Reads",
                                "ylab": "Unique Reads"}]
        }

        # initialize data structures and variabe;s
        data = [{}, {}]

        for key in json.keys():

            data[0][key] = {}
            data[1][key] = {}

            for item in json[key]["Sd_Saturation"]:

                perc = (item[0] / json[key]["Sd_Total_Reads"]) * 100 
                perc_dup = (item[1] / json[key]["Sd_Total_Reads"]) * 100

                data[0][key][perc] = perc_dup
                data[1][key][item[0]] = item[0] - item[1]

        return linegraph.plot(data, config)

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        for key in json.keys():

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Percent_Duplicates": json[key]["Fragment"]["duplicate"] / json[key]["Fragment"]["in"],
                "Percent_Ignored": json[key]["Fragment"]["ignored"] / json[key]["Fragment"]["in"],
            }

            # sample instance in ordered dict
            stats_json[key] = {
                "Sd_Total_Reads": json[key]["Fragment"]["in"],
                "Sd_Duplicates": json[key]["Fragment"]["duplicate"],
                "Sd_Saturation": json[key]["Fragment"]["duplicate_saturation"],
            }

        # output dictionary, keys are section, value is function called for figure generation
        section = {
            "Duplicate Saturation": self.linegraph(stats_json, index),
            "Overview": overview_dict,
        }

        return section
