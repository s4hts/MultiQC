from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" CutTrim submodule for HTStream charts and graphs """

#################################################


class CutTrim:

    ########################
    # Info about App
    def __init__(self):
        self.info = "Trims a fixed number of bases from the 5' and/or 3' end of each read."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def table(self, json, index):

        # Table constructor. Just like the MultiQC docs.
        headers = OrderedDict()

        headers["Ct_%_BP_Lost" + index] = {
            "title": "% Bp Lost",
            "namespace": "% Bp Lost",
            "description": "Percentage of Input bps (SE and PE) trimmed.",
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "RdPu",
        }

        return table.plot(json, headers)

    ########################
    # MainFunction
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}
        overall_pe = 0
        overall_se = 0

        for key in json.keys():

            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            perc_bp_lost = (total_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100
            total_r1 = (
                json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]
            )
            total_r2 = (
                json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]
            )
            total_se = json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": (total_bp_lost / json[key]["Fragment"]["basepairs_in"]),
            }

            # sample dictionary entry
            stats_json[key] = {
                "Ct_%_BP_Lost" + index: perc_bp_lost,
                "Ct_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
                "Ct_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
                "Ct_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
                "Ct_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
                "Ct_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
                "Ct_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"],
            }

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, index),
            "Overview": overview_dict,
        }

        return section
