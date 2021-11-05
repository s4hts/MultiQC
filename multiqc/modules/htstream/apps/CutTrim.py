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

        headers["Ct_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)

    ########################
    # Bargraph Function
    def bargraph(self, json, index):

        # config dict for bar graph
        config = {
            "title": "HTStream: CutTrim Trimmed Basepairs Bargraph",
            "id": "htstream_cuttrim_bargraph_" + index,
            "ylab": "Basepairs",
            "cpswitch_c_active": False,
            "data_labels": [{"name": "Read 1"}, {"name": "Read 2"}, {"name": "Single End"}],
        }

        # Header
        html = "<h4> CutTrim: Trimmed Basepairs Composition </h4>\n"
        html += "<p>Plots the number of basepairs cut from paired end and single end reads</p>"

        r1_data = {}
        r2_data = {}
        se_data = {}

        # Create Dictionarys for multidataset bargraphs
        for key in json:

            r1_data[key] = {"LT_R1": json[key]["Ct_Left_Trimmed_R1"], "RT_R1": json[key]["Ct_Right_Trimmed_R1"]}

            r2_data[key] = {"LT_R2": json[key]["Ct_Left_Trimmed_R2"], "RT_R2": json[key]["Ct_Right_Trimmed_R2"]}

            se_data[key] = {"LT_SE": json[key]["Ct_Left_Trimmed_SE"], "RT_SE": json[key]["Ct_Right_Trimmed_SE"]}

        # Categories for multidataset bargraphs
        cats = [OrderedDict(), OrderedDict(), OrderedDict()]
        cats[0]["LT_R1"] = {"name": "Left Trimmmed"}
        cats[0]["RT_R1"] = {"name": "Right Trimmmed"}
        cats[1]["LT_R2"] = {"name": "Left Trimmmed"}
        cats[1]["RT_R2"] = {"name": "Right Trimmmed"}
        cats[2]["LT_SE"] = {"name": "Left Trimmmed"}
        cats[2]["RT_SE"] = {"name": "Right Trimmmed"}

        # Create bargraph
        html += bargraph.plot([r1_data, r2_data, se_data], cats, config)

        return html

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
                "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"]
                + json[key]["Paired_end"]["Read2"]["basepairs_out"],
                "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
                "Fraction_Bp_Lost": (total_bp_lost / json[key]["Fragment"]["basepairs_in"]),
            }

            # sample dictionary entry
            stats_json[key] = {
                "Ct_%_BP_Lost" + index: perc_bp_lost,
                "Ct_Notes" + index: json[key]["Program_details"]["options"]["notes"],
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
            "Trimmed Bp Composition Bargraph": self.bargraph(stats_json, index),
            "Overview": overview_dict,
        }

        return section
