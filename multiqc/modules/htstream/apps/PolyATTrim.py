from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" PolyATTrim submodule for HTStream charts and graphs """

#################################################


class PolyATTrim:

    ########################
    # Info about App
    def __init__(self):
        self.info = "Attempts to trim poly-A and poly-T sequences from the end of reads."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def table(self, json, overall, zeroes, index):

        # Table construction. Taken from MultiQC docs.

        # If no polyAT trimmerd, no need for table
        if (overall) == 0:
            return ""

        headers = OrderedDict()

        # If values are small, use raw counts
        if zeroes == False:
            headers["Pt_%_BP_Lost" + index] = {
                "title": "% Bp Lost",
                "namespace": "% Bp Lost",
                "description": "Percentage of Input bps (SE and PE) trimmed.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "Greens",
            }
        else:
            headers["Pt_BP_Lost" + index] = {
                "title": "Total Bp Lost",
                "namespace": "Total Bp Lost",
                "description": "Total input bps (SE and PE) trimmed.",
                "format": "{:,.0f}",
                "scale": "Greens",
            }

        # If values are small, use raw counts
        if zeroes == False:
            headers["Pt_Avg_BP_Trimmed" + index] = {
                "title": "Avg. Bps Trimmed",
                "namespace": "Avg. Bps Trimmed",
                "description": "Average Number of Basepairs Trimmed per Read",
                "format": "{:,.2f}",
                "scale": "Blues",
            }

        headers["Pt_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)


    ########################
    # Bargraph Function
    def bargraph(self, json, bps):

         # config dict for bar graph
        config = {
            "title": "HTStream: PolyATTrim Trimmed Basepairs Bargraph",
            "id": "htstream_polyattrim_bargraph",
            "ylab": "Basepairs",
            "cpswitch_c_active": False,
            "data_labels": [{"name": "Read 1"}, {"name": "Read 2"}, {"name": "Single End"}],
        }

        # Header
        html = "<h4> PolyATTrim: Trimmed Basepairs Composition </h4>\n"
        html += "<p>Plots the number of N bases trimmed from ends of paired end and single end reads.</p>"

        # returns nothing if no reads were trimmed.
        if bps == 0:
            html = '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            return html

        r1_data = {}
        r2_data = {}
        se_data = {}

        # Create dictionaries for multidataset bargraphs
        for key in json:

            r1_data[key] = {"LT_R1": json[key]["Pt_Left_Trimmed_R1"], "RT_R1": json[key]["Pt_Right_Trimmed_R1"]}

            r2_data[key] = {"LT_R2": json[key]["Pt_Left_Trimmed_R2"], "RT_R2": json[key]["Pt_Right_Trimmed_R2"]}

            se_data[key] = {"LT_SE": json[key]["Pt_Left_Trimmed_SE"], "RT_SE": json[key]["Pt_Right_Trimmed_SE"]}

        # Create categores for multidatatset bragraphs
        cats = [OrderedDict(), OrderedDict(), OrderedDict()]
        cats[0]["LT_R1"] = {"name": "R1 Left Trimmmed"}
        cats[0]["RT_R1"] = {"name": "R1 Right Trimmmed"}
        cats[1]["LT_R2"] = {"name": "R2 Left Trimmmed"}
        cats[1]["RT_R2"] = {"name": "R2 Right Trimmmed"}
        cats[2]["LT_SE"] = {"name": "SE Left Trimmmed"}
        cats[2]["RT_SE"] = {"name": "SE Right Trimmmed"}

        # create bargraphs
        html += bargraph.plot([r1_data, r2_data, se_data], cats, config)

        return html

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        # accumulator variable. Used to prevent empty bargraphs
        overall = 0
        zeroes = False

        for key in json.keys():

            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            # If no bps lost, prevent zero division
            if total_bp_lost == 0:
                perc_bp_lost = 0

            else:
                perc_bp_lost = (total_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100
                

            # If values are small, use raw counts
            if perc_bp_lost < 0.01 and zeroes == False:
                zeroes = True

            # Overview stats
            overview_dict[key] = {
                "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"]
                + json[key]["Paired_end"]["Read2"]["basepairs_out"],
                "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
                "Fraction_Bp_Lost": total_bp_lost / json[key]["Fragment"]["basepairs_in"],
            }

            # sample entry in stats dictionary
            stats_json[key] = {
                "Pt_%_BP_Lost" + index: perc_bp_lost,
                "Pt_BP_Lost" + index: total_bp_lost,
                "Pt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
                "Pt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
                "Pt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
                "Pt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
                "Pt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
                "Pt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"],
                "Pt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
                "Pt_Notes" + index: json[key]["Program_details"]["options"]["notes"],
            }

            # Accumulate totals
            overall += total_bp_lost

        # section and figure function calls
        section = {"Table": self.table(stats_json, overall, zeroes, index), 
                   "Trimmed Bp Composition Bargraph": self.bargraph(stats_json, overall),
                   "Overview": overview_dict}

        return section
