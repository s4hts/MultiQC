from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" QWindowTrim submodule for HTStream charts and graphs """

#################################################


class QWindowTrim:

    ########################
    # Info about App
    def __init__(self):
        self.info = "Uses a sliding window approach to remove the low quality ends of reads."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def table(self, json, total, index):

        # Table construction. Taken from MultiQC docs.

        # If no SE and BP lost, dont add table
        if (total) == 0:
            return ""

        headers = OrderedDict()

        headers["Qt_%_BP_Lost" + index] = {
            "title": "% Bp Lost",
            "namespace": "% Bp Lost",
            "description": "Percentage of Input bps (SE and PE) trimmed.",
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "Greens",
        }

        headers["Qt_Avg_BP_Trimmed" + index] = {
            "title": "Avg. Bps Trimmed",
            "namespace": "Avg. Bpss Trimmed",
            "description": "Average Number of Basepairs Trimmed per Read",
            "format": "{:,.2f}",
            "scale": "Blues",
        }
        headers["Qt_%_Discarded" + index] = {
            "title": "% Discarded",
            "namespace": "% Discarded",
            "description": "Percentage of Reads (SE and PE) Discarded",
            "suffix": "%",
            "format": "{:,.2f}",
            "scale": "Oranges",
        }

        headers["Qt_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)

    ########################
    # Bargraphs Function
    def bargraph_1(self, json, reads_trimmed, index):

        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Read Composition of Bps Trimmed Bargraph",
            "id": "htstream_qwindowtrimmer_bargraph_1" + index,
            "ylab": "Reads",
            "cpswitch_c_active": False,
        }

        # Title
        html = "<h4> QWindowTrim: Read Composition of Bps Trimmed </h4>\n"
        html += "<p>Read Composition of basepairs trimmed.</p>"

        # if no overlaps at all are present, return nothing
        if reads_trimmed == 0:
            html += (
                '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from samples. </div>'
            )
            return html

        # bargraph dictionary. Exact use of example in MultiQC docs.
        categories = OrderedDict()

        # Colors for sections
        categories["Qt_R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
        categories["Qt_R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
        categories["Qt_SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}

        # Create bargrpah
        html += bargraph.plot(json, categories, config)

        return html

    ########################
    # Bargraph Function
    def bargraph_2(self, json, bps, index):

        # config dict for bar graph
        config = {
            "title": "HTStream: Trimmed Basepairs Bargraph",
            "id": "htstream_qwindowtrimmer_bargraph_2_" + index,
            "ylab": "Basepairs",
            "cpswitch_c_active": False,
            "data_labels": [{"name": "Read 1"}, {"name": "Read 2"}, {"name": "Single End"}],
        }

        # Header
        html = "<h4> QWindowTrim: Trimmed Basepairs Composition </h4>\n"
        html += "<p>Plots the number of low quality basepairs trimmed from ends of paired end and single end reads.</p>"

        # If too many samples, don't add bargraph
        if len(json.keys()) > 150:
            html += '<div class="alert alert-warning"> <strong>Warning:</strong> Too many samples for bargraph. </div>'
            return html

        r1_data = {}
        r2_data = {}
        se_data = {}

        # Construct data for multidataset bargraph
        for key in json:

            r1_data[key] = {"LT_R1": json[key]["Qt_Left_Trimmed_R1"], "RT_R1": json[key]["Qt_Right_Trimmed_R1"]}

            r2_data[key] = {"LT_R2": json[key]["Qt_Left_Trimmed_R2"], "RT_R2": json[key]["Qt_Right_Trimmed_R2"]}

            se_data[key] = {"LT_SE": json[key]["Qt_Left_Trimmed_SE"], "RT_SE": json[key]["Qt_Right_Trimmed_SE"]}

        # returns nothing if no reads were trimmed.
        if bps == 0:
            html += '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            return html

        # Create categories for multidataset bargraph
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
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        # accumular variable that prevents empty bar graphs
        overall_se = 0
        overall_pe = 0

        for key in json.keys():

            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]

            # If no bps lost, prevent zero division
            if total_bp_lost == 0:
                perc_bp_lost = 0
                total_r1 = 0
                total_r2 = 0
                total_se = 0
                total_pe = 0

            else:
                perc_bp_lost = (total_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100

                total_r1 = (
                    json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]
                )
                total_r2 = (
                    json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]
                )

                left_pe_trimmed = (
                    json[key]["Paired_end"]["Read1"]["leftTrim"] + json[key]["Paired_end"]["Read2"]["leftTrim"]
                )
                right_pe_trimmed = (
                    json[key]["Paired_end"]["Read1"]["rightTrim"] + json[key]["Paired_end"]["Read2"]["rightTrim"]
                )
                total_pe = left_pe_trimmed + right_pe_trimmed

                total_se = json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]
                left_se_trimmed = json[key]["Single_end"]["leftTrim"]
                right_se_trimmed = json[key]["Single_end"]["rightTrim"]

            bp_in = json[key]["Fragment"]["basepairs_in"]

            # overview data
            overview_dict[key] = {
                "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"]
                + json[key]["Paired_end"]["Read2"]["basepairs_out"],
                "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
                "Fraction_R1_Bp_Trimmed_Left": json[key]["Paired_end"]["Read1"]["leftTrim"] / bp_in,
                "Fraction_R1_Bp_Trimmed_Right": json[key]["Paired_end"]["Read1"]["rightTrim"] / bp_in,
                "Fraction_R2_Bp_Trimmed_Left": json[key]["Paired_end"]["Read2"]["leftTrim"] / bp_in,
                "Fraction_R2_Bp_Trimmed_Right": json[key]["Paired_end"]["Read2"]["rightTrim"] / bp_in,
                "Fraction_SE_Bp_Trimmed_Left": json[key]["Single_end"]["leftTrim"] / bp_in,
                "Fraction_SE_Bp_Trimmed_Right": json[key]["Single_end"]["rightTrim"] / bp_in,
            }

            # sample dictionary entry
            stats_json[key] = {
                "Qt_%_BP_Lost" + index: perc_bp_lost,
                "Qt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
                "Qt_Notes" + index: json[key]["Program_details"]["options"]["notes"],
                "Qt_R1_lost": total_r1,
                "Qt_R2_lost": total_r2,
                "Qt_SE_lost": total_se,
                "Qt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
                "Qt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
                "Qt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
                "Qt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
                "Qt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
                "Qt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"],
            }

            # total basepairs accumlation
            overall_pe += total_pe
            overall_se += total_se

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, (overall_pe + overall_se), index),
            "Trimmed Composition": self.bargraph_1(stats_json, (overall_pe + overall_se), index),
            "Trimmed Basepairs": self.bargraph_2(stats_json, (overall_pe + overall_se), index),
            "Overview": overview_dict,
        }

        return section
