from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" NTrimmer submodule for HTStream charts and graphs """

#################################################


class NTrimmer:

    ########################
    # Info about App
    def __init__(self):
        self.info = "Trims reads to the longest subsequence that contains no N's."
        self.type = "bp_reducer"

    ########################
    # Table Function
    def table(self, json, total, zeroes, index):

        # Table construction. Taken from MultiQC docs.

        # If no PE and SE removed, return nothin
        if total == 0:
            return ""

        headers = OrderedDict()

        # IF values sufficiently small, use raw values
        if zeroes == True:

            decimals = "{:,.0f}"

            headers["Nt_BP_Lost" + index] = {
                "title": "Bp Lost",
                "namespace": "Bp Lost",
                "description": "Input bps (SE and PE) trimmed.",
                "format": decimals,
                "scale": "Greens",
            }

        else:

            decimals = "{:,.2f}"

            headers["Nt_%_BP_Lost" + index] = {
                "title": "% Bp Lost",
                "namespace": "% Bp Lost",
                "description": "Percentage of Input bps (SE and PE) trimmed.",
                "suffix": "%",
                "format": decimals,
                "scale": "Greens",
            }


        headers["Nt_Avg_BP_Trimmed" + index] = {
            "title": "Avg. Bps Trimmed",
            "namespace": "Avg. Bps Trimmed",
            "description": "Average Number of Basepairs Trimmed per Read",
            "format": "{:,.2f}",
            "scale": "Blues",
        }

        headers["Nt_%_Discarded" + index] = {
            "title": "% Discarded",
            "namespace": "% Discarded",
            "description": "Percentage of Reads (SE and PE) Discarded",
            "suffix": "%",
            "max": 100,
            "format": "{:,.2f}",
            "scale": "Oranges",
        }

        return table.plot(json, headers)

    ########################
    # Table Function
    def bargraph(self, json, bps, index):

        # config dict for bar graph
        config = {
            "title": "HTStream: NTrimmer Trimmed Basepairs Bargraph",
            "id": "htstream_ntrimmer_bargraph_" + index,
            "ylab": "Basepairs",
            "cpswitch_c_active": False,
            "data_labels": [{"name": "Read 1"}, {"name": "Read 2"}, {"name": "Single End"}],
        }

        # Header
        html = "<h4> NTrimmer: Trimmed Basepairs Composition </h4>\n"
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

            r1_data[key] = {"LT_R1": json[key]["Nt_Left_Trimmed_R1"], "RT_R1": json[key]["Nt_Right_Trimmed_R1"]}

            r2_data[key] = {"LT_R2": json[key]["Nt_Left_Trimmed_R2"], "RT_R2": json[key]["Nt_Right_Trimmed_R2"]}

            se_data[key] = {"LT_SE": json[key]["Nt_Left_Trimmed_SE"], "RT_SE": json[key]["Nt_Right_Trimmed_SE"]}

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
        overall_pe = 0
        overall_se = 0
        zeroes = False

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

            # number ofreads discarded
            discarded_reads = json[key]["Single_end"]["discarded"] + json[key]["Paired_end"]["discarded"]

            # are values very small?
            if perc_bp_lost < 0.01 and zeroes == False:
                zeroes = True

            # overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": total_bp_lost / json[key]["Fragment"]["basepairs_in"],
            }

            # sample entry in stats dictionary
            stats_json[key] = {
                "Nt_%_BP_Lost" + index: perc_bp_lost,
                "Nt_BP_Lost" + index: total_bp_lost,
                "Nt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
                "Nt_%_Discarded" + index: (discarded_reads / json[key]["Fragment"]["in"]) * 100,
                "Nt_R1_lost": total_r1,
                "Nt_R2_lost": total_r2,
                "Nt_SE_lost": total_se,
                "Nt_Left_Trimmed_R1": json[key]["Paired_end"]["Read1"]["leftTrim"],
                "Nt_Right_Trimmed_R1": json[key]["Paired_end"]["Read1"]["rightTrim"],
                "Nt_Left_Trimmed_R2": json[key]["Paired_end"]["Read2"]["leftTrim"],
                "Nt_Right_Trimmed_R2": json[key]["Paired_end"]["Read2"]["rightTrim"],
                "Nt_Left_Trimmed_SE": json[key]["Single_end"]["leftTrim"],
                "Nt_Right_Trimmed_SE": json[key]["Single_end"]["rightTrim"],
            }

            # accumulatr totals
            overall_pe += total_r1 + total_r2
            overall_se += total_se

        # section and figure function calls
        section = {
            "Table": self.table(stats_json, (overall_pe + overall_se), zeroes, index),
            "Trimmed Bps": self.bargraph(stats_json, (overall_pe + overall_se), index),
            "Overview": overview_dict,
        }

        return section
