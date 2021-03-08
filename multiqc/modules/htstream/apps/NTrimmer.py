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
    def table(self, json, overall_pe, overall_se, zeroes, index):

        # Table construction. Taken from MultiQC docs.

        # If no PE and SE removed, return nothin
        if (overall_pe + overall_se) == 0:
            return ""

        headers = OrderedDict()

        # IF values sufficiently small, use raw values
        if zeroes == False:
            headers["Nt_%_BP_Lost" + index] = {
                "title": "% Bp Lost",
                "namespace": "% Bp Lost",
                "description": "Percentage of Input bps (SE and PE) trimmed.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "Greens",
            }
        else:
            headers["Nt_BP_Lost" + index] = {
                "title": "Total Bp Lost",
                "namespace": "Total Bp Lost",
                "description": "Total input bps (SE and PE) trimmed.",
                "format": "{:,.0f}",
                "scale": "Greens",
            }

        # IF PE data, add columns
        if overall_pe != 0:
            headers["Nt_%_R1_BP_Lost" + index] = {
                "title": "% R1 of Bp Lost",
                "namespace": "% Bp Lost from R1",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "RdPu",
            }
            headers["Nt_%_R2_BP_Lost" + index] = {
                "title": "% R2 of Bp Lost",
                "namespace": "% Bp Lost from R2",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "Greens",
            }

        # If SE data, add columns
        if overall_se != 0:
            headers["Nt_%_SE_BP_Lost" + index] = {
                "title": "% SE of Bp Lost",
                "namespace": "% Bp Lost from SE",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "RdPu",
            }

        # IF data is large enough, include avg.
        if zeroes == False:
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

        headers["Nt_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)

    ########################
    # Table Function
    def bargraph(self, json, bps):

        # config dict for bar graph
        config = {
            "title": "HTStream: NTrimmer Trimmed Basepairs Bargraph",
            "id": "htstream_ntrimmer_bargraph",
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

        # If too many samples, get outta here
        if len(json.keys()) > 150:
            html = '<div class="alert alert-warning"> <strong>Notice:</strong> Too many samples for bargraph. </div>'
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
        cats[0]["LT_R1"] = {"name": "Left Trimmmed"}
        cats[0]["RT_R1"] = {"name": "Right Trimmmed"}
        cats[1]["LT_R2"] = {"name": "Left Trimmmed"}
        cats[1]["RT_R2"] = {"name": "Right Trimmmed"}
        cats[2]["LT_SE"] = {"name": "Left Trimmmed"}
        cats[2]["RT_SE"] = {"name": "Right Trimmmed"}

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

            # IF total lost is zero, avoid division by zero
            if total_bp_lost == 0:
                perc_bp_lost = 0
                total_r1 = 0
                total_r2 = 0
                total_se = 0

            else:
                perc_bp_lost = (total_bp_lost / json[key]["Fragment"]["basepairs_in"]) * 100

                total_r1 = (
                    (
                        json[key]["Paired_end"]["Read1"]["basepairs_in"]
                        - json[key]["Paired_end"]["Read1"]["basepairs_out"]
                    )
                    / total_bp_lost
                ) * 100
                total_r2 = (
                    (
                        json[key]["Paired_end"]["Read2"]["basepairs_in"]
                        - json[key]["Paired_end"]["Read2"]["basepairs_out"]
                    )
                    / total_bp_lost
                ) * 100
                total_se = (
                    (json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]) / total_bp_lost
                ) * 100

            # number ofreads discarded
            discarded_reads = json[key]["Single_end"]["discarded"] + json[key]["Paired_end"]["discarded"]

            # are values very small?
            if perc_bp_lost < 0.01 and zeroes == False:
                zeroes = True

            # overview stats
            overview_dict[key] = {
                "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"]
                + json[key]["Paired_end"]["Read2"]["basepairs_out"],
                "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
                "Fraction_Bp_Lost": total_bp_lost / json[key]["Fragment"]["basepairs_in"],
            }

            # sample entry in stats dictionary
            stats_json[key] = {
                "Nt_%_BP_Lost" + index: perc_bp_lost,
                "Nt_BP_Lost" + index: total_bp_lost,
                "Nt_%_R1_BP_Lost" + index: total_r1,
                "Nt_%_R2_BP_Lost" + index: total_r2,
                "Nt_%_SE_BP_Lost" + index: total_se,
                "Nt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
                "Nt_%_Discarded" + index: (discarded_reads / json[key]["Fragment"]["in"]) * 100,
                "Nt_Notes" + index: json[key]["Program_details"]["options"]["notes"],
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
            "Table": self.table(stats_json, overall_pe, overall_se, zeroes, index),
            "Trimmed Reads": self.bargraph(stats_json, (overall_pe + overall_se)),
            "Overview": overview_dict,
        }

        return section
