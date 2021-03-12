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
    def table(self, json, overall_pe, overall_se, zeroes, index):

        # Table construction. Taken from MultiQC docs.

        # If no polyAT trimmerd, no need for table
        if (overall_pe + overall_se) == 0:
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

        # If PE, add columns
        if overall_pe != 0:
            headers["Pt_%_R1_BP_Lost" + index] = {
                "title": "% R1 of Bp Lost",
                "namespace": "% Bp Lost from R1",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "RdPu",
            }
            headers["Pt_%_R2_BP_Lost" + index] = {
                "title": "% R2 of Bp Lost",
                "namespace": "% Bp Lost from R2",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "Greens",
            }

        # If SE, add collumns
        if overall_se != 0:
            headers["Pt_%_SE_BP_Lost" + index] = {
                "title": "% SE of Bp Lost",
                "namespace": "% Bp Lost from SE",
                "description": "Percentage of total trimmed bps.",
                "suffix": "%",
                "format": "{:,.2f}",
                "scale": "RdPu",
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

            # If no bps lost, prevent zero division
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
                "Pt_%_R1_BP_Lost" + index: total_r1,
                "Pt_%_R2_BP_Lost" + index: total_r2,
                "Pt_%_SE_BP_Lost" + index: total_se,
                "Pt_Avg_BP_Trimmed" + index: total_bp_lost / json[key]["Fragment"]["in"],
                "Pt_Notes" + index: json[key]["Program_details"]["options"]["notes"],
            }

            # Accumulate totals
            overall_pe += total_r1 + total_r2
            overall_se += total_se

        # section and figure function calls
        section = {"Table": self.table(stats_json, overall_pe, overall_se, zeroes, index), "Overview": overview_dict}

        return section
