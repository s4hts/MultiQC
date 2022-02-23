from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" AdapterTrimmer submodule for HTStream charts and graphs """

#################################################


class AdapterTrimmer:

    ########################
    # Info about App
    def __init__(self):
        self.info = (
            "Trims adapters which are sequenced when the fragment insert length is shorter than the read length."
        )
        self.type = "bp_reducer"

    ########################
    # Table Function
    def table(self, json, total, zeroes, index):

        # Table constructor. Just like the MultiQC docs.

        # If total is zero, no need for a plot
        if total == 0:
            return ""

        headers = OrderedDict()

        # Some columns have SUPER small values, add raw counts instead of percentages
        if zeroes == False:
            decimals = "{:,.2f}"

            headers["At_%_BP_Lost" + index] = {
                "title": "% Bp Lost",
                "namespace": "Bp Lost",
                "description": "Percentage of Input bps (SE and PE) trimmed.",
                "scale": "RdPu",
                "suffix": "%",
                "format": decimals,
            }
            headers["At_%_Adapters" + index] = {
                "title": "% Adapters",
                "namespace": "Adapters",
                "description": "Percentage of Reads (SE and PE) with an Adapter",
                "scale": "Blues",
                "suffix": "%",
                "format": decimals,
            }

        else:
            decimals = "{:,.0f}"

            headers["At_BP_Lost" + index] = {
                "title": "Bp Lost",
                "namespace": "Bp Lost",
                "description": "Input bps (SE and PE) trimmed.",
                "scale": "RdPu",
                "format": decimals,
            }
            headers["At_Adapters" + index] = {
                "title": "Adapters",
                "namespace": "Adapters",
                "description": "Reads (SE and PE) with an Adapter",
                "scale": "Blues",
                "format": decimals,
            }

        # More columns
        headers["At_Avg_BP_Trimmed" + index] = {
            "title": "Avg. Bps Trimmed",
            "namespace": "Avg. Bps Trimmed",
            "description": "Average Number of basepairs trimmed from reads",
            "format": "{:,.2f}",
            "scale": "Oranges",
        }

        return table.plot(json, headers)

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}
        total = 0
        zeroes = False

        for key in json.keys():

            frag_in = json[key]["Fragment"]["in"]
            bp_in = json[key]["Fragment"]["basepairs_in"]

            # calculations for reads with adapters and bps trimmed
            adapter_reads = (
                json[key]["Single_end"]["adapterTrim"]
                + json[key]["Paired_end"]["Read1"]["adapterTrim"]
            )  # total reads trimmed
            bp_trimmed = (
                json[key]["Single_end"]["adapterBpTrim"]
                + json[key]["Paired_end"]["Read1"]["adapterBpTrim"]
                + json[key]["Paired_end"]["Read2"]["adapterBpTrim"]
            )  # total basepairs trimmed
            perc_bp_lost = (
                (json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]) / bp_in
            ) * 100

            # if adapter trim is zero, so is the percentage and the avg basepair trimmed.
            #   This prevents division by zero error.
            if adapter_reads == 0:
                perc_adapters = 0
                avg_bp_trimmed = 0

            else:
                perc_adapters = (adapter_reads / frag_in) * 100
                avg_bp_trimmed = bp_trimmed / adapter_reads

            # Accumulate avg bp trimmed
            total += avg_bp_trimmed

            # If percentages are small, use raw counts
            if perc_bp_lost < 0.01 and zeroes == False:
                zeroes = True

            # Overview stats
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Bp_Lost": (bp_in - json[key]["Fragment"]["basepairs_out"]) / bp_in,
                "Fraction_PE_Bp_Trimmed": (
                    json[key]["Paired_end"]["Read1"]["adapterBpTrim"]
                    + json[key]["Paired_end"]["Read2"]["adapterBpTrim"]
                )
                / bp_in,
                "Fraction_PE_Read_Trimmed": (
                    json[key]["Paired_end"]["Read1"]["adapterTrim"] + json[key]["Paired_end"]["Read2"]["adapterTrim"]
                )
                / frag_in,
                "Fraction_SE_Bp_Trimmed": json[key]["Single_end"]["adapterBpTrim"] / bp_in,
                "Fraction_SE_Read_Trimmed": json[key]["Single_end"]["adapterTrim"] / frag_in,
            }

            # sample dictionary entry
            stats_json[key] = {
                "At_%_BP_Lost" + index: perc_bp_lost,
                "At_%_Adapters" + index: perc_adapters,
                "At_BP_Lost" + index: bp_trimmed,
                "At_Adapters" + index: adapter_reads,
                "At_Avg_BP_Trimmed" + index: avg_bp_trimmed,
                "At_R1": json[key]["Paired_end"]["Read1"]["adapterBpTrim"],
                "At_R2": json[key]["Paired_end"]["Read2"]["adapterBpTrim"],
                "At_SE": json[key]["Single_end"]["adapterBpTrim"],
            }

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, total, zeroes, index),
            "Overview": overview_dict,
        }

        return section
