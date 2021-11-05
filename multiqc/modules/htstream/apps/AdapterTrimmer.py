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
        headers["At_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)

    ########################
    # Bargraphs Function
    def bargraph(self, json, avg_bp_trimmed, index):

        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Trimmed Bp Composition Bargraph",
            "id": "htstream_adaptertrimmer_bargraph_" + index,
            "ylab": "Basepairs",
            "cpswitch_c_active": True,
        }

        # Title
        html = "<h4> AdapterTrimmer: Trimmed Basepairs Composition </h4>\n"
        html += "<p>Composition of basepairs trimmed from the ends of paired end and single end reads.</p>"

        # if no overlaps at all are present, return nothing
        if avg_bp_trimmed == 0:
            html += (
                '<div class="alert alert-info"> <strong>Notice:</strong> No adapters were trimmed from samples. </div>'
            )
            return html

        # bargraph dictionary. Exact use of example in MultiQC docs.
        categories = OrderedDict()

        # Colors for sections
        categories["At_R1"] = {"name": "Read 1", "color": "#779BCC"}
        categories["At_R2"] = {"name": "Read 2", "color": "#C3C3C3"}
        categories["At_SE"] = {"name": "Single End", "color": "#D1ADC3"}

        # Create bargrpah
        html += bargraph.plot(json, categories, config)

        return html

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
                + json[key]["Paired_end"]["Read2"]["adapterTrim"]
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
                "PE_Output_Bps": json[key]["Paired_end"]["Read1"]["basepairs_out"]
                + json[key]["Paired_end"]["Read2"]["basepairs_out"],
                "SE_Output_Bps": json[key]["Single_end"]["basepairs_out"],
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
                "At_Notes" + index: json[key]["Program_details"]["options"]["notes"],
                "At_R1": json[key]["Paired_end"]["Read1"]["adapterBpTrim"],
                "At_R2": json[key]["Paired_end"]["Read2"]["adapterBpTrim"],
                "At_SE": json[key]["Single_end"]["adapterBpTrim"],
            }

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, total, zeroes, index),
            "Bp Composition Bargraph": self.bargraph(stats_json, total, index),
            "Overview": overview_dict,
        }

        return section
