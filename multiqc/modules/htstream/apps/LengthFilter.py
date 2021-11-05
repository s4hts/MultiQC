from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, bargraph

#################################################

""" LengthFilter submodule for HTStream charts and graphs """

#################################################


class LengthFilter:

    ########################
    # Info about App
    def __init__(self):
        self.info = "Discards reads below a minimum length threshold."
        self.type = "read_reducer"

    ########################
    # Table Function
    def table(self, json, pe_total_loss, se_total_loss, pe_orphaned_total, index):

        # Basic table constructor. See MultiQC docs.
        headers = OrderedDict()

        if pe_total_loss == 0 and se_total_loss == 0 and pe_orphaned_total == 0:
            html = '<div class="alert alert-info"> No reads discarded in any sample. </div>'
            return html

        if pe_total_loss != 0:
            headers["Lf_PE_loss" + index] = {
                "title": "% PE Lost",
                "namespace": "% PE Lost",
                "description": "Percentage of Paired End Reads Lost",
                "format": "{:,.2f}",
                "suffix": "%",
                "scale": "Greens",
            }

        if pe_orphaned_total != 0:
            headers["Lf_PE_Orphaned" + index] = {
                "title": "% PE Orphaned",
                "namespace": "% PE Orphaned",
                "description": "Percentage of Paired End Reads Orphaned (Now Single End)",
                "format": "{:,.2f}",
                "suffix": "%",
                "scale": "RdPu",
            }

        if se_total_loss != 0:
            headers["Lf_SE_loss" + index] = {
                "title": "% SE Lost",
                "namespace": "% SE Lost",
                "description": "Percentage of Single End Reads Lost",
                "format": "{:,.2f}",
                "suffix": "%",
                "scale": "Blues",
            }

        headers["Lf_Notes" + index] = {"title": "Notes", "namespace": "Notes", "description": "Notes"}

        return table.plot(json, headers)

    ########################
    # Bargraphs Function
    def bargraph(self, json, reads_trimmed, index):

        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Composition of Reads Lost Bargraph",
            "id": "htstream_lengthfilter_bargraph_" + index,
            "ylab": "Reads",
        }

        # Title
        html = "<h4> LengthFilter: Composition of Reads Lost </h4>\n"
        html += "<p>Composition of reads removed.</p>"

        # if no overlaps at all are present, return nothing
        if reads_trimmed == 0:
            html += '<div class="alert alert-info"> <strong>Notice:</strong> No reads were removed from samples. </div>'
            return html

        # bargraph dictionary. Exact use of example in MultiQC docs.
        categories = OrderedDict()

        # Colors for sections
        categories["Lf_R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
        categories["Lf_R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
        categories["Lf_SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}

        # Create bargrpah
        html += bargraph.plot(json, categories, config)

        return html

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        # Accumulator vars
        pe_total_loss = 0
        se_total_loss = 0
        pe_orphaned_total = 0

        for key in json.keys():

            # try to calculate reads lost, prevent division by zero
            try:
                reads_lost = (json[key]["Fragment"]["in"] - json[key]["Fragment"]["out"]) / json[key]["Fragment"]["in"]

            except:
                reads_lost = 0

            # try to calculate %  PE reads lost and PE orphaned, prevent division by zero
            try:
                pe_perc_loss = (
                    (json[key]["Paired_end"]["in"] - json[key]["Paired_end"]["out"]) / json[key]["Paired_end"]["in"]
                ) * 100
                pe_orphaned = (
                    (json[key]["Paired_end"]["Read1"]["discarded"] + json[key]["Paired_end"]["Read2"]["discarded"])
                    / json[key]["Paired_end"]["in"]
                ) * 100

            except:
                pe_perc_loss = 0
                pe_orphaned = 0

            # Calculates % SE lost, prevents zero division
            if json[key]["Single_end"]["in"] == 0:
                se_perc_loss = 0

            else:
                se_perc_loss = (json[key]["Single_end"]["discarded"] / json[key]["Single_end"]["in"]) * 100

            # Accumulators accumulating, lol
            pe_total_loss += pe_perc_loss
            se_total_loss += se_perc_loss
            pe_orphaned_total += pe_orphaned

            overview_dict[key] = {
                "PE_Output_Reads": json[key]["Paired_end"]["out"],
                "SE_Output_Reads": json[key]["Single_end"]["out"],
                "Fraction_Reads_Lost": reads_lost,
            }

            # sample entry for stats dictionary
            stats_json[key] = {
                "Lf_PE_loss" + index: pe_perc_loss,
                "Lf_PE_Orphaned" + index: pe_orphaned,
                "Lf_SE_loss" + index: se_perc_loss,
                "Lf_R1_lost": json[key]["Paired_end"]["Read1"]["discarded"],
                "Lf_R2_lost": json[key]["Paired_end"]["Read2"]["discarded"],
                "Lf_SE_lost": json[key]["Single_end"]["discarded"],
                "Lf_Notes" + index: json[key]["Program_details"]["options"]["notes"],
            }

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, pe_total_loss, se_total_loss, pe_orphaned_total, index),
            "Read Composition Bargraph": self.bargraph(stats_json, (pe_total_loss + se_total_loss), index),
            "Overview": overview_dict,
        }

        return section
