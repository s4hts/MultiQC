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
    def table(self, json, pe_total_loss, se_total_loss, index):

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

        if se_total_loss != 0:
            headers["Lf_SE_loss" + index] = {
                "title": "% SE Lost",
                "namespace": "% SE Lost",
                "description": "Percentage of Single End Reads Lost",
                "format": "{:,.2f}",
                "suffix": "%",
                "scale": "Blues",
            }

        if pe_total_loss != 0 and se_total_loss != 0:
            headers["Lf_PE_SE_ratio" + index] = {
                "title": "PE/SE Ratio",
                "namespace": "PE/SE Ratio",
                "description": "Ratio of Paired End Reads lost to Single Edn Reads Lost",
                "format": "{:,.2f}",
                "scale": "Oranges",
            }

        return table.plot(json, headers)

    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        # Accumulator vars
        pe_total_loss = 0
        se_total_loss = 0

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

            except:
                pe_perc_loss = 0

            # Calculates % SE lost, prevents zero division
            if json[key]["Single_end"]["in"] == 0:
                se_perc_loss = 0

            else:
                se_perc_loss = (json[key]["Single_end"]["discarded"] / json[key]["Single_end"]["in"]) * 100

            # Calculate PE / SE ratio, prevents zero division
            try:
                pe_se_ratio = (json[key]["Paired_end"]["discarded"] / json[key]["Single_end"]["discarded"])
            except:
                pe_se_ratio = 0 

            # Accumulators accumulating, lol
            pe_total_loss += pe_perc_loss
            se_total_loss += se_perc_loss

            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_Reads_Lost": reads_lost,
            }

            # sample entry for stats dictionary
            stats_json[key] = {
                "Lf_PE_loss" + index: pe_perc_loss,
                "Lf_SE_loss" + index: se_perc_loss,
                "Lf_PE_SE_ratio" + index: pe_se_ratio,
                "Lf_R1_lost": json[key]["Paired_end"]["Read1"]["discarded"],
                "Lf_R2_lost": json[key]["Paired_end"]["Read2"]["discarded"],
                "Lf_SE_lost": json[key]["Single_end"]["discarded"],
            }

        # sections and figure function calls
        section = {
            "Table": self.table(stats_json, pe_total_loss, se_total_loss, index),
            "Overview": overview_dict,
        }

        return section
