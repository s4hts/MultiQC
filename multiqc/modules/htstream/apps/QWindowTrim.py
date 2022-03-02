from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import bargraph

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
    # Bargraphs Function
    def bargraph(self, json, bps_trimmed, index):

        # configuration dictionary for bar graph
        config = {
            "title": "HTStream: Read Composition of Bps Trimmed Bargraph",
            "id": "htstream_qwindowtrimmer_bargraph_1" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [{"name": "Percentage of Total", "ylab": "Percentage of Total Basepairs"}, 
                            {"name": "Raw Counts", "ylab": "Basepairs"}],
        }

        # Title
        html = ""

        # if no overlaps at all are present, return nothing
        if bps_trimmed == 0:
            html += (
                '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from samples. </div>'
            )
            return html


        perc_data = {}
        read_data = {}

        # Construct data for multidataset bargraph
        for key in json:

            perc_data[key] = {"Perc_R1_lost": json[key]["Qt_Perc_R1_lost"], 
                              "Perc_R2_lost": json[key]["Qt_Perc_R2_lost"], 
                              "Perc_SE_lost": json[key]["Qt_Perc_SE_lost"]}
            read_data[key] = {"R1_lost": json[key]["Qt_R1_lost"],
                              "R2_lost": json[key]["Qt_R2_lost"], 
                              "SE_lost": json[key]["Qt_SE_lost"]}


        # bargraph dictionary. Exact use of example in MultiQC docs.
        categories = [OrderedDict(), OrderedDict()]

        # Colors for sections
        categories[0]["Perc_R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
        categories[0]["Perc_R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
        categories[0]["Perc_SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}
        categories[1]["R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
        categories[1]["R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
        categories[1]["SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}

        # Create bargrpah
        html += bargraph.plot([perc_data, read_data], categories, config)

        return html


    ########################
    # Main Function
    def execute(self, json, index):

        stats_json = OrderedDict()
        overview_dict = {}

        overall_trim = 0

        for key in json.keys():

            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]
            overall_trim += total_bp_lost

            # If no bps lost, prevent zero division
            if total_bp_lost == 0:
                total_r1 = 0
                total_r2 = 0
                total_se = 0
                total_pe = 0

            else:
                total_r1 = (
                    json[key]["Paired_end"]["Read1"]["basepairs_in"] - json[key]["Paired_end"]["Read1"]["basepairs_out"]
                )
                total_r2 = (
                    json[key]["Paired_end"]["Read2"]["basepairs_in"] - json[key]["Paired_end"]["Read2"]["basepairs_out"]
                )

               
                total_se = json[key]["Single_end"]["basepairs_in"] - json[key]["Single_end"]["basepairs_out"]

           
            bp_in = json[key]["Fragment"]["basepairs_in"]

            # overview data
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_R1_Bp_Trimmed_Left": json[key]["Paired_end"]["Read1"]["leftTrim"] / bp_in,
                "Fraction_R1_Bp_Trimmed_Right": json[key]["Paired_end"]["Read1"]["rightTrim"] / bp_in,
                "Fraction_R2_Bp_Trimmed_Left": json[key]["Paired_end"]["Read2"]["leftTrim"] / bp_in,
                "Fraction_R2_Bp_Trimmed_Right": json[key]["Paired_end"]["Read2"]["rightTrim"] / bp_in,
                "Fraction_SE_Bp_Trimmed_Left": json[key]["Single_end"]["leftTrim"] / bp_in,
                "Fraction_SE_Bp_Trimmed_Right": json[key]["Single_end"]["rightTrim"] / bp_in,
            }

            # sample dictionary entry
            stats_json[key] = {
                "Qt_Perc_R1_lost": (total_r1 / bp_in) * 100,
                "Qt_Perc_R2_lost": (total_r2 / bp_in) * 100,
                "Qt_Perc_SE_lost": (total_se / bp_in) * 100,
                "Qt_R1_lost": total_r1,
                "Qt_R2_lost": total_r2,
                "Qt_SE_lost": total_se,
            }


        # sections and figure function calls
        section = {
            "Trimmed Composition": self.bargraph(stats_json, overall_trim, index),
            "Overview": overview_dict,
        }

        return section
