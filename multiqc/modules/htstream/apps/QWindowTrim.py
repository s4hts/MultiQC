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

    # ########################
    # # Bargraphs Function
    # def bargraph(self, json, bps_trimmed, index):

    #     # configuration dictionary for bar graph
    #     config = {
    #         "title": "HTStream: Read Composition of Bps Trimmed Bargraph",
    #         "id": "htstream_qwindowtrimmer_bargraph_" + index,
    #         "ylab": "Percentage of Total Basepairs",
    #         "cpswitch": False,
    #         "data_labels": [
    #             {"name": "Percentage of Total", "ylab": "Percentage of Total Basepairs"},
    #             {"name": "Raw Counts", "ylab": "Basepairs"},
    #         ],
    #     }

    #     # Title
    #     html = ""

    #     # if no overlaps at all are present, return nothing
    #     if bps_trimmed == 0:
    #         html += (
    #             '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from samples. </div>'
    #         )
    #         return html

    #     perc_data = {}
    #     read_data = {}

    #     # Construct data for multidataset bargraph
    #     for key in json:

    #         perc_data[key] = {
    #             "Perc_R1_lost": json[key]["Qt_Perc_R1_lost"],
    #             "Perc_R2_lost": json[key]["Qt_Perc_R2_lost"],
    #             "Perc_SE_lost": json[key]["Qt_Perc_SE_lost"],
    #         }
    #         read_data[key] = {
    #             "R1_lost": json[key]["Qt_R1_lost"],
    #             "R2_lost": json[key]["Qt_R2_lost"],
    #             "SE_lost": json[key]["Qt_SE_lost"],
    #         }

    #     # bargraph dictionary. Exact use of example in MultiQC docs.
    #     categories = [OrderedDict(), OrderedDict()]

    #     # Colors for sections
    #     categories[0]["Perc_R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
    #     categories[0]["Perc_R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
    #     categories[0]["Perc_SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}
    #     categories[1]["R1_lost"] = {"name": "Read 1", "color": "#779BCC"}
    #     categories[1]["R2_lost"] = {"name": "Read 2", "color": "#C3C3C3"}
    #     categories[1]["SE_lost"] = {"name": "Single End", "color": "#D1ADC3"}

    #     # Create bargrpah
    #     html += bargraph.plot([perc_data, read_data], categories, config)

    #     return html

    ########################
    # Table Function
    def bargraph(self, json, bps, index):

        # config dict for bar graph
        config = {
            "title": "HTStream: QWindowTrim Trimmed Basepairs Bargraph",
            "id": "htstream_qwindowtrim_bargraph_" + index,
            "ylab": "Percentage of Total Basepairs",
            "cpswitch": False,
            "data_labels": [
                {"name": "Read 1", "ylab": "Percentage of Total Basepairs"},
                {"name": "Read 2", "ylab": "Percentage of Total Basepairs"},
                {"name": "Single End", "ylab": "Percentage of Total Basepairs"},
            ],
        }

        # Header
        html = ""

        # returns nothing if no reads were trimmed.
        if bps == 0:
            html = '<div class="alert alert-info"> <strong>Notice:</strong> No basepairs were trimmed from any sample. </div>'
            return html

        r1_data = {}
        r2_data = {}
        se_data = {}

        # Create dictionaries for multidataset bargraphs
        for key in json:

            r1_data[key] = {"LT_R1": json[key]["Qt_Left_Trimmed_R1"], "RT_R1": json[key]["Qt_Right_Trimmed_R1"]}

            r2_data[key] = {"LT_R2": json[key]["Qt_Left_Trimmed_R2"], "RT_R2": json[key]["Qt_Right_Trimmed_R2"]}

            se_data[key] = {"LT_SE": json[key]["Qt_Left_Trimmed_SE"], "RT_SE": json[key]["Qt_Right_Trimmed_SE"]}

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

        overall_trim = 0

        for key in json.keys():

            total_bp_lost = json[key]["Fragment"]["basepairs_in"] - json[key]["Fragment"]["basepairs_out"]
            overall_trim += total_bp_lost

            bp_in = json[key]["Fragment"]["basepairs_in"]

            if bp_in != 0:
                fract_r1_bp_left = json[key]["Paired_end"]["Read1"]["leftTrim"] / bp_in
                fract_r1_bp_right = json[key]["Paired_end"]["Read1"]["rightTrim"] / bp_in
                fract_r2_bp_left = json[key]["Paired_end"]["Read2"]["leftTrim"] / bp_in
                fract_r2_bp_right = json[key]["Paired_end"]["Read2"]["rightTrim"] / bp_in
                fract_se_bp_left = json[key]["Single_end"]["leftTrim"] / bp_in
                fract_se_bp_right = json[key]["Single_end"]["rightTrim"] / bp_in

            else:
                fract_r1_bp_left = 0
                fract_r1_bp_right = 0
                fract_r2_bp_left = 0
                fract_r2_bp_right = 0
                fract_se_bp_left = 0
                fract_se_bp_right = 0

                log = logging.getLogger(__name__)
                report = "HTStream: Zero Reads or Basepairs Reported for " + key + "."
                log.error(report)

            # overview data
            overview_dict[key] = {
                "Output_Reads": json[key]["Fragment"]["out"],
                "Output_Bps": json[key]["Fragment"]["basepairs_out"],
                "Fraction_R1_Bp_Trimmed_Left": fract_r1_bp_left,
                "Fraction_R1_Bp_Trimmed_Right": fract_r1_bp_right,
                "Fraction_R2_Bp_Trimmed_Left": fract_r2_bp_left,
                "Fraction_R2_Bp_Trimmed_Right": fract_r2_bp_right,
                "Fraction_SE_Bp_Trimmed_Left": fract_se_bp_left,
                "Fraction_SE_Bp_Trimmed_Right": fract_se_bp_right,
            }

            # sample dictionary entry
            stats_json[key] = {
                "Qt_Left_Trimmed_R1": fract_r1_bp_left * 100,
                "Qt_Right_Trimmed_R1": fract_r1_bp_right * 100,
                "Qt_Left_Trimmed_R2": fract_r2_bp_left * 100,
                "Qt_Right_Trimmed_R2": fract_r2_bp_right * 100,
                "Qt_Left_Trimmed_SE": fract_se_bp_left * 100,
                "Qt_Right_Trimmed_SE": fract_se_bp_right * 100,
            }

        # sections and figure function calls
        section = {"Trimmed Composition": self.bargraph(stats_json, overall_trim, index), "Overview": overview_dict}

        return section
