from collections import OrderedDict
import logging
import numpy as np
import math

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, linegraph, scatter


class OverviewStats:
    def composition_and_reduction(self, json, app_list, data_type):

        line_config_1 = {
            "id": "htstream_overview_" + data_type + "_reduction",
            "smooth_points_sumcounts": False,
            "logswitch": False,
            "categories": True,
            "xlab": "Tool",
            "ylab": "Counts",
            "tt_decimals": "{point.y:.0f}'",
        }

        line_config_2 = {
            "id": "htstream_overview_" + data_type + "composition",
            "smooth_points_sumcounts": False,
            "logswitch": True,
            "categories": True,
            "xlab": "Tool",
            "ylab": "Counts",
            "tt_decimals": "{point.y:.0f}'",
            "colors": {"SE Reads": "#1EC2D0", "PE Reads": "#E8961B", "SE Bps": "#1EC2D0", "PE Bps": "#E8961B"},
            "data_labels": [],
        }

        # Define variables so code is less messy later
        if data_type == "read":
            html_title = "Fragment Reduction"
            line_config_1["title"] = "HTStream: " + html_title
            line_config_2["title"] = "HTStream: Fragment Composition"
            reducers = json["details"]["read_reducer"]
            index = "Reads"
            notice = "No Read Reducing Apps were found."

        else:
            html_title = "Basepair Reduction"
            line_config_1["title"] = "HTStream: " + html_title
            line_config_2["title"] = "HTStream: Basepair Composition"
            reducers = json["details"]["bp_reducer"]
            index = "Bps"
            notice = "No Read Reducing Apps were found."

        # Containers for line graph data
        reducing_line = {}
        composition_data_list = []

        # Initialize lists for sample and app order
        samples = list(json["Pipeline Input"].keys())
        app_list = ["Pipeline Input"] + app_list  # preserves order of elements
        app_subset = ["Pipeline Input"]

        for samp in samples:

            # dictionaries for line graphs
            reducing_line[samp] = {}
            composition_line = {"SE " + index: {}, "PE " + index: {}}

            # Iterate through app list, if desired app is found,
            #   grab total read counts and se/pe compisiton
            #   and add them to line graphs.

            for app in app_list:

                include = False

                if app[4:-2] in reducers:
                    total = json[app][samp]["PE_Output_" + index] + json[app][samp]["SE_Output_" + index]
                    app_subset.append(app)
                    prefix = "Output_"
                    include = True

                elif app == "Pipeline Input":
                    total = json[app][samp]["PE_Input_" + index] + json[app][samp]["SE_Input_" + index]
                    prefix = "Input_"
                    include = True

                if include == True:

                    reducing_line[samp][app] = total

                    try:
                        composition_line["SE " + index][app] = json[app][samp]["SE_" + prefix + index]
                    except:
                        composition_line["SE " + index][app] = 0

                    try:
                        composition_line["PE " + index][app] = json[app][samp]["PE_" + prefix + index]
                    except:
                        composition_line["PE " + index][app] = 0

            composition_data_list.append(composition_line)
            line_config_2["data_labels"].append({"name": samp, "ylab": "Counts", "xlab": "Tool"})

        # Construct html sections
        header = "<h4> {t} </h4>".format(t=html_title)
        header += """<p> Provides scaled statistics collected throughout the preprocessing pipeline, highlighting variable statistics across experiment. </p>"""

        btn_label_1 = "Reduction"
        btn_label_2 = "Composition"

        line_1_id = "htstream_comp_table_{b}".format(b=data_type)
        line_2_id = "htstream_comp_line_{b}".format(b=data_type)

        # if no apps found in section, create alert div, otherwise, create plots
        if len(app_subset) < 2:
            html = title + "\n<br>"
            html = '<div class="alert alert-info">{n}</div>'.format(n=notice)
            return html

        else:
            line_1_html = linegraph.plot(reducing_line, line_config_1)
            line_2_html = linegraph.plot(composition_data_list, line_config_2)

        # add html
        html = htstream_utils.multi_plot_html(
            header, btn_label_1, btn_label_2, line_1_id, line_2_id, line_1_html, line_2_html
        )

        return html

    def hts_line(self, json):

        line_config = {
            "id": "htstream_overview_linechart",
            "title": "HTStream: Preprocessing Statistics",
            "smooth_points_sumcounts": False,
            "categories": True,
            "yCeiling": 1,
            "xlab": "Tool",
            "ylab": "Value",
            "tt_decimals": "{point.y:.3f}'",
            "data_labels": [],
        }

        keys = list(json.keys())
        samples_list = list(json[keys[0]].keys())
        row_length = len(samples_list)

        exclude_list = ["PE_Output_Reads", "PE_Output_Bps", "SE_Output_Reads", "SE_Output_Bps"]

        data = [[] for x in range(row_length)]

        stats_order = []
        stats_bool = True

        # iterate through data and find usable apps
        for x in range(row_length):

            sample = samples_list[x]

            for key in keys:

                if key != "Pipeline Input" and key != "details":
                    sample_json = json[key][sample]
                    temp = []

                    for k, v in sample_json.items():
                        if k not in exclude_list:
                            temp.append(v)

                            if stats_bool == True:
                                stats_order.append(str(key + ": " + k))

                    data[x] += temp

            stats_bool = False

        # prepe matrix
        data = np.array(data).T

        # normalize
        data, stats_order, raw_data = htstream_utils.normalize(data, samples_list, stats_order)

        data_dict = {}

        # add data points for each sample
        for x in range(len(samples_list)):

            samp = samples_list[x]
            data_dict[samp] = {}
            values = data[:, x]

            for y in range(len(values)):
                data_dict[samp][stats_order[y]] = values[y]

            line_config["data_labels"].append({"name": samp, "ylab": "Value", "xlab": "Tool"})

        # add html
        html = "<hr><h4> Preprocessing Statistics </h4>\n"
        html += linegraph.plot(data_dict, line_config) + "\n<br>"

        return html, raw_data

    def execute(self, json, app_list):

        read_html = self.composition_and_reduction(json, app_list, "read") + "\n<br>"
        bps_html = self.composition_and_reduction(json, app_list, "bp")
        line_html, line_data = self.hts_line(json)

        html = line_html + read_html + bps_html

        return html, line_data
