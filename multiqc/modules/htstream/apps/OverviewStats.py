from collections import OrderedDict
import logging
import numpy as np
import math

from . import htstream_utils
from multiqc import config
from multiqc.plots import table, linegraph, scatter


class OverviewStats:
    def read_and_basepair_reduction(self, json, app_list):

        line_config = {
            "id": "htstream_overview_reduction",
            "smooth_points_sumcounts": False,
            "categories": True,
            "tt_decimals": "{point.y:.0f}'",
            "title": "HTStream: Read and Basepair Reduction",
            "ylab": "Counts",
            "data_labels": [
                {"name": "Reads", "ylab": "Counts", "xlab": "Tool"},
                {"name": "Basepairs", "ylab": "Counts", "xlab": "Tool"},
            ],
        }

        data = [{}, {}]

        # Initialize lists for sample and app order
        samples = list(json["Pipeline Input"].keys())
        app_list = ["Pipeline Input"] + app_list  # preserves order of elements
        read_reducers = json["details"]["read_reducer"]
        bp_reducers = json["details"]["bp_reducer"]

        for samp in samples:

            # initilize read and bp line dicts
            data[0][samp] = {}
            data[1][samp] = {}

            # Iterate through app list, if desired app is found,
            #   grab total read counts and bp counts
            #   and add them to line graphs.

            for app in app_list:

                if app == "Pipeline Input":
                    io = "Input"

                else:
                    io = "Output"

                total_reads = json[app][samp][io + "_Reads"]
                data[0][samp][app] = total_reads

                total_bps = json[app][samp][io + "_Bps"]
                data[1][samp][app] = total_bps

        # Construct html sections
        header = "<h4> {t} </h4>".format(t="Fragment and Basepair Reduction")
        header += """<p> Plots reduction of reads and basepairs across the preprocessing pipeline </p>"""

        # if no apps found in section, create alert div, otherwise, create plots
        if len(app_list) < 2:
            html = title + "\n<br>"
            html = '<div class="alert alert-info">{n}</div>'.format(n=notice)
            return html

        html = header + linegraph.plot(data, line_config)

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
        html = "<hr><h4> Provides scaled statistics collected throughout the preprocessing pipeline, highlighting variable statistics across experiment. </h4>\n"
        html += linegraph.plot(data_dict, line_config) + "\n<br>"

        return html, raw_data

    def execute(self, json, app_list):

        # read_html = self.composition_and_reduction(json, app_list, "read") + "\n<br>"
        # bps_html = self.composition_and_reduction(json, app_list, "bp")
        reduction_html = self.read_and_basepair_reduction(json, app_list)
        line_html, line_data = self.hts_line(json)

        # html = line_html + read_html + bps_html

        html = line_html + reduction_html

        return html, line_data
