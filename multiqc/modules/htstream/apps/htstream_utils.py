import json, math
import numpy as np
import logging

#################################################

""" Utilities for HTStream submodule """

#################################################

# Logger Initialization
log = logging.getLogger(__name__)

###################################
# convert json
def resolve(pairs):

    resolved_dict = {}
    index_dict = {}

    # iterates through json key value pairs, resolves key conflits
    for k, v in pairs:

        if k in index_dict.keys() and "hts_" in k:
            resolved_dict[k + "_" + str(index_dict[k])] = v
            index_dict[k] += 1

        elif "hts_" in k:
            resolved_dict[k + "_1"] = v
            index_dict[k] = 2

        else:
            resolved_dict[k] = v

    return resolved_dict


###################################
# Json and stats parsing functions
def parse_json(name, f):

    app_dict = {}
    apps = json.loads(f)

    # Will fail if old format is usef
    try:

        # Allows for multiple instances of app, just adds number suffix
        for a in apps:
            i = 1
            app_name = a["Program_details"]["program"] + "_" + str(i)

            if app_name in app_dict.keys():

                while app_name in app_dict.keys():
                    i += 1
                    app_name = a["Program_details"]["program"] + "_" + str(i)

            app_dict[app_name] = a

    except:

        # Used to parse older json files. Will likely be removed in future.
        app_dict = json.loads(f, object_pairs_hook=resolve)
        log.warning("Sample " + name + " uses old json format. Please update to a newer version of HTStream.")

    return app_dict


###################################
# prints keys in a pretty way
def key_print(dictionary):

    string = ""

    for key in dictionary.keys():
        string += key + ", "

    string = string[:-2] + "."

    return string


###################################
# Checks if read lengths are uniform
def uniform(json, read):

    midpoint = 0

    # Check if read lengths are uniform across all samples
    for key in json.keys():

        temp = json[key][read][0]["shape"][-1] * 2

        if midpoint == 0:
            midpoint = temp

        elif midpoint == temp:
            midpoint = midpoint

        else:
            midpoint = -1
            break

    return midpoint


###################################
# Multiplot html formatter
def multi_plot_html(header, btn_1, btn_2, id_1, id_2, graph_1, graph_2, exempt=True):

    # section header
    wrapper_html = header

    # Buttons
    wrapper_html += '<div class="btn-group hc_switch_group {}">\n'.format("htstream_exempt")
    wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this, \'{t}\')" id="{i}_btn">{b}</button>\n'.format(
        i=id_1, t=id_2, b=btn_1
    )
    wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this, \'{t}\')" id="{i}_btn">{b}</button>\n'.format(
        i=id_2, t=id_1, b=btn_2
    )
    wrapper_html += "</div>\n"

    # this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
    # line graph div
    wrapper_html += '<div id="{b}" class="htstream_fadein">'.format(b=id_1)
    wrapper_html += graph_1 + "</div>"

    # this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
    # line graph div
    wrapper_html += '<div id="{b}" class="htstream_fadein" style="display:none;"><br>'.format(b=id_2)
    wrapper_html += graph_2 + "</div>"

    return wrapper_html


###################################
# scale overview linegraph plot
def normalize(data, samples_list, stats_order):

    # Lists for processing data
    n, m = data.shape  # rows, col
    to_delete = []
    raw_data = {}
    factor_dict = {}

    include_list = ["hts_Stats"]

    # format dictionary for output pca stats (raw data)
    for x in range(len(samples_list)):
        raw_data[samples_list[x]] = dict(zip(stats_order, data[:, x]))

    for x in range(n):

        # Get App Name and Statistic
        app = "_".join(stats_order[x].split(": ")[0].split("_")[:-1])
        stat = stats_order[x].split(": ")[-1]
        stat = app + "_" + stat

        row = data[x, :]

        # remove rows with no variation or all zero
        if np.all(row == row[0]) and len(samples_list) > 1:
            to_delete.append(x)
            continue

        elif np.all(row == 0):
            to_delete.append(x)
            continue

        # Get Stat and scale so max value is be between 0.1 and 1
        if stat in factor_dict.keys() and app in include_list:
            scale = factor_dict[stat]

        elif max(row) < 0.1 or max(row) > 1:
            scale = math.floor(-1 * math.log10(max(row)))

            if app in include_list:
                factor_dict[stat] = scale

        else:
            continue

        # Scale row
        row = row * (10 ** (scale))
        factor = " x 10^" + str(scale)
        stats_order[x] = stats_order[x] + factor

        # reaplce row with processed data
        data[x, :] = row

    # remove indeterminant columns
    to_delete = sorted(to_delete, reverse=True)
    for x in to_delete:
        data = np.delete(data, x, 0)
        stats_order.remove(stats_order[x])

    return data, stats_order, raw_data
