#################################################

""" Utilities for HTStream submodule """

#################################################

# json key collision resolving function

def resolve(pairs):

	resolved_dict = {}

	# iterates through json key value pairs
	for k, v in pairs:

		# if key is stats, return both entries are added to list
		if k == "hts_Stats":

			try:
				resolved_dict[k].append(v)

			except:
				resolved_dict[k] = []
				resolved_dict[k].append(v)

		# if not stats, business as usual
		else:
			resolved_dict[k] = v

	return  resolved_dict

#######################################

# sample status div creator

def sample_status(samples):

	# color status dictinoary
	color_dict = {
				  "PASS": {"background": "#c3e6c3", "text": "#196F3D"},
				  "QUESTIONABLE": {"background": "#e6dcc3", "text": "#946A04"},
				  "FAIL": {"background": "#e6c3c3", "text": "#C15F5F"}
				  }

	# wrapper divs
	html = '<div class="hts_status_header">Sample Checks:</div>'
	html += '<div style="display: inline-block;">\n'

	# initilize important variables
	index = 0
	lim = len(samples.keys()) - 1

	for sample, status in samples.items():

		# border radius formatting
		if index == 0:
			if index == lim:
				style = 'border-top-left-radius: 5px; border-bottom-left-radius: 5px; border-top-right-radius: 5px; border-bottom-right-radius: 5px;'
			else:
				style = 'border-top-left-radius: 5px; border-bottom-left-radius: 5px;'
		elif index == lim:
			style = 'border-top-right-radius: 5px; border-bottom-right-radius: 5px; margin-left: -4px;'
		else:
			style = "margin-left: -4px;"

		# background and text colors
		color, text = color_dict[status]["background"], color_dict[status]["text"]

		html += '<div class="htstream_status" style="background-color: {c}; color: {t}; {r}">{s}</div>\n'.format(s=sample, c=color, t=text, r=style)

		index += 1

	# close divs
	html += "</div>\n"
	html += "<br>\n"

	# embed in alert div
	notice = '<div class="alert alert-info">{n}</div>'.format(n = html)	

	return notice

#######################################

# Quality by Base html formatter

def qual_by_cycle_html(read, status_div, line_plot, btn_id, button_list, heatmap):

	read_header  = " ".join(read.split("_")[1:3])

	# id of switch buttun, named after read type.
	btn_id = "-".join(read.split("_")[:3]).lower()

	# section header
	wrapper_html = '<h4> Quality by Cycle: ' + read_header + '</h4>'

	wrapper_html += status_div 

	wrapper_html += '<div class="btn-group hc_switch_group">\n'
	wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this)" id="htstream_qbc_line_{r}_btn">Linegraph</button>\n'.format(r=btn_id)
	wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this)" id="htstream_qbc_heat_{r}_btn">Heatmaps</button></div>\n'.format(r=btn_id)
	wrapper_html += "<br></br>"

	# this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
	# line graph div
	wrapper_html += '<div id="htstream_qbc_line_{r}">'.format(r=btn_id)
	wrapper_html += line_plot + "</div>"

	# The heatmaps of this section occur on a per sample basis, meaning we need another subset of buttons to switch between the samples
	heatmap_html = '<div class="btn-group hc_switch_group">\n'

	for buttons in button_list:
		heatmap_html += buttons

	heatmap_html += '</div>\n\n<br></br>\n\n'
	heatmap_html += heatmap

	# heatmap div
	wrapper_html += '<div id="htstream_qbc_heat_{r}" style="display:none;">'.format(r=btn_id)
	wrapper_html += heatmap_html + "</div>"

	final_html = wrapper_html 

	return final_html 


