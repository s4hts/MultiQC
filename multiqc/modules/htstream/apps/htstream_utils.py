#################################################

""" Utilities for HTStream submodule """

#################################################

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
