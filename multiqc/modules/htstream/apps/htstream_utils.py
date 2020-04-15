#################################################

""" Utilities for HTStream submodule """

#################################################

# sample status div creator

def sample_status(samples):

	html = '<br><div class="hts_status_header">Sample Status:</div>'
	html += '<div style="display: inline-block;">\n'

	index = 0
	lim = len(samples.keys()) - 1
	for sample, color in samples.items():

		if index == 0:
			if index == lim:
				style = 'border-top-left-radius: 5px; border-bottom-left-radius: 5px; border-top-right-radius: 5px; border-bottom-right-radius: 5px;'
			else:
				style = 'border-top-left-radius: 5px; border-bottom-left-radius: 5px;'
		elif index == lim:
			style = 'border-top-right-radius: 5px; border-bottom-right-radius: 5px; margin-left: -4px;'
		else:
			style = "margin-left: -4px;"

		html += '<div class="htstream_status" style="background-color: {c}; {r}">{s}</div>\n'.format(s=sample, c=color, r=style)

		index += 1

	html += "</div>\n"
	html += "<br><hr>\n"

	return html

#######################################

# json key collision resolving function

def resolve(pairs):

	resolved_dict = {}

	for k, v in pairs:

		if k == "hts_Stats":

			try:
				resolved_dict[k].append(v)

			except:
				resolved_dict[k] = []
				resolved_dict[k].append(v)

		else:
			resolved_dict[k] = v

	return  resolved_dict
