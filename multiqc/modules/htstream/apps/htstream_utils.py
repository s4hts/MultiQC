#################################################

""" Utilities for HTStream submodule """

#################################################

def sample_status(samples):

	html = '<br><div class="hts_status_header">Sample Status:</div>'
	html += '<div style="display: inline-block;">\n'

	index = 0
	lim = len(samples.keys()) - 1
	for sample, color in samples.items():

		if index == 0:
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