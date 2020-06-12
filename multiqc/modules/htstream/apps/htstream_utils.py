import json
import numpy as np

#################################################

""" Utilities for HTStream submodule """

#################################################

# convert json

def resolve(pairs):


	resolved_dict = {}
	index_dict = {}

	# iterates through json key value pairs
	for k, v in pairs:

		if k in index_dict.keys() and "hts_" in k:
			resolved_dict[k + " (" + str(index_dict[k]) + ")"] = v
			index_dict[k] += 1

		elif "hts_" in k:
			resolved_dict[k + " (1)"] = v
			index_dict[k] = 2

		else:
			resolved_dict[k] = v

	return  resolved_dict


#######################################

# prints keys in a pretty way

def key_print(dictionary):

	string = ""

	for key in dictionary.keys():
		string += key + ", "

	string = string[:-2] + "."

	return  string 


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
	html = '<div class="hts_status_header" style="display: inline-block; margin-bottom: 8px;">Sample Checks: </div>'
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
				style = 'border-top-left-radius: 5px; border-bottom-left-radius: 5px; margin-left: -4px;'
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

def qual_by_cycle_html(read, status_div, line_plot, unique_id, button_list, heatmap, index):

	read_header  = " ".join(read.split("_")[1:3])

	# id of switch buttun, named after read type.
	btn_id = "-".join(read.split("_")[:3]).lower()
	

	# section header
	wrapper_html = '<h4> Quality by Cycle: ' + read_header + '</h4>'

	wrapper_html += status_div 

	wrapper_html += '<div class="btn-group hc_switch_group">\n'
	wrapper_html += '<button class="btn btn-default btn-sm active" onclick="htstream_plot_switch(this, {i})" id="htstream_qbc_line_{b}_{u}_btn">Linegraph</button>\n'.format(b=btn_id, i=index, u=unique_id)
	wrapper_html += '<button class="btn btn-default btn-sm " onclick="htstream_plot_switch(this, {i})" id="htstream_qbc_heat_{b}_{u}_btn">Heatmaps</button></div>\n'.format(b=btn_id, i=index, u=unique_id)
	wrapper_html += "<br></br>"

	# this is where the previous html is added to the wrapper html (two separate divs that can be toggled for each graph)
	# line graph div
	wrapper_html += '<div id="htstream_qbc_line_{b}_{u}" class="htstream_fadein">'.format(b=btn_id, u=unique_id)
	wrapper_html += line_plot + "</div>"

	# The heatmaps of this section occur on a per sample basis, meaning we need another subset of buttons to switch between the samples
	heatmap_html = '<div class="btn-group hc_switch_group">\n'

	for buttons in button_list:
		heatmap_html += buttons

	heatmap_html += '</div>\n\n<br></br>\n\n'
	heatmap_html += heatmap

	# heatmap div
	wrapper_html += '<div id="htstream_qbc_heat_{b}_{u}" class="htstream_fadein" style="display:none;">'.format(b=btn_id, u=unique_id)
	wrapper_html += heatmap_html + "</div>"

	final_html = wrapper_html 

	return final_html 


#######################################

# Quality by Base html formatter

def stats_histogram_html(read, data, unique_id, button_list, notice):


	header_dict = {"PE": "Paried End",
				   "R1": "Read 1",
				   "R2": "Read 2",
				   "SE": "Single End"}


	read_header = header_dict[" ".join(read.split("_")[1:2])]


	html = '<h4>Read Length Histogram: '+ read_header + '</h4>'

	if data != {}:

		html += '''<div class="mqc_hcplot_plotgroup">'''
		html += '<div class="btn-group hc_switch_group">\n'

		for buttons in button_list:
			html += buttons

		html += '</div>\n\n'

		data = "["+ json.dumps(data) +"]"

		html += '''<div id="htstream_histogram_{r}_{u}" class="hc-plot"></div></div>'''.format(r = read, u = unique_id)
		html += '''<script type="text/javascript" class="htstream_histogram_content_{r}_{u}">{d}</script>'''.format(r = read, u = unique_id, d = data) 

	html += notice

	return html


#######################################

# pca plot

def pca(matrix, stats_order):


	n, m = matrix.shape # rows, col
	
	to_delete = []

	# normalize 
	for x in range(n):

		row = matrix[x,:]

		# ensure all numbers are between zero and one
		if row[0] > 1.0:
			norm = np.linalg.norm(row)
			row = row / norm

		# remove rows with no variation, also, mean center and normalize variance
		if np.all(row == row[0]):
			to_delete.append(x)
		else:
			matrix[x,:] = (row - np.mean(row) ) / np.std(row)

	# remove indeterminant columns
	to_delete = sorted(to_delete, reverse=True)
	for x in to_delete:
		matrix = np.delete(matrix, x, 0)
		stats_order.remove(stats_order[x])

	n, m = matrix.shape # rows, col

	# sample cov
	cov_mat = np.cov(matrix)

	# for some reason, gives complex numbers, this ensures eigenvalues are real
	eig_val_cov, eig_vec_cov = np.linalg.eig(cov_mat)
	eig_val_cov = eig_val_cov.real
	eig_vec_cov = eig_vec_cov.real

	# creat pair list
	eig_pairs = [(np.abs(eig_val_cov[i]), eig_vec_cov[:,i]) for i in range(len(eig_val_cov))]

	# Sort the (eigenvalue, eigenvector) tuples from high to low
	eig_pairs, stats_order = (list(t) for t in zip(*sorted(zip(eig_pairs, stats_order), key=lambda x: x[0][0], reverse=True)))
	#eig_pairs.sort(key=lambda x: x[0], reverse=True)

	# pc percentages
	eig_sum = sum([eig_pairs[x][0] for x in range(len(eig_pairs))])
	pc_perc = [(eig_pairs[0][0] / eig_sum) * 100, (eig_pairs[1][0] / eig_sum) * 100]

	# eigan vector matrix
	matrix_w = np.hstack((eig_pairs[0][1].reshape(n,1), eig_pairs[1][1].reshape(n,1)))

	# transform
	transformed = matrix_w.T.dot(matrix)

	# VALIDATION
	# VERIFIED using PCA from sklearn.decomposition
	# pca = PCA(n_components=2)
	# transformed = pca.fit_transform(matrix.T)

	return transformed, stats_order, pc_perc







