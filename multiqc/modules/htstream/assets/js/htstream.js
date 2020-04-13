function htstream_div_switch(ele) {

var plot_id = ele.id.split("_b")[0].concat("_heatmap");
var plot_div = ele.parentNode.parentNode.querySelector('.hc-plot');

plot_div.id = plot_id;
plot_div.className = "hc-plot not_rendered hc-heatmap"; 
plot_graph(plot_id);

}

function htstream_plot_switch(ele) {

var plot_id = ele.id.split("_btn")[0];
var read = plot_id.split("_").pop();

if(plot_id.includes('htstream_qbc_line')) {

	var on = document.getElementById(plot_id);
	var off = document.getElementById("htstream_qbc_heat_" + read);
	on.style.display = 'block';
	off.style.display = 'none';

} else {

	console.log("htstream_qbc_line_" + read)

	var on = document.getElementById(plot_id);
	var off = document.getElementById("htstream_qbc_line_" + read);
	on.style.display = 'block';
	off.style.display = 'none';

	var plot_div = on.querySelector('.hc-plot');
	plot_graph(plot_div.id);

}

}
