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

if (plot_id.includes('htstream_qbc_line')) {

	var on = document.getElementById(plot_id);
	var off = document.getElementById("htstream_qbc_heat_" + read);
	on.style.display = 'block';
	off.style.display = 'none';

} else {
	var on = document.getElementById(plot_id);
	var off = document.getElementById("htstream_qbc_line_" + read);
	on.style.display = 'block';
	off.style.display = 'none';
	var plot_div = on.querySelector('.hc-plot');
	plot_graph(plot_div.id);

}
}


function htstream_histogram(read, sample) {

var text = document.getElementsByClassName("htstream_histogram_content_" + read)[0].innerText;
var data = JSON.parse(text)[0];

var container = "htstream_histogram_".concat(read);
var title_read = read.split("_")[1];
var title = "Read Length Histogram (" + title_read + "): " + sample;


Highcharts.chart(container, {
  chart: {
    type: 'column'
  },
  title: {
    text: title
  },
  subtitle: {
    text: ''
  },
  xAxis: {
    categories: data[sample]["bins"],
    crosshair: true,
     title: {
      text: 'Read Length'
    }
  },
  yAxis: {
    min: 0,
    title: {
      text: ''
    }
  },
  tooltip: {
    headerFormat: '<span style="font-size:12px">{point.key} bp</span><table>',
    pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
      '<td style="padding:0"><b>{point.y:.3f}</b></td></tr>',
    footerFormat: '</table>',
    shared: true,
    useHTML: true
  },
  plotOptions: {
    column: {
      pointPadding: 0,
      borderWidth: 0.10,
      groupPadding: 0.10,
      shadow: false
    }
  },
  series: [{
    name: 'log10 Reads',
    data: data[sample]["vals"]

  }]
});

}

$("document").ready(function() {
   $('.active.hist_btn').trigger( "click" );
});



