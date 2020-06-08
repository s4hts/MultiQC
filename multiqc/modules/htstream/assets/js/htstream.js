function htstream_div_switch(ele) {

var plot_id = ele.id.split("_b")[0].concat("_heatmap");
var plot_div = ele.parentNode.parentNode.querySelector('.hc-plot');


plot_div.id = plot_id;
plot_div.className = "hc-plot not_rendered hc-heatmap";
plot_graph(plot_id);

}

function htstream_plot_switch(ele) {

var plot_id = ele.id.split("_btn")[0];
var read = plot_id.split("heat_").pop();

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


console.log(data);
var container = "htstream_histogram_".concat(read);
var title_read = read.split("_")[1];
var title_sample = sample.split("_")[0] + "_" + sample.split("_")[1]
var title = "Read Length Histogram (" + title_read + "): " + title_sample;


// Single end is handled differently, histogram parameters must be set differently
if (title.includes('SE')) {

  var bins = data[sample]["bins"];

  var series_var = [{
        name: 'Single End',
        type: 'histogram',
        xAxis: 1,
        yAxis: 0,
        stack: 0,
        data: data[sample]["vals"],
        type: 'column',
        color: '#84A7CA'
    }]; 

    var padidng_var = {
                      pointPadding: 0,
                      borderWidth: 0.1,
                      groupPadding: 0.1,
                      shadow: false
                      }; 
  
} else {

  var bins = data[sample][0]["bins"];

  var series_var = [{
        name: 'Read 1',
        type: 'histogram',
        xAxis: 1,
        yAxis: 0,
        stack: 0,
        data: data[sample][0]["vals"],
        type: 'column',
        color: '#84A7CA',
    }, {
        name: 'Read 2',
        xAxis: 1,
        yAxis: 0,
        stack: 0,
        data: data[sample][1]["vals"],
        type: 'column',
        color: '#E19CC6',
        showInLegend: true
    }]; 

  var padidng_var = {
                    pointPadding: 0,
                    borderWidth: 0.05,
                    groupPadding: 0.05,
                    shadow: false
                    }; 
}

// high chart histogram function
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
  xAxis: [{
    title: {
      text: 'Data'
    },
    visible: false
  }, {
    title: {
      text: ''
    },
    categories: bins,
    type: 'linear',
    tickPosition: 'outside'
  }],
  yAxis: [{
    min: 0,
    title: {
      text: 'Reads'
    },
  }, {
    visible: false,
    type: 'linear'
  }],

  tooltip: {
    headerFormat: '<span style="font-size:12px">{point.key} bp</span><table>',
    pointFormat: '<tr><td style="color:{series.color};padding:0">{series.name}: </td>' +
      '<td style="padding:0"><b>{point.y:.3f}</b></td></tr>',
    footerFormat: '</table>',
    shared: true,
    useHTML: true
  },
  plotOptions: {
    column: padidng_var
  },
  series: series_var,
});

}

$("document").ready(function() {
   $('.active.hist_btn').trigger( "click" );
});


