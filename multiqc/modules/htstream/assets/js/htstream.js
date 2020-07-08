function htstream_div_switch(ele, suffix) {

  var plot_id = ele.id.split("_b")[0].concat(suffix);
  var parent_node = $(ele).closest('.htstream_fadein');
  var plot_div = parent_node.find('.hc-plot');

  plot_div.attr("id", plot_id);
  plot_div.attr("class", "hc-plot not_rendered hc-heatmap");

  plot_graph(plot_id);

}


function htstream_plot_switch(ele, target) {

  var plot_id = ele.id.split("_btn")[0];
  var on = $("#" + plot_id);
  var off = $("#" + target);

  off.css('display', 'none');
  on.css('display', 'block');


  if (plot_id.includes('htstream_qbc_heat') || plot_id.includes('htstream_comp_line_')) {

    var plot_div = on.find('.hc-plot');
    plot_graph(plot_div.attr('id'));

  }


}


function htstream_histogram(read, sample) {

  var data = JSON.parse($(".htstream_histogram_content_" + read).text())[0];

  var container = "htstream_histogram_".concat(read);
  var title_read = read.split("_")[1];
  var temp_array = sample.split("_");
  var title = "Read Length Histogram (" + title_read + "): " + temp_array.splice(0, temp_array.length - 1).join(by="_");


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
      type: 'column',
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


$(document).on('mqc_hidesamples', function(e, f_texts, regex_mode){

  // console.log(regex_mode);
  
  f_texts.push("Base: A", "Base: C", "Base: G", "Base: T", "Base: N");

  // $("document").ready(function() {

  //   // Base by Cycle Break 
  //   var div = $('*[id*=Base_by_Cycle]').parent().parent();
  //   div.prev('.samples-hidden-warning').remove();
  //   div.css('display', 'block');

  //   var div = $('*[id*=htstream_qbc_heat]');
  //   div.find('.samples-hidden-warning').remove();
  //   div.find('.mqc_hcplot_plotgroup').css('display', 'block');
  //   div.find('.hc-plot-wrapper').css('display', 'block');


  // });

});


$("document").ready(function() {
  $('.active.hist_btn').trigger( "click" );
});



