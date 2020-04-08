function htstream_switch(ele) {

var plot_id = ele.id.split("_b")[0].concat("_heatmap");
var plot_div = ele.parentNode.parentNode.querySelector('.hc-plot');

plot_div.id = plot_id;
plot_div.className = "hc-plot not_rendered hc-heatmap"; 
plot_graph(plot_id);

}
