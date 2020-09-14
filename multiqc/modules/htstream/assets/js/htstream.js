// Functions 

//////////////////////////////////////////////////
// Button Status Switcher

function btn_disable(mode, regex, hide_list, user_hide_list) {

  var exempt_list = ["Base: A", "Base: C", "Base: G", "Base: T", "Base: N", 
                     "PE Reads", "SE Reads", "PE Bps", "SE Bps", "Read 1", "Read 2", "Single End"]; 

  var parent_div = $("#mqc-module-section-htstream"); 
  var unfiltered_btn_divs = parent_div.find("[data-action='set_data'], *[id*=htstream_]");
  var btn_divs = unfiltered_btn_divs.filter(":button").filter(":not(*[onClick*=htstream_plot_switch])");

  $.each(btn_divs.get().reverse(), function(x, value) {

    var btn_text = $(value).text();

    if (mode == "show") {

      if (regex) {

        var show = false;

        for (i = 0; i < user_hide_list.length; i++) { 

          if (btn_text.match(user_hide_list[i])) {
            $(value).prop("disabled", false);
            show = true;
            break;

          }
        }  
        if (show == false  && !exempt_list.includes(btn_text)) {
          $(value).prop("disabled", true);

        }
      } else {

        if (hide_list.indexOf(btn_text) == -1 && exempt_list.indexOf(btn_text) == -1) {

          $(value).prop("disabled", true);

        } else {

          $(value).prop("disabled", false);

        }
      }
    } else {

      if (regex) {

         show = true;

        for (i = 0; i < user_hide_list.length; i++) { 

          if (btn_text.match(user_hide_list[i]) && !exempt_list.includes(btn_text)) { 
            $(value).prop("disabled", true);
            show = false;
            break;

          }
        }

        if (show == true) {
          $(value).prop("disabled", false);

        }  
      } else {

        if (hide_list.indexOf(btn_text) == -1 || exempt_list.indexOf(btn_text) != -1) {

          $(value).prop("disabled", false);

        } else {

          $(value).prop("disabled", true);

        }
      }
    }
  });
}

function btn_activator() {

  var parent_div = $("#mqc-module-section-htstream"); 
  var btn_groups = parent_div.find(".hc_switch_group").filter(":not(*[class*=htstream_exempt])");

  $.each(btn_groups, function(x, value) {

    var btns = $(value).find(".btn");
    var first = true;

    for (i = 0; i < btns.length; i++) { 

      var attr = $(btns[i]).attr('disabled');

      if (typeof attr == typeof undefined && first == true) {
        $(btns[i]).addClass('active');
        first = false;

      } else {
        $(btns[i]).removeClass('active');

      }
    }
  });
}


//////////////////////////////////////////////////
// Div and Plot Switches

function htstream_div_switch(ele, suffix) {

  var plot_id = ele.id.split("_btn")[0];
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

  var plot_div = on.find('.hc-plot');
  plot_graph(plot_div.attr('id'));

}


//////////////////////////////////////////////////
// Get Primers "Samples"

function get_primers() {

  var primers_heatmap_samples = [];
  var temp = [];
  var primers_heatmap_id = Object.keys(mqc_plots).filter(s => s.includes('primers'));

  for (i = 0; i < primers_heatmap_id.length; i++) { 
      temp = mqc_plots[primers_heatmap_id[i]]["xcats"].concat(mqc_plots[primers_heatmap_id[i]]["ycats"]);
      primers_heatmap_samples = primers_heatmap_samples.concat(temp);
    
  }

  return primers_heatmap_samples;
}

//////////////////////////////////////////////////
// Hide, Rename, and Highlight Samples Handlers

var global_f_add = ["Base: A", "Base: C", "Base: G", "Base: T", "Base: N", 
                    "PE Reads", "SE Reads", "PE Bps", "SE Bps"];

var global_on_colors = ["#B62612", "#82A7E0", "#0B8E0B", "#DE7D00", "#000000",
                        "#1EC2D0", "#EA8645", "#1EC2D0",  "#EA8645"];

var sample_num;

// Hide 
$(document).on('mqc_hidesamples', function(e, f_texts, regex_mode){

  mode = $('.mqc_hidesamples_showhide:checked').val() == 'show' ? 'show' : 'hide';
  regex = regex_mode;

  var hide_list = []; 
  var user_hide_list = [];

  if (mode == "hide") {

    var f_add = [];
    hide_list = f_add; 
  
  } else {

    var f_add = global_f_add;
    primers_heatmap_samples = get_primers();

    f_add = f_add.concat(primers_heatmap_samples);

  }

  if (f_texts.length != 0) {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", true);

    user_hide_list = f_texts.slice();

    for (i = 0; i < f_add.length; i++) { 
      f_texts.push(f_add[i]);
    }

    hide_list = f_texts.slice();

  } else {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", false);

  }

  btn_disable(mode, regex, hide_list, user_hide_list);
  btn_activator();
  

});

// Highlight
$(document).on('mqc_highlights', function(e, f_texts, f_cols, regex_mode){

  
  if (f_texts.length != 0 && f_texts.length != (1 + sample_num)) {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", true);

  } else {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", false);

  }
 
  var always_on = global_f_add; 
  var always_on_colors = global_on_colors;


  for (i = 0; i <  always_on.length; i++) { 
      f_texts.push(always_on[i]);
      f_cols.push(always_on_colors[i]);
    }

});


// Rename
$(document).on('mqc_renamesamples', function(e, f_texts, t_texts, regex_mode){

  if (f_texts.length != 0) {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", true);

  } else {

    $('*[id*=htstream_qbc_line]').filter(":button").click();
    $('*[id*=htstream_qbc_heat]').filter(":button").prop("disabled", false);

  }
    
  var primers_on = get_primers();
  var always_on = global_f_add.concat(primers_on);
  var always_on_transform = always_on;

  for (i = 0; i <  always_on.length; i++) { 
      f_texts.push(always_on[i]);
      t_texts.push(always_on_transform[i]);
    }

});


//////////////////////////////////////////////////
// Page Load Magic

$("document").ready(function() {

  // $('*[class*=btn]').each(function(x, ele) {
  //   console.log(ele);
  //   //plot_graph(ele.id, undefined, 100000);
  // });

  var data = JSON.parse($("#htstream_config").text());
  var samples = Object.keys(data["sample_colors"]);
  var colors = data["sample_colors"]
  
  sample_num = data["htstream_number_of_samples"];

  if (samples.length != 0) {  

    for (i = 0; i < samples.length; i++) {
      $('#mqc_col_filters').append('<li style="color:'+colors[samples[i]]+';" id="'+samples[i]+'"><span class="hc_handle"><span></span><span></span></span><input class="f_text" value="'+samples[i]+'" /><button type="button" class="close" aria-label="Close"><span aria-hidden="true">&times;</span></button></li>');
    
    }
   
    $("#mqc_cols_apply").click();
  }

});

