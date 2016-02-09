var ar;
var extra;
var graphs_data = {};
var command_list = [];
var graphs = {};
var chi2_hist = [];
var p_hist = [];


$(document).on("pagecontainershow", function(event, ui) {

  if (ui.toPage[0].id === "main_page") {
    update = new Update();
    update.update(['start', 'tree', 'model', 'update_graphs'], {}, event);
  };

  if (ui.toPage[0].id === "landing_page") {
    $(".panel_button").on('click', function(e) {
      $('.panel').toggle();
    });
    $.ajax({
      dataType: 'json',
      url: ($SCRIPT_ROOT + '/_filetree'),
      success: function(data) {
        html = create_filetree(data.tree);
        $('.landing_page_content').html(html);
        $('.landing_page_content a').on('click', function(e) {
          $.mobile.loading('show');
          var url_data = {};
          url_data.model = e.target.name;
          url_data.nFinePoints = $('#nFinePoints').val();
          url_data.round = $('#round').val();
          $.ajax({
            dataType: 'json',
            data: url_data,
            url: ($SCRIPT_ROOT + '/_start'),
            success: function(data) {
              $.mobile.loading('hide');
              $.mobile.pageContainer.pagecontainer("change", "#main_page", {
                transition: "fade"
              });
              location.reload();
            }

          });
        });
      }
    });
  };
});

// updater object/class
function Update() {
  this.complete = 1;
  this.update = function(options, url_data, event, async) {
    options = typeof options !== 'undefined' ? options : false;
    url_data = typeof url_data !== 'undefined' ? url_data : {};
    async = typeof async !== 'undefined' ? async : true;
    event = event || null;

    if (event === 'slidestop') {
      this.complete = 1;
    };

    if (this.complete === 1) {
      this.complete = 0;

      var options_string = '';
      for (var key in options) {
        options_string += options[key] + ';';
      };
      url_data.options = options_string;


      $.ajax({
        context: this,
        dataType: 'json',
        async: async,
        url: ($SCRIPT_ROOT + '/_update'),
        data: url_data,
        success: function(data) {
          ar = $.extend(ar, data.ar);
          extra = $.extend(extra, data.extra);
          this.complete = 1;
          this.success(options, event);
          $.mobile.loading('hide');
        }
      });
    };
  };

  this.success = function(options, event) {

    if (options.indexOf('start') > -1) {
      create_content();
      create_events();
    };

    if (options.indexOf('model') > -1) {

      $('#model_title').html(ar.model.name);
      var description = '';
      for (var _ in ar.model.description) {
        description += ar.model.description[_] + '<br>';
      };
      $('#model_description').html(description);
      $('#model_svg').html(extra.svg);


      $('#filetree').html(create_filetree(extra.tree));
    };
    $('#save, #filetree a').on('click', function(e) {
      var url_data = {};
      var command = ['read'];
      if (e.target.id === 'save') {
        command = ['write'];
        $('#console_editor .CodeMirror').each(function(i, e) {
          url_data.content = e.CodeMirror.getValue();
        });
      };
      url_data.filename = e.target.name;
      update.update(command, url_data, e);
    });

    if (options.indexOf('fit') > -1) {
      if (event.target.name !== 'fit_auto') {
        create_graphs_data(['fit']);
      } else {
        chi2_hist = ar.fit.chi2_hist[0];
        p_hist = ar.fit.p_hist;
      };

    };

    if (options.indexOf('fit_auto') > -1) {
      chi2 = chi2_hist.shift();
      p = p_hist.shift();

      if ((chi2 !== 'NaN') && (chi2 !== undefined)) {
        url_data = {};
        for (var key in ar.pLabel[0]) {
          url_data[ar.pLabel[0][key]] = p[key];
        };
            create_graphs_data(['fit', 'fit_auto', 'observables']);
            this.update(['fit_auto', 'set_sliders', 'update_graphs'], url_data);
      };
    };

    if (options.indexOf('fit_clear') > -1) {
      graphs_data.fit.parameters = [
        [0]
      ];
      graphs_data.fit.chi2plot = [
        [0, ar.chi2fit]
      ];
      for (var i in ar.p[0]) {
        graphs_data.fit.parameters[0].push(ar.p[0][i]);
      };
    };

    if (options.indexOf('read') > -1 || options.indexOf('write') > -1) {
      if ($('#console_editor .CodeMirror').length === 0) {
        CodeMirror.fromTextArea(document.getElementById(
          'editor_area'), {
          theme: 'default',
          lineNumbers: true,
          readOnly: false
        });
      };
      $('#save').attr('name', extra.filename);
      $('#editor_title').html(extra.filename.replace(/^.*(\\|\/|\:)/, ''));
      content = extra.editor_data.content;
      $('#console_editor .CodeMirror').each(function(i, e) {
        e.CodeMirror.refresh();
        e.CodeMirror.setValue(content);
        e.CodeMirror.setSize('auto', 'auto');
      });
    };

    if (options.indexOf('console') > -1) {
      $('#console_out .CodeMirror').each(function(i, e) {
        e.CodeMirror.refresh();
        e.CodeMirror.setValue(extra.console);
        e.CodeMirror.scrollTo(null,
          e.CodeMirror.heightAtLine(e.CodeMirror.getDoc().lineCount()));
      });

      $('#command_multiline .CodeMirror').each(function(i, e) {
        e.CodeMirror.refresh();
      });

    };

    if (options.indexOf('set_sliders') > -1) {
      for (var key in ar.pLabel[0]) {
        $('#' + ar.pLabel[0][key] + '_slider').attr('min', ar.lb[0][key]);
        $('#' + ar.pLabel[0][key] + '_slider').attr('max', ar.ub[0][key]);
        $('#' + ar.pLabel[0][key] + '_slider').val(ar.p[0][key])
      }
    };

    if (options.indexOf('update_graphs') > -1) {
      // check that it wont update the graphs after full fit when auto fitting
      if ( ( event !== null && event.target.name !== 'fit_auto' ) ||
          ( event === null ) ) {
            create_graphs_data(['inputs', 'observables', 'variables']);

            for (var type in graphs) {
              for (var key in graphs[type]) {
                for (var i = 0; i < $('.g_' + type +
                    '_' + key).length; i++) {
                  graphs[type][key][i].updateOptions({
                    'file': graphs_data[type][key]
                  });
                };
              };
            };
      };

    };

    if (options.indexOf('create_graphs') > -1) {
      $('#' + ar.pLabel[0][key] + '_slider').slider('refresh');
      create_graphs_data(['inputs', 'observables', 'variables']);
      create_graphs();
    };

  };
};


function create_content() {

  for (var key in extra.size) {

    var name;
    if (key === 'M') name = 'Model';
    if (key === 'D') name = 'Data';
    if (key === 'C') name = 'Condition';

    var div = $('<div />', {
      'class': 'ui-field-contain',
      text: ''
    });

    var select = $('<select />', {
      name: key,
      class: 'select_mdc'
    });
    var selected;
    for (var i = 1; i <= extra.size[key]; i++) {
      if (i === extra[key]) selected = 'selected';
      var option = $('<option />', {
        value: i - 1,
        text: name + ' ' + i,
        selected: selected
      });
      select.append(option)
    };
    $('#' + key).html(div.append(select));
  };

  var cm_console_out = CodeMirror.fromTextArea(document.getElementById(
    'cm_console_out'), {
    theme: 'default',
    lineNumbers: true,
    readOnly: true
  });
  cm_console_out.setSize('auto', 'auto');

  var command_multiline = CodeMirror.fromTextArea(
    document.getElementById('command_multiline_input'), {
      theme: 'default',
      lineNumbers: true,
      readOnly: false
    });
  command_multiline.setSize('auto', 'auto');

  for (var key in ar.model.u[0]) {
    var div = $('<div />', {
      'class': 'g_container_items g_inputs_' +
        ar.model.u[0][key],
      id: 'g_inputs_' + ar.model.u[0][key],
      text: ar.model.u[0][key]
    });
    $('.g_inputs').append(div);
  };

  for (var key in ar.model.xNames[0]) {
    var div = $('<div />', {
      'class': 'g_container_items g_variables_' +
        ar.model.xNames[0][key],
      id: 'g_variables_' + ar.model.xNames[0][key],
      text: ar.model.xNames[0][key]
    });
    $('.g_variables').append(div);
  };

  for (var key in ar.data.yNames[0]) {
    var div = $('<div />', {
      'class': 'g_container_items g_observables_' +
        ar.data.yNames[0][key],
      id: 'g_observables_' + ar.data.yNames[0][key],
      text: ar.data.yNames[0][key]
    });
    $('.g_observables').append(div);
    div = $('<div />', {
      'class': 'g_container_items g_observables_' +
        ar.data.yNames[0][key],
      id: 'g_observables_' + ar.data.yNames[0][key],
      text: ar.data.yNames[0][key]
    });
    $('.g_fit').append(div);
  };

  var label = $('<label />', {
      for: 'max_iter',
      text: 'Max. fit iterations:'
    });
    var input = $('<input />', {
      'data-type': 'range',
      'data-theme': 'a',
      name: 'max_iter',
      id: 'max_iter',
      min: 1,
      max: 1000,
      value: ar.config.MaxIter,
      step: 1
    });
    $('#max_iter').append(label, input);

  pv = ar.model.pv.map(function(value,index) { return value[0]; });
  pu = ar.model.pu.map(function(value,index) { return value[0]; });
  py = ar.model.py.map(function(value,index) { return value[0]; });
  pystd = ar.model.pystd.map(function(value,index) { return value[0]; });
  pcond = ar.model.pcond.map(function(value,index) { return value[0]; });

  for (var key in ar.pLabel[0]) {
    var label = $('<label />', {
      for: ar.pLabel[0][key] + '_slider',
      text: ar.pLabel[0][key]
    });
    var input = $('<input />', {
      'data-type': 'range',
      'data-theme': 'a',
      name: ar.pLabel[0][key] + '_slider',
      id: ar.pLabel[0][key] + '_slider',
      min: ar.lb[0][key],
      max: ar.ub[0][key],
      value: ar.p[0][key],
      step: 0.1,
      text: ar.pLabel[0][key]
    });
    if (pv.indexOf(ar.pLabel[0][key]) > -1){
      $('#sliders_variables').append(label, input).show();
    } else if (pu.indexOf(ar.pLabel[0][key]) > -1) {
      $('#sliders_inputs').append(label, input).show();
    } else if (py.indexOf(ar.pLabel[0][key]) > -1) {
        $('#sliders_observables').append(label, input).show();
    } else if (pystd.indexOf(ar.pLabel[0][key]) > -1) {
        $('#sliders_observables_std').append(label, input).show();
    } else if (pcond.indexOf(ar.pLabel[0][key]) > -1) {
        $('#sliders_conditions').append(label, input).show();
    };
  };
  $('.panel').trigger('create');
};


function create_events(options) {

  $(".select_mdc").change(function(e) {
    var url_data = {};
    url_data.name = e.target.name
    url_data.value = e.target.value;
    update.update(['change_mdc', 'model', 'create_graphs'], url_data, e)
  });

  for (var key in ar.pLabel[0]) {
    $('#' + ar.pLabel[0][key] + '_slider').on('slidestop change',
      function(e) {
        url_data = {};
        for (var key in ar.pLabel[0]) {
          url_data[ar.pLabel[0][key]] =
            $('#' + ar.pLabel[0][key] + '_slider').val();
        };
        update.update(['update_graphs'], url_data, e);
      });
  };

  $('#max_iter').on('slidestop', function(e) {
    update.update(['max_iter'], {
      max_iter: e.target.value
    }, e);
  });

  $("input[name='toggle_graphs']").change(function(e) {
    $('#graphs_' + e.target.id).toggle();
    create_graphs();
  });
  $("#autoscale").change(function(e) {
      if ( e.target.checked === true){
          for (var type in graphs) {
            for (var key in graphs[type]) {
              for (var i = 0; i < $('.g_' + type +
                  '_' + key).length; i++) {
                graphs[type][key][i].updateOptions({
                  'valueRange': [null, null]
                });
              };
            };
          };
      } else {
          for (var type in graphs) {
            for (var key in graphs[type]) {
              for (var i = 0; i < $('.g_' + type +
                  '_' + key).length; i++) {
                graphs[type][key][i].updateOptions({
                  'valueRange': graphs[type][key][i].yAxisRange()
                });
              };
            };
          };
      };

  });
  $(".panel_button").on('click', function(e) {
    $('.panel').toggle();
  });

  $("button[name='simu_data']").on('click', function(e) {
    update.update(['simu_data', 'create_graphs'], {}, e);
  });

  $("button[name='initial_guess']").on('click', function(e) {
      var url_data = {};
      for (var key in ar.pLabel[0]) {
        var mean = ar.p[0][key];
        var sigma = (ar.ub[0][key] - ar.lb[0][key]) / 100;
        var x =
            Math.random() * (ar.ub[0][key] - ar.lb[0][key]) + ar.lb[0][key];
        x = x * sigma + mean;

        if (x < ar.lb[0][key]) { x = ar.lb[0][key] };
        if (x > ar.ub[0][key]) { x = ar.ub[0][key] };
        url_data[ar.pLabel[0][key]] = x;
      };
      console.log(url_data);
      create_graphs_data(['fit', 'fit_auto']);
    update.update(['set_sliders', 'update_graphs'], url_data, e);
  });

  $("button[name='clear']").on('click', function(e) {
    $.mobile.loading('show');
    update.update(['fit_clear', 'update_graphs'], {}, e);
  });

  $("button[name='fit']").on('click', function(e) {
    $.mobile.loading('show');
    update.update(['fit', 'set_sliders', 'update_graphs'], {}, e);
  });

  $("button[name='fit_auto']").on('click', function(e) {
    $.mobile.loading('show');
    update.update(['fit', 'fit_auto', 'update_graphs'], {}, e);
  });

  $('#tab_model').click(function(e) {
    update.update(['update_graphs', 'model'], {}, e);
  });
  $('#tab_plots').click(function(e) {
    $.mobile.loading('show');
    update.update(['set_sliders', 'create_graphs'], {}, e);
  });
  $('#tab_fit').click(function(e) {
    $.mobile.loading('show');
    update.update(['chi2', 'create_graphs'], {}, e);
  });
  $('#tab_editor').one('click', function(e) {
    var url_data = {};
    url_data.filename = e.target.name;
    update.update(['read'], url_data, e);
  });

  $('#tab_console').on('click', function(e) {
    update.update(['console'], {}, e);
  });

  $('#command_multiline .CodeMirror').each(function(i, e) {
    command_multiline = e.CodeMirror;
  });
  var command = '';
  command_multiline.on('keyHandled', function(name,
    event) {
    if (event === 'Up' &&
      command_multiline.getDoc().lineCount() ==
      1) {
      var end = command_list.pop();
      if (typeof end === 'undefined') {
        end = '';
      } else {
        command_multiline.getDoc().setValue(
          end);
        command_list.unshift(end);
      };
    };
    if (event === 'Down' &&
      command_multiline.getDoc().lineCount() ==
      1) {
      var start = command_list.shift();
      if (typeof start === 'undefined') {
        start = '';
      } else {
        command_multiline.getDoc().setValue(
          command_list[0]);
        command_list.push(start);
      };
    };
    if (event === 'Enter') {
      var command = command_multiline.getDoc()
        .getValue();
      if (command_multiline.getDoc().lineCount() ==
        2) {
        command = command.replace(
          /\n/g, '');
        if (command != '' &&
          command_list[command_list.length -
            1] != command) {
          command_list.push(command);
        };
      };
      command_multiline.getDoc().setValue('');
      $.mobile.loading('show');
      update.update(['console'], {
        command: command
      }, event);
    };
  });

};


function create_graphs_data(options) {

  // create first fit graphs entries if not already existent
  if (!graphs_data.hasOwnProperty('fit')) {
    graphs_data.fit = {};
    graphs_data.fit.parameters = [
      [0]
    ];
    graphs_data.fit.chi2plot = [
      [0, ar.chi2fit]
    ];
    for (var i in ar.p[0]) {
      graphs_data.fit.parameters[0].push(ar.p[0][i]);
    };
  };


  if (options.indexOf('fit') > -1 || options.indexOf('create_graphs') > -1) {
    var fit_step = 1;
    var chi2 = ar.chi2fit;
    if ((chi2_hist[0] !== undefined) && (chi2_hist[0] !== 'NaN')) {
      chi2 = chi2_hist[0];
    };

    if (options.indexOf('fit_auto') === -1) {
      fit_step = parseInt(ar.fit.iter_count);
      chi2 = ar.chi2fit;
    };
    graphs_data.fit.chi2plot.push([
      graphs_data.fit.chi2plot.slice(-1).pop()[0] + fit_step, chi2
    ]);
    arr = [graphs_data.fit.parameters.slice(-1).pop()[0] + fit_step];
    for (var i in ar.pLabel[0]) {
      arr.push(ar.p[0][i]);
    };
    graphs_data.fit.parameters.push(arr);
  };
  if (options.indexOf('observables') > -1) {
    graphs_data['observables'] = {};
    for (var i = 0; i < ar.data.yNames[0].length; i++) {
        graphs_data.observables[ar.data.yNames[0][i]] = ar.plots.observables[i];
    }
  };

  if (options.indexOf('inputs') > -1) {
    graphs_data['inputs'] = {};
    for (var i = 0; i < ar.condition.uFineSimu[0].length; i++) {
      graphs_data.inputs[ar.model.u[0][i]] = [];
      for (var _ in ar.condition.tFine) {
        graphs_data.inputs[ar.model.u[0][i]].push(
          [ar.condition.tFine[_][0], ar.condition.uFineSimu[_][i]]
        );
      };
    };
  };

  if (options.indexOf('variables') > -1) {
    graphs_data['variables'] = {};
    for (var i = 0; i < ar.condition.xFineSimu[0].length; i++) {
      graphs_data.variables[ar.model.xNames[0][i]] = [];
      for (var _ in ar.condition.tFine) {
        graphs_data.variables[ar.model.xNames[0][i]].push(
          [ar.condition.tFine[_][0], ar.condition.xFineSimu[_][i]]
        );
      };
    };

    return graphs_data
  };
};

function create_graphs() {
  var g_settings = {};
  g_settings['global'] = {
    valueRange: [null, null],
    strokeWidth: 2,
    highlightCircleSize: 5,
    highlightSeriesBackgroundAlpha: 0.8,
    highlightSeriesOpts: {
      strokeWidth: 3,
      strokeBorderWidth: 0,
      highlightCircleSize: 5
    }
  };
  labels_fit_parameters = ['n'];
  for (var key in ar.pLabel[0]) {
    labels_fit_parameters.push(ar.pLabel[0][key]);
  };
  for (var plot in graphs_data) {
    graphs[plot] = {};
    for (var key in graphs_data[plot]) {
      graphs[plot][key] = [];
      labels_obs = ['x'];
      for (var i = 0; i < ar.plot.dLink[0].length; i++) {
          labels_obs = labels_obs.concat(["Condition" + i, "Condition" + i]);
      };
      g_settings['observables'] = {
        labels: labels_obs,
        errorBars: true,
        title: key
      };
      g_settings['observables'][key + ' exp'] = {
        fillAlpha: 1,
        strokeWidth: 3,
        drawPoints: true,
        pointSize: 3,
        plotter: singlePointPlotter
      };
      g_settings['fit'] = {
        labels: ['x', key, key + ' exp'],
        errorBars: true,
        title: key
      };
      g_settings['fit'][key + ' exp'] = {
        fillAlpha: 1,
        strokeWidth: 3,
        drawPoints: true,
        pointSize: 3,
        plotter: singlePointPlotter
      };
      g_settings['inputs'] = {
        labels: ['t', key],
        xlabel: ar.model.fu[0],
        title: key
      };
      g_settings['variables'] = {
        labels: ['x', key],
        title: key
      };
      g_settings['chi2plot'] = {
        labels: ['n', '-2*log(L)'],
        errorBars: false,
        title: '-2*log(L)'
      };
      g_settings['parameters'] = {
        labels: labels_fit_parameters,
        errorBars: false,
        title: 'Fit parameters'
      };
      for (var i = 0; i < $('.g_' + plot + '_' + key).length; i++) {
        graphs[plot][key][i] = new Dygraph(
          $('.g_' + plot + '_' + key)[i],
          graphs_data[plot][key],
          $.extend({},
            g_settings['global'], g_settings[plot], g_settings[key]));
      };
    };
  };
};


function create_filetree(tree) {
  var html = '<ul>';
  for (var key in tree['children']) {
    if ('children' in tree['children'][key]) {
      html += '<li>' + tree['children'][key]['name'].replace(/^.*(\\|\/|\:)/, '') + '</li>';
      html += create_filetree(tree['children'][key]);
    } else {
      html += '<li><a name=\'' + tree['children'][key]['name'] +
        '\'  href=\'#\'>' + tree['children'][key]['name'].replace(/^.*(\\|\/|\:)/, '') + '</a></li>';
    };
  };
  html += '</ul>';
  return html
};

// single point error bars for dygraphs
function singleErrorPlotter(e) {
  var ctx = e.drawingContext;
  var points = e.points;
  var g = e.dygraph;
  var color = e.color;
  ctx.save();
  ctx.strokeStyle = e.color;
  for (var i = 0; i < points.length; i++) {
    var p = points[i];
    var center_x = p.canvasx;
    if (isNaN(p.y_bottom)) continue;
    var low_y = g.toDomYCoord(p.yval_minus),
      high_y = g.toDomYCoord(p.yval_plus);
    ctx.beginPath();
    ctx.moveTo(center_x, low_y);
    ctx.lineTo(center_x, high_y);
    ctx.stroke();
  }
  ctx.restore();
};
// single point plotter (without connecting line) for dygraphs
function singlePointPlotter(e) {
  var ctx = e.drawingContext;
  var points = e.points;
  var g = e.dygraph;
  var color = e.color;
  var pointSize = e.dygraph.user_attrs_[e.setName].pointSize;
  var drawPointCallback = Dygraph.Circles.DEFAULT;
  for (var i = 0; i < points.length; i++) {
    var p = points[i];
    ctx.save();
    ctx.beginPath();
    ctx.fillStyle = color;
    ctx.arc(p.canvasx, p.canvasy, pointSize, 0, 2 * Math.PI, false);
    ctx.fill();
  };
  ctx.restore();
};
