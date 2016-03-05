var ar;
var extra;
var graphs_data = {};
var fit_data = {};
var command_list = [];
var graphs = {};
var chi2_hist = [];
var p_hist = [];
var cpal = [
    "#5DA5DA",
    "#FAA43A",
    "#60BD68",
    "#F17CB0",
    "#B2912F",
    "#B276B2",
    "#DECF3F",
    "#F15854",
    "#4D4D4D",
    "#33315B",
    "#76515B",
    "#306DCE",
];
$(document).on("pagecontainershow", function(event, ui) {
    if (ui.toPage[0].id === "main_page") {
        update = new Update();
        update.update([
            'start',
            'tree',
            'model',
            'update_graphs',
        ], {}, event);
    };
    if (ui.toPage[0].id === "landing_page") {
        $(".panel_button").on('click', function(e) {
            $('.panel').toggle();
        });
        $.ajax({
            dataType: 'json',
            url: ($SCRIPT_ROOT + '/_filetree'),
            success: function(data) {
                html = select_model(data.tree);
                $('.choose_model').html(html);
                $('.choose_model').enhanceWithin();
                $("#start").on('click', function(e) {
                    if ($( "#select_model" ).val() === 'Select model...'){
                        alert("Please select a model first.");
                        return
                    };
                    $.mobile.loading('show');
                    var url_data = {};
                    url_data.model = $( "#select_model" ).val();
                    url_data.nFinePoints = $('#nFinePoints').val();
                    url_data.round = $('#round').val();
                    url_data.compile = $('#compile').val();
                    $.ajax({
                        dataType: 'json',
                        data: url_data,
                        url: ($SCRIPT_ROOT + '/_start'),
                        success: function(data) {
                            $.mobile.loading('hide');
                            if (data.status.arSimu !== null) {
                                alert("Simulation failed. Try to use higher " +
                                "values for nFinePoints.");
                                $.mobile.pageContainer.pagecontainer("change", "#landing_page", {transition: "fade"});
                                return
                            };
                            if (data.status.nFinePoints_min) {
                                alert("nFinePoints set too low for data " +
                                "simulation. nFinePoints has been set to " +
                                data.status.nFinePoints_min + " instead.");
                            };
                            $.mobile.pageContainer.pagecontainer("change", "#main_page", {transition: "fade"});
                            location.reload();
                        },
                    });
                });
            },
        });
    };
});
// updater object/class
function Update() {
    this.complete = 1;
    this.update = function(options, url_data, event, async) {
        options = typeof options !== 'undefined'
            ? options
            : false;
        url_data = typeof url_data !== 'undefined'
            ? url_data
            : {};
        async = typeof async !== 'undefined'
            ? async
            : true;
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
                    if (data.status.arSimuData === false) {
                        alert("arSimuData failed - no data has been set up.");
                    };
                    if (data.status.arSimu === 'p_reset') {
                        alert("Simulation failed. Probably caused by wrong " +
                              "paraemter settings. Parameters will be " +
                              "reseted from ar.d2d_presenter.p and the "+
                              "content reloaded.");
                        location.reload();
                        return
                    } else if (data.status.arSimu === 'arLoad') {
                        alert("Simulation failed. Probably caused by wrong " +
                        "paraemter settings. The last saved state from " +
                        "ar.config.savepath will be used and the content " +
                        "reloaded.");
                        location.reload();
                    }
                    else if(data.status.arSimu === 'broken'){
                        alert("Simulation failed. Probably caused by wrong " +
                        "paraemter settings. No working configuration has " +
                        "been found, please rerun the setup.");
                    };
                    this.complete = 1;
                    this.success(options, event);
                    $.mobile.loading('hide');
                },
            });
        };
    };
    this.success = function(options, event) {

        if (options.indexOf('start') > -1) {
            create_content();
            create_events();
        };
        if (options.indexOf('model') > -1) {
            $('#model_title').html(ar['model.name']);
            var description = '';
            for (var _ in ar['model.description']) {
                description += ar['model.description'][_] + '<br>';
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
            chi2_hist = ar['fit.chi2_hist'];
            p_hist = ar['fit.p_hist'];
          };

        };

        if (options.indexOf('fit_auto') > -1) {
          chi2 = chi2_hist.shift()[0];
          p = [];
          for (var i in p_hist) {
              p.push(p_hist[i].shift());
          };
          if (chi2 !== null) {
            url_data = {};
            for (var key in ar['pLabel']) {
                url_data[ar['pLabel'][key]] = p[key];
            };
                create_graphs_data(['fit_auto', 'observables']);
                this.update(['chi2', 'fit_auto', 'set_sliders', 'update_graphs'], url_data);
          };
        };
        if (options.indexOf('fit_clear') > -1) {
          graphs_data.fit.parameters = [
            [0]
          ];
          for (var i in ar['p']) {
            graphs_data.fit.parameters[0].push(ar['p'][i]);
          };
          graphs_data.fit.chi2plot = [
            [0, ar['fit.chi2_hist']]
          ];

        };

        if (options.indexOf('read') > -1 || options.indexOf('write') > -1) {
            if ($('#console_editor .CodeMirror').length === 0) {
                CodeMirror.fromTextArea(document.getElementById('editor_area'), {
                    theme: 'default',
                    lineNumbers: true,
                    readOnly: false,
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
                e.CodeMirror.scrollTo(null, e.CodeMirror.heightAtLine(e.CodeMirror.getDoc().lineCount()));
            });
            $('#command_multiline .CodeMirror').each(function(i, e) {
                e.CodeMirror.refresh();
            });
        };
        if (options.indexOf('set_sliders') > -1) {
            for (var key in ar['pLabel']) {
                $('#' + ar['pLabel'][key] + '_slider').attr('min', ar['lb'][key][0]);
                $('#' + ar['pLabel'][key] + '_slider').attr('max', ar['ub'][key][0]);
                $('#' + ar['pLabel'][key] + '_slider').val(ar['p'][key][0]);
                $('#' + ar['pLabel'][key] + '_slider').attr('alt', 'refresh'); // hack to inform the event that this was a refresh
                $('#' + ar['pLabel'][key] + '_slider').slider("refresh");
            }
        };
        if (options.indexOf('create_graphs') > -1) {
            $('#' + ar['pLabel'][key] + '_slider').slider('refresh');
            create_graphs_data(['inputs', 'observables', 'variables']);
            create_graphs();
        };
        if (options.indexOf('update_graphs') > -1) {
          // check that it wont update the graphs after full fit when auto fitting
          if ( ( event !== null && event.target.name !== 'fit_auto' ) ||
              ( event === null ) ) {
                  if ($("ul li.ui-state-active")[0].children[0].id === 'tab_fit') {
                      create_graphs_data(['observables', 'fit']);
                      var type = 'fit';
                        for (var key in graphs[type]) {
                          for (var i = 0; i < $('.g_' + type +
                              '_' + key).length; i++) {
                            graphs[type][key][i].updateOptions({
                              'file': graphs_data[type][key]
                            });
                          };
                        };
                        var type = 'observables';
                        for (var key in graphs[type]) {
                          for (var i = 0; i < $('.g_' + type +
                              '_' + key).length; i++) {
                            graphs[type][key][i].updateOptions({
                              'file': graphs_data[type][key]
                            });
                          };
                        };
                    } else if ($("ul li.ui-state-active")[0].children[0].id === 'tab_plots') {
                     if ($('#inputs')[0].checked) {

                         create_graphs_data(['inputs']);
                         var type = 'inputs';
                         for (var key in graphs[type]) {
                           for (var i = 0; i < $('.g_' + type +
                               '_' + key).length; i++) {
                             graphs[type][key][i].updateOptions({
                               'file': graphs_data[type][key]
                             });
                           };
                         };
                    };
                    if ($('#variables')[0].checked) {
                        create_graphs_data(['variables']);
                        var type = 'variables';
                        for (var key in graphs[type]) {
                          for (var i = 0; i < $('.g_' + type +
                              '_' + key).length; i++) {
                            graphs[type][key][i].updateOptions({
                              'file': graphs_data[type][key]
                            });
                          };
                        };
                   };
                   if ($('#observables')[0].checked) {
                       create_graphs_data(['observables']);
                       var type = 'observables';
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
        };
    };
};
function create_content() {

        var div = $('<label />', {
            'for': 'MODEL',
            text: 'Model:',
        });
        var select = $('<select />', {
            name: 'MODEL',
            class: 'select_mdc',
        });
        for (var i = 1; i <= extra.size['MODEL']; i++) {
            if (i === extra['MODEL']) {
                var option = $('<option />', {
                    value: i,
                    text: extra.size['MODELNAMES'][i-1],
                    selected: '',
                });
            } else {
                var option = $('<option />', {
                    value: i,
                    text: extra.size['MODELNAMES'][i-1],
                });
            };
            select.append(option)
        };
        $('#MODEL').html(div.append(select));
        var div = $('<label />', {
            'for': 'DSET',
            text: 'Data:',
        });
        var select = $('<select />', {
            name: 'DSET',
            class: 'select_mdc',
        });
        for (var i = 1; i <= extra.size['DSET']; i++) {
            if (i === extra['DSET']) {
                var option = $('<option />', {
                    value: i,
                    text: extra.size['model.plot.name'][i-1],
                    selected: '',
                });
            } else {
                var option = $('<option />', {
                    value: i,
                    text: extra.size['model.plot.name'][i-1],
                });
            };
            select.append(option)
        };
        $('#DSET').html(div.append(select));

    var cm_console_out = CodeMirror.fromTextArea(document.getElementById('cm_console_out'), {
        theme: 'default',
        lineNumbers: true,
        readOnly: true,
    });
    cm_console_out.setSize('auto', 'auto');
    var command_multiline = CodeMirror.fromTextArea(document.getElementById('command_multiline_input'), {
        theme: 'default',
        lineNumbers: true,
        readOnly: false,
    });
    command_multiline.setSize('auto', 'auto');
    for (var key in ar['model.u']) {
        var div = $('<div />', {
            'class': 'g_container_items g_inputs_' + ar['model.u'][key],
            id: 'g_inputs_' + ar['model.u'][key],
        });
        $('.g_inputs').append(div);
    };
    for (var key in ar['model.xNames']) {
        var div = $('<div />', {
            'class': 'g_container_items g_variables_' + ar['model.xNames'][key],
            id: 'g_variables_' + ar['model.xNames'][key],
        });
        $('.g_variables').append(div);
    };
    for (var key in ar['model.data.yNames'][0]) {
        var div = $('<div />', {
            'class': 'g_container_items g_observables_' + ar['model.data.yNames'][0][key],
            id: 'g_observables_' + ar['model.data.yNames'][0][key],
        });
        $('.g_observables').append(div);
        div = $('<div />', {
            'class': 'g_container_items g_observables_' + ar['model.data.yNames'][0][key],
            id: 'g_observables_' + ar['model.data.yNames'][0][key],
        });
        $('.g_fit').append(div);
    };
    var label = $('<label />', {
        for: 'max_iter',
        text: 'Max. fit iterations:',
    });
    var input = $('<input />', {
        'data-type': 'range',
        'data-theme': 'a',
        name: 'max_iter',
        id: 'max_iter',
        min: 1,
        max: 1000,
        value: ar['config.optim.MaxIter'],
        step: 1,
    });
    $('#max_iter').append(label, input);
    pv = ar['model.pv'].map(function(value, index) {
        return value;
    });
    pu = ar['model.pu'].map(function(value, index) {
        return value;
    });
    py = ar['model.py'].map(function(value, index) {
        return value;
    });
    pystd = ar['model.pystd'].map(function(value, index) {
        return value;
    });
    pcond = ar['model.pcond'].map(function(value, index) {
        return value;
    });
    pc = ar['model.pc'].map(function(value, index) {
        return value;
    });
    for (var key in ar['pLabel']) {
        var label = $('<label />', {for: ar['pLabel'][key] + '_slider',
            text: ar['pLabel'][key],
        });
        var input = $('<input />', {
            'data-type': 'range',
            'data-theme': 'a',
            name: ar['pLabel'][key] + '_slider',
            id: ar['pLabel'][key] + '_slider',
            min: ar['lb'][key][0],
            max: ar['ub'][key][0],
            value: ar['p'][key][0],
            step: 0.00001,
            text: ar['pLabel'][key],
        });
        if (pc.indexOf(ar['pLabel'][key]) > -1) {
            $('#sliders_compartments').append(label, input).show();
        } else if (pv.indexOf(ar['pLabel'][key]) > -1) {
            $('#sliders_variables').append(label, input).show();
        } else if (pu.indexOf(ar['pLabel'][key]) > -1) {
            $('#sliders_inputs').append(label, input).show();
        } else if (py.indexOf(ar['pLabel'][key]) > -1) {
            $('#sliders_observables').append(label, input).show();
        } else if (pystd.indexOf(ar['pLabel'][key]) > -1) {
            $('#sliders_observables_std').append(label, input).show();
        } else if (pcond.indexOf(ar['pLabel'][key]) > -1) {
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
        update.update([
            'change_mdc', 'model',
        ], url_data, e)
    });
    for (var key in ar['pLabel']) {
        $('#' + ar['pLabel'][key] + '_slider').on('slidestop change', function(e) {
            if (e.target.alt === "refresh") {
                for (var key in ar['pLabel']) {
                    $('#' + ar['pLabel'][key] + '_slider').attr('alt', '');
                };
            } else {
                url_data = {};
                for (var key in ar['pLabel']) {
                    url_data[ar['pLabel'][key]] = $('#' + ar['pLabel'][key] + '_slider').val();
                };
                update.update(['update_graphs'], url_data, e);
            };
        });
    };
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
    $('#max_iter').on('slidestop', function(e) {
        update.update(['max_iter'], {
            max_iter: e.target.value
        }, e);
    });
    $("input[name='toggle_graphs']").change(function(e) {
        $('#graphs_' + e.target.id).toggle();
        if ($('#graphs_' + e.target.id).is(":visible")) {
        update.update(['create_graphs'], {}, e);
    };
    });
    $(".panel_button").on('click', function(e) {
        $('.panel').toggle();
    });
    $("button[name='simu_data']").on('click', function(e) {
        update.update(['simu_data', 'update_graphs'], {}, e);
    });
    $("button[name='initial_guess']").on('click', function(e) {
        var url_data = {};
        for (var key in ar['p']) {
            var mean = ar['p'][key][0];
            var sigma = (ar['ub'][key][0] - ar['lb'][key][0]) / 100;
            var x = Math.random() * (ar['ub'][key][0] - ar['lb'][key][0]) + ar['lb'][key][0];
            x = x * sigma + mean;
            if (x < ar['lb'][key][0]) {
                x = ar['lb'][key][0]
            };
            if (x > ar['ub'][key][0]) {
                x = ar['ub'][key][0]
            };
            url_data[ar['pLabel'][key]] = x;
        };
        update.update([
            'chi2', 'set_sliders', 'update_graphs',
        ], url_data, e);
    });
    $("button[name='clear']").on('click', function(e) {
        $.mobile.loading('show');
        update.update([
            'fit_clear', 'update_graphs',
        ], {}, e);
    });
    $("button[name='fit']").on('click', function(e) {
        $.mobile.loading('show');
        update.update([
            'fit', 'set_sliders', 'update_graphs',
        ], {}, e);
    });
    $("button[name='setup']").on('click', function(e) {
        $.mobile.loading('show');
        update.update([
            'setup', 'set_sliders', 'update_graphs',
        ], {}, e);
    });
    $("button[name='fit_auto']").on('click', function(e) {
        $.mobile.loading('show');
        update.update([
            'fit', 'fit_auto', 'update_graphs',
        ], {}, e);
    });
    $('#tab_model').click(function(e) {
        update.update(['model',
        ], {}, e);
    });
    $('#tab_plots').click(function(e) {
        $.mobile.loading('show');
        update.update([
            'set_sliders', 'create_graphs',
        ], {}, e);
    });
    $('#tab_fit').click(function(e) {
        $.mobile.loading('show');
        update.update([
            'chi2', 'create_graphs', 'update_graphs',
        ], {}, e);
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
    command_multiline.on('keyHandled', function(name, event) {
        if (event === 'Up' && command_multiline.getDoc().lineCount() == 1) {
            var end = command_list.pop();
            if (typeof end === 'undefined') {
                end = '';
            } else {
                command_multiline.getDoc().setValue(end);
                command_list.unshift(end);
            };
        };
        if (event === 'Down' && command_multiline.getDoc().lineCount() == 1) {
            var start = command_list.shift();
            if (typeof start === 'undefined') {
                start = '';
            } else {
                command_multiline.getDoc().setValue(command_list[0]);
                command_list.push(start);
            };
        };
        if (event === 'Enter') {
            var command = command_multiline.getDoc().getValue();
            if (command_multiline.getDoc().lineCount() == 2) {
                command = command.replace(/\n/g, '');
                if (command != '' && command_list[command_list.length - 1] != command) {
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
      [0, ar['chi2fit']]
    ];
    for (var i in ar['p']) {
      graphs_data.fit.parameters[0].push(ar['p'][i][0]);
    };
  }
  else{

  if (options.indexOf('fit') > -1 || options.indexOf('create_graphs') > -1) {
    var fit_step = 1;
    var chi2 = ar['chi2fit'];
    if ((chi2_hist !== undefined) && (chi2_hist !== 'NaN')) {
      chi2 = chi2_hist;
    };

    if (options.indexOf('fit_auto') === -1) {
        fit_step = parseInt(ar['fit.iter_count']) || 1;
    };
    graphs_data.fit.chi2plot.push([
      graphs_data.fit.chi2plot.slice(-1).pop()[0] + fit_step, ar['chi2fit']
    ]);

    var arr = [graphs_data.fit.parameters.slice(-1).pop()[0] + fit_step];
    for (var i in ar['pLabel']) {
      arr.push(ar['p'][i][0]);
    };
    graphs_data.fit.parameters.push(arr);
  };
  };
  if (options.indexOf('observables') > -1) {
    graphs_data['observables'] = {};
    for (var i = 0; i < ar['model.data.yNames'][0].length; i++) {
        graphs_data.observables[ar['model.data.yNames'][0][i]] = ar.plots.observables[i];
    }
  };

  if ((options.indexOf('inputs') > -1) && (ar['model.u'] != undefined)) {
    graphs_data['inputs'] = {};
    for (var i = 0; i < ar['model.u'].length; i++) {
        graphs_data.inputs[ar['model.u'][i]] = ar.plots.inputs[i];
    };
  };

  if (options.indexOf('variables') > -1) {
    graphs_data['variables'] = {};
    for (var i = 0; i < ar['model.xNames'].length; i++) {
        graphs_data.variables[ar['model.xNames'][i]] = ar.plots.variables[i];
    };
  };
      return graphs_data
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

  for (var key in ar['pLabel']) {
    labels_fit_parameters.push(ar['pLabel'][key]);
  };

  cp_parameters = cpal.slice();
  for (var plot in graphs_data) {
      graphs[plot] = {};

    for (var key in graphs_data[plot]) {
        graphs[plot][key] = [];
        cp = cpal.slice();

        if (plot === 'fit'){
            g_settings['fit'] = {
              labels: ['x', key, key + ' exp'],
              errorBars: true,
              title: key
            };
            g_settings['fit']['series'] = {};
            g_settings['fit']['series'][key + ' exp'] = {
              fillAlpha: 1,
              strokeWidth: 3,
              drawPoints: true,
              pointSize: 3,
              plotter: singlePointPlotter
            };
            g_settings['chi2plot'] = {
              labels: ['n', 'Chi2'],
              errorBars: false,
              title: 'Chi2'
            };
            g_settings['parameters'] = {
              labels: labels_fit_parameters,
              errorBars: false,
              title: 'Fit parameters'
            };
        };
        if (plot === 'variables'){
            labels_var = ['x'];
            for (var i = 0; i < (ar.plots.variables[0][0].length-1); i++) {
                labels_var = labels_var.concat(["Condition " + i]);
            };
            g_settings['variables'] = {
              labels: labels_var,
              title: key
            };
        };

        if (plot === 'observables'){

            series_settings_obs = {};
            labels_obs = ['x'];
            for (var i = 0; i < (ar.plots.variables[0][0].length-1); i++) {

                label = "Condition " + i;
                label_exp = "Condition " + i + " exp";
                labels_obs = labels_obs.concat(["Condition " + i, "Condition " + i + " exp"]);
                opts = {};
                opts['series'] = {};
                color = cp[cp.push(cp.shift())-1];
                opts['series'][label_exp] = {
                    color: color,
                    drawPoints: true,
                    connectSeparatedPoints: false,
                    strokeWidth: 0
                };
                opts['series'][label] = {
                    color: color,
                    drawPoints: false,
                    connectSeparatedPoints: true,
                };
                series_settings_obs = $.extend(series_settings_obs, opts['series']);
            };
            g_settings['observables'] = {
              connectSeparatedPoints: true,
              labels: labels_obs,
              errorBars: true,
              title: key,
              series: series_settings_obs
            };
        };

        if (plot === 'inputs'){
            g_settings['inputs'] = {
              labels: ['t', key],
              xlabel: ar['model.fu'][0],
              title: key
            };
      };
      for (var i = 0; i < $('.g_' + plot + '_' + key).length; i++) {
        graphs[plot][key][i] = new Dygraph(
          $('.g_' + plot + '_' + key)[i],
          graphs_data[plot][key],
          $.extend({}, g_settings[plot], g_settings[key], g_settings['global']));
      };
    };
  };
};

function select_model(tree) {
    var html = '<form><div class="ui-field-contain"><select id="select_model" name="select_model"><option>Select model...</option>';
    for (var key in tree['children']) {
        html += '<option value=\'' +  tree['children'][key]['children'][0]['name'] + '\'>' + tree['children'][key]['name'].replace(/^.*(\\|\/|\:)/, '') + '</option>';
    };
    html += '</select></div></form>';
    return html
};

function create_filetree(tree) {
    var html = '<ul>';
    for (var key in tree['children']) {
        if ('children' in tree['children'][key]) {
            html += '<li>' + tree['children'][key]['name'].replace(/^.*(\\|\/|\:)/, '') + '</li>';
            html += create_filetree(tree['children'][key]);
        } else {
            html += '<li><a name=\'' + tree['children'][key]['name'] + '\'  href=\'#\'>' + tree['children'][key]['name'].replace(/^.*(\\|\/|\:)/, '') + '</a></li>';
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
