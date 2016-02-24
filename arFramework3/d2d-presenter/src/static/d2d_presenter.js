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
            'create_fit_data',
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
                    this.complete = 1;
                    this.success(options, event);
                    $.mobile.loading('hide');
                },
            });
        };
    };
    this.success = function(options, event) {
        if (options.indexOf('create_fit_data') > -1) {
            create_fit_data(['fit']);
        };
        if (options.indexOf('start') > -1) {
            create_content();
            create_events();
            create_graphs();
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
                create_fit_data(['fit']);
            } else {
                chi2_hist = ar['fit.chi2_hist'];
                p_hist = ar['fit.p_hist'];
            };
        };
        if (options.indexOf('fit_auto') > -1) {
            chi2 = chi2_hist.shift()[0];
            p = [];
            for (i in p_hist) {
                p.push(p_hist[i].shift());
            };
            if ((chi2 !== null)) {
                url_data = {};
                for (var key in ar['pLabel']) {
                    url_data[ar['pLabel'][key]] = p[key];
                };
                create_fit_data(['fit', 'fit_auto',]);
                this.update([
                    'fit_auto', 'set_sliders', 'update_graphs',
                ], url_data);
            };
        };
        if (options.indexOf('fit_clear') > -1) {
            fit_data = {};
            create_fit_data(['fit']);
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
        if (options.indexOf('resize_graphs') > -1) {
            for (key in ar['plots']) {
                for (var i in ar['plots'][key]) {
                    elm = $('.g_' + key + '_' + ar['plots'][key][i]['layout']['title']);
                    Plotly.Plots.resize(elm[0]);
                };
            };
            for (var i in ar['plots']['observables']) {
                elm = $('.g_observables_fit_' + ar['plots']['observables'][i]['layout']['title']);
                Plotly.Plots.resize(elm[0]);
            };
            Plotly.Plots.resize($('.g_fit_chi2plot')[0]);
            Plotly.Plots.resize($('.g_fit_parameters')[0]);
        };
        if (options.indexOf('update_graphs') > -1) {
            // check that it wont update the graphs after full fit when auto fitting
            if ((event !== null && event.target.name !== 'fit_auto') || (event === null)) {
                if ($("ul li.ui-state-active")[0].children[0].id === 'tab_fit') {
                    $('.g_fit_chi2plot')[0].data[0]['x'] = fit_data.chi2plot.x;
                    $('.g_fit_chi2plot')[0].data[0]['y'] = fit_data.chi2plot.y;
                    Plotly.redraw($('.g_fit_chi2plot')[0]);
                    for (var i in ar['p']) {
                        $('.g_fit_parameters')[0].data[i]['x'] = fit_data.chi2plot.x;
                        $('.g_fit_parameters')[0].data[i]['y'] = fit_data.parameters.y[i];
                    };
                    Plotly.redraw($('.g_fit_parameters')[0]);
                    for (var i in ar['plots']['observables']) {
                        elm = $('.g_observables_fit_' + ar['plots']['observables'][i]['layout']['title']);
                        elm[0].data = ar['plots']['observables'][i]['data'];
                        Plotly.redraw(elm[0]);
                    };
                } else if ($("ul li.ui-state-active")[0].children[0].id === 'tab_plots') {
                    if ($('#inputs')[0].checked) {
                        for (var i in ar['plots']['inputs']) {
                            elm = $('.g_inputs_' + ar['plots']['inputs'][i]['layout']['title']);
                            elm[0].data[0].y = ar['plots']['inputs'][i]['data'][0]['y'];
                            Plotly.redraw(elm[0]);
                        };
                    };
                    if ($('#variables')[0].checked) {
                        for (var i in ar['plots']['variables']) {
                            elm = $('.g_variables_' + ar['plots']['variables'][i]['layout']['title']);
                            elm[0].data = ar['plots']['variables'][i]['data'];
                            Plotly.redraw(elm[0]);
                        };
                    };
                    if ($('#observables')[0].checked) {
                        for (var i in ar['plots']['observables']) {
                            elm = $('.g_observables_' + ar['plots']['observables'][i]['layout']['title']);
                            elm[0].data = ar['plots']['observables'][i]['data'];
                            Plotly.redraw(elm[0]);
                        };
                    };
                };
            };
        };
        if (options.indexOf('create_graphs') > -1) {
            $('#' + ar['pLabel'][key] + '_slider').slider('refresh');
            create_graphs();
        };
    };
};
function create_content() {
    for (var key in extra.size) {
        var name;
        if (key === 'MODEL')
            name = 'Model';
        if (key === 'DSET')
            name = 'Data';

        var div = $('<div />', {
            'class': 'ui-field-contain',
            text: '',
        });
        var select = $('<select />', {
            name: key,
            class: 'select_mdc',
        });
        for (var i = 1; i <= extra.size[key]; i++) {
            if (i === extra[key]) {
                var option = $('<option />', {
                    value: i,
                    text: name + ' ' + i,
                    selected: '',
                });
            } else {
                var option = $('<option />', {
                    value: i,
                    text: name + ' ' + i,
                });
            };
            select.append(option)
        };
        $('#' + key).html(div.append(select));
    };
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
            'class': 'g_container_items g_observables_fit_' + ar['model.data.yNames'][0][key],
            id: 'g_observables_fit_' + ar['model.data.yNames'][0][key],
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
            step: 0.1,
            text: ar['pLabel'][key],
        });
        if (pv.indexOf(ar['pLabel'][key]) > -1) {
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
            'change_mdc', 'model', 'create_graphs',
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
                update.update([
                    'resize_graphs', 'update_graphs',
                ], url_data, e);
            };
        });
    };
    $('#max_iter').on('slidestop', function(e) {
        update.update(['max_iter'], {
            max_iter: e.target.value
        }, e);
    });
    $("input[name='toggle_graphs']").change(function(e) {
        $('#graphs_' + e.target.id).toggle();
        if ($('#graphs_' + e.target.id).is(":visible")) {
            for (key in ar['plots']) {
                for (var i in ar['plots'][key]) {
                    elm = $('.g_' + key + '_' + ar['plots'][key][i]['layout']['title']);
                    Plotly.Plots.resize(elm[0]);
                };
            };
        };
    });
    $(".panel_button").on('click', function(e) {
        $('.panel').toggle();
    });
    $("button[name='simu_data']").on('click', function(e) {
        update.update([
            'create_fit_data', 'simu_data', 'update_graphs',
        ], {}, e);
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
            'create_fit_data', 'set_sliders', 'update_graphs',
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
        update.update([
            'update_graphs', 'model',
        ], {}, e);
    });
    $(window).on('resize', function(e) {
        update.update(['resize_graphs'], {}, e);
    });
    $('#tab_plots').click(function(e) {
        $.mobile.loading('show');
        update.update([
            'set_sliders', 'resize_graphs',
        ], {}, e);
    });
    $('#tab_fit').click(function(e) {
        $.mobile.loading('show');
        update.update([
            'chi2', 'resize_graphs',
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
function create_fit_data(options) {
    if (!fit_data.hasOwnProperty('parameters')) {
        fit_data.parameters = {};
        fit_data.parameters.y = [];
        for (var i in ar['p']) {
            fit_data.parameters.y.push([ar['p'][i][0]]);
        };
        fit_data.chi2plot = {};
        fit_data.chi2plot.x = [0];
        fit_data.chi2plot.y = [ar['chi2fit']];
    } else {
        if (options.indexOf('fit') > -1 || options.indexOf('create_graphs') > -1) {
            var fit_step = 1;
            if (options.indexOf('fit_auto') === -1) {
                fit_step = ar['fit.iter_count'] || 1;
            };
            fit_data.chi2plot.x.push(fit_data.chi2plot.x[fit_data.chi2plot.x.length - 1] + fit_step);
            fit_data.chi2plot.y.push(ar['chi2fit']);
            for (var i in ar['p']) {
                fit_data.parameters.y[i].push(ar['p'][i][0]);
            };
        };
    };
    return fit_data
};
function create_graphs() {
    labels_fit_parameters = ['n'];
    for (var key in ar['pLabel']) {
        labels_fit_parameters.push(ar['pLabel'][key]);
    };
    Plotly.newPlot('g_fit_chi2plot', [
        {
            x: fit_data.chi2plot.x,
            y: fit_data.chi2plot.y,
        }
    ], {
        title: '-2*log(L)'
    }, {scrollZoom: true});
    param_data = [];
    for (var i in ar['p']) {
        param_data.push({name: ar['pLabel'][i], x: fit_data.chi2plot.x, y: fit_data.parameters.y[i],})
    };
    Plotly.newPlot('g_fit_parameters', param_data, {
        title: 'Fit parameters'
    }, {scrollZoom: true});
    for (key in ar['plots']) {
        for (var i in ar['plots'][key]) {
            elm = $('.g_' + key + '_' + ar['plots'][key][i]['layout']['title']);
            Plotly.newPlot(elm[0], // the ID of the div
                    ar['plots'][key][i]['data'], ar['plots'][key][i]['layout'], {scrollZoom: true});
        };
    };
    for (var i in ar['plots']['observables']) {
        elm = $('.g_observables_fit_' + ar['plots']['observables'][i]['layout']['title']);
        Plotly.newPlot(elm[0], // the ID of the div
                ar['plots']['observables'][i]['data'], ar['plots']['observables'][i]['layout'], {scrollZoom: true});
    };
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
