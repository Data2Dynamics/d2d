import os
import re
import uuid
import time
import plotly
import numpy as np
from threading import Thread
from flask import Flask, jsonify, render_template, request, session, Markup
import pyd2d

__author__ = 'Clemens Blank'
app = Flask(__name__)

try:
    from flask.ext.compress import Compress
    app.config['COMPRESS_LEVEL'] = 4
    Compress(app)
except:
    pass

FIELDS = {
    'pLabel', 'p', 'lb', 'ub', 'chi2fit', 'config.optim.MaxIter',
    'config.nFinePoints', 'model.u', 'model.pu', 'model.py',
    'model.pystd', 'model.pv', 'model.pcond', 'model.fu',
    'model.name', 'model.xNames', 'model.z', 'model.description',
    'model.condition.tFine', 'model.condition.uFineSimu',
    'model.data.yNames', 'model.data.tFine', 'model.data.tExp',
    'model.plot.dLink', 'model.data.yFineSimu',
    'model.data.ystdFineSimu', 'model.data.yExp',
    'model.data.yExpStd', 'model.condition.xFineSimu',
    'model.condition.zFineSimu'
}

FIELDS_FIT = {
    'fit.iter', 'fit.iter_count', 'fit.improve', 'fit.chi2_hist',
    'fit.p_hist'
}

SESSION_LIFETIME = 3000
DEBUG = True

d2d_instances = {}


@app.before_request
def before_request():

    # create session
    if not session.get('uid'):
        session.permanent = True
        session['uid'] = uuid.uuid4()

    # update instance lifetime on any activity
    elif session.get('uid') in d2d_instances:
        d2d_instances[session.get('uid')]['alive'] = 1


@app.route('/_filetree', methods=['GET'])
def filetree():
    tree = create_filetree(path=os.path.join(os.getcwd(), 'models'),
                           max_depth=2)
    return jsonify(tree=tree)


@app.route('/', methods=['GET'])
@app.route('/d2d_presenter', methods=['GET'])
def d2d_presenter():

    return render_template('d2d_presenter.html')


@app.route('/_start', methods=['GET'])
def start():

    ROUND = int(request.args.get('round'))
    if ROUND is 0:
        ROUND = False

    nFinePoints = int(request.args.get('nFinePoints'))

    d2d_instances.update(
        {
            session['uid']: {'d2d': pyd2d.d2d(), 'alive': 1,
                             'model': os.path.join(
                                os.getcwd(),
                                request.args.get('model')),
                             'nFinePoints': nFinePoints,
                             'ROUND': ROUND,
                             'MODEL': 1,
                             'DSET': 1
                             }
        })

    d2d_instances[session['uid']]['d2d'].load_model(
        d2d_instances[session['uid']]['model'], load=1)

    d2d_instances[session['uid']]['d2d'].set(
        {'ar.config.nFinePoints': nFinePoints})

    d2d_instances[session['uid']]['d2d'].eval("arLink;")
    d2d_instances[session['uid']]['d2d'].simu()
    d2d_instances[session['uid']]['d2d'].eval("arChi2;")

    # thread will shut down the matlab instance on inactivity
    t = Thread(target=d2d_close_instance, args=(session['uid'], ))
    t.start()

    return jsonify({})


@app.route('/_update', methods=['GET'])
def update():

    if not session['uid'] in d2d_instances:
        return jsonify(session_expired=True)

    data = {}
    extra = {}

    d2d = d2d_instances[session['uid']]['d2d']

    if (request.args.get('filename') is None or
            request.args.get('filename') == "" or
            request.args.get('filename') == 'undefined'):
        extra['filename'] = os.path.join(
            d2d.path.replace(os.getcwd() +
                             os.path.sep, ''), d2d.filename)
    else:
        extra['filename'] = request.args.get('filename')

    options = request.args.get('options').split(';')

    if 'console' in options:
        command = request.args.get('command')
        d2d.eval(command)

        extra['console'] = d2d.output_total

    if 'chi2' in options:
        d2d.simu()
        d2d.eval("arChi2")

    if 'setup' in options:
        d2d.load_model(d2d_instances[session['uid']]['model'], load=None)

    if 'change_mdc' in options:

        if request.args.get('name') == 'MODEL':
            d2d_instances[session['uid']]['MODEL'] =\
                int(request.args.get('value'))
        elif request.args.get('name') == 'DSET':
            d2d_instances[session['uid']]['DSET'] =\
                int(request.args.get('value'))

    if 'simu_data' in options:
        d2d.simu('false', 'false', 'false')
        d2d.eval(
            'arSimuData(' +
            str(d2d_instances[session['uid']]['MODEL']) + ',' +
            str(d2d_instances[session['uid']]['DSET']) + '); arChi2;'
        )

    if 'model' in options:
        try:
            extra['svg'] = Markup(
                editor('read', os.path.join(
                        d2d.path, 'Models',
                        d2d.get(
                            {'model.name'},
                            d2d_instances[session['uid']]['MODEL'],
                            d2d_instances[session['uid']]['DSET']
                        ) + '.svg'))['content'])
        except:
            extra['svg'] = ''
        extra['size'] = {}

        try:
            extra['MODEL'] = d2d_instances[session['uid']]['MODEL']
            extra['size']['MODEL'] = d2d.eval('length(ar.model)', 1)
        except:
            pass
        try:
            extra['DSET'] = d2d_instances[session['uid']]['DSET']
            extra['size']['DSET'] = d2d.eval(
                "length(ar.model(" +
                str(d2d_instances[session['uid']]['MODEL']) +
                ").plot)", 1)
        except:
            pass

    if 'max_iter' in options:
        max_iter = int(request.args.get('max_iter'))
        d2d.set({'ar.config.optim.MaxIter': max_iter})

    if 'fit' in options:
        d2d.eval("arFit")

    if 'tree' in options:
        extra['tree'] = create_filetree(path=d2d.path)

    if 'read' in options:
        extra['editor_data'] = editor('read', extra['filename'])

    if 'write' in options:
        extra['editor_data'] = editor('write', extra['filename'],
                                      request.args.get('content'))

    if 'update_graphs' in options:  # set new parameters
        d2d.set_pars_from_dict(request.args.to_dict())

    if 'update_graphs' in options or 'create_graphs' in options:
        d2d.simu('false', 'true', 'false')
        data = select_data(d2d_instances[session['uid']], options)
        # converting data for use in plotly.js
        data = create_plotly_data(data)

    return jsonify(ar=data, extra=extra)


@app.route('/_console', methods=['GET'])
def console():
    if not session['uid'] in d2d_instances:
        return jsonify(session_expired=True)

    if request.args.get('command') != '':
        command = request.args.get('command')
        d2d_instances[session['uid']]['d2d'].eval(command)

    console = d2d_instances[session['uid']]['d2d'].output_total

    return jsonify(console=console)


def select_data(d2d_instance, options):

    global FIELDS, FIELDS_FIT

    if 'fit' in options:
        fields = FIELDS.copy()
        fields.update(FIELDS_FIT)
    else:
        fields = FIELDS

    return d2d_instance['d2d'].get(fields, d2d_instance['MODEL'],
                                   d2d_instance['DSET'],
                                   d2d_instance['ROUND'], True)


def create_plotly_data(data):

    data['plots'] = {}
    data['plots']['observables'] = []

    for i in range(len(data['model.data.yNames'][0])):
        data['plots']['observables'].append({'data': [], 'layout': dict(
            title=data['model.data.yNames'][0][i]
        )})

        for j in range(len(data['model.data.tExp'])):
            yFine_err = np.append(
                np.add(
                    data['model.data.yFineSimu'][j][i],
                    data['model.data.ystdFineSimu'][j][i]),
                np.subtract(
                    data['model.data.yFineSimu'][j][i],
                    data['model.data.ystdFineSimu'][j][i])[::-1])

            data['plots']['observables'][i]['data'].append(dict(
                name='Condition ' + str(j),
                x=np.append(
                    data['model.data.tFine'][j][0],
                    data['model.data.tFine'][j][0][::-1]),
                y=yFine_err,
                type='scatter',
                fill='tozerox',
                fillcolor='rgba(0,100,80,0.2)',
                line=dict(color='transparent'),
                showlegend=False,
                hoverinfo='none',
                )
            )
            data['plots']['observables'][i]['data'].append(dict(
                name='Condition ' + str(j),
                x=data['model.data.tFine'][j][0],
                y=data['model.data.yFineSimu'][j][i],
                type='scatter',
                mode='lines',
                )
            )
            # experimental data
            data['plots']['observables'][i]['data'].append(dict(
                name='Condition ' + str(j),
                x=data['model.data.tExp'][j][0],
                y=data['model.data.yExp'][j][i],
                type='scatter',
                showlegend=False,
                mode='markers',
                error_y=dict(
                    type='data',
                    array=data['model.data.yExpStd'][j][i],
                    visible=True,
                    thickness=1,
                    width=5
                    )
                )
            )

    data['plots']['variables'] = []

    if len(data['model.z']) > 0:

        data['model.xNames'] = np.append(
            data['model.xNames'], data['model.z'])

        data['model.condition.xFineSimu'] = np.append(
            data['model.condition.xFineSimu'],
            data['model.condition.zFineSimu'], axis=1)

    for i in range(len(data['model.xNames'])):
        data['plots']['variables'].append({'data': [], 'layout': dict(
            title=data['model.xNames'][i]
        )})
        for j in range(len(data['model.condition.tFine'])):

            data['plots']['variables'][i]['data'].append(dict(
                name='Condition ' + str(j),
                x=data['model.condition.tFine'][j][0],
                y=data['model.condition.xFineSimu'][j][i],
                type='scatter',
                )
            )

    data['plots']['inputs'] = []

    if len(data['model.u']) > 0:
        if isinstance(data['model.u'][0], str):
            data['model.condition.uFineSimu'] =\
                [data['model.condition.uFineSimu']]
        for i in range(len(data['model.u'])):
            data['plots']['inputs'].append({'data': [], 'layout': dict(
                title=data['model.u'][i]
            )})

            for j in range(len(data['model.condition.uFineSimu'])):

                data['plots']['inputs'][i]['data'].append(dict(
                    name='Condition ' + str(j),
                    x=data['model.condition.tFine'][j][0],
                    y=data['model.condition.uFineSimu'][j][0][0],
                    type='scatter',
                    )
                )

    # remove unecessary data to lower the traffic
    data.pop('model.data.tFine', None)
    data.pop('model.data.yFineSimu', None)
    data.pop('model.data.ystdFineSimu', None)
    data.pop('model.data.tExp', None)
    data.pop('model.data.yExp', None)
    data.pop('model.data.yExpStd', None)
    data.pop('model.condition.tFine', None)
    data.pop('model.condition.uFineSimu', None)
    data.pop('model.condition.xFineSimu', None)
    data.pop('model.condition.zFineSimu', None)
    data.pop('model.condition.z', None)

    return data


def editor(option, filename, content=None):

    file_content = {}

    if option == 'read':
        try:
            file = open(filename, 'r')
            file_content = {"content": file.read()}
            file.close()
        except:
            print("Couldn't read file.")

    elif option == 'write':

        try:
            file = open(filename, 'w')
            file.write(content)
            file.close()
            file_content = {"content": content}

        except:
            print("Couldn't write file.")

    return file_content


def create_filetree(path=None, depth=0, max_depth=0):

    tree = None

    if max_depth == 0 or depth < max_depth:
        if path is None:
            path = os.getcwd()

        tree = dict(name=os.path.basename(path), children=[])

        try:
            lst = os.listdir(path)
        except OSError:
            pass  # ignore errors
        else:
            for name in lst:
                fn = os.path.join(path, name)
                if (os.path.isdir(fn) and
                        re.match('^.*(Compiled)$', fn) is None):
                    child = create_filetree(fn, depth + 1, max_depth)
                    if child is not None:
                        tree['children'].append(child)
                elif re.match('^.*\.(m|def|txt|csv)$', fn) is not None:
                    tree['children'].append(dict(name=fn.replace(
                        os.getcwd() + os.path.sep, "")))

    return tree


def d2d_close_instance(uid):

    while True:
        try:
            if d2d_instances[uid]['alive'] == 0:
                # clean up the global dicts, deletes the instances
                del(d2d_instances[uid])
                break

            d2d_instances[uid]['alive'] = 0
            time.sleep(SESSION_LIFETIME)
        except:
            print('Unable to shutdown thread ' + str(uid))
            break


if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.permanent_session_lifetime = SESSION_LIFETIME
    app.debug = DEBUG
    app.json_encoder = plotly.utils.PlotlyJSONEncoder
    app.run(threaded=True,
            host="0.0.0.0",
            port=int("5000")
            )
