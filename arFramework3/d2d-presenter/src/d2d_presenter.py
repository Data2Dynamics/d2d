import os
import re
import uuid
import time
import numpy
from threading import Thread
from flask import Flask, jsonify, render_template, request, session, Markup
from simplejson import JSONEncoder
import pyd2d
import logging
logger = logging.getLogger('ftpuploader')

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
                        )['model.name'] + '.svg'))['content'])
        except:
            extra['svg'] = ''

        extra['size'] = {}

        try:
            extra['MODEL'] = d2d_instances[session['uid']]['MODEL']
            extra['size']['MODEL'] = d2d.eval('length(ar.model)', 1)
            extra['size']['MODELNAMES'] = []

            for i in range(int(extra['size']['MODEL'])):
                extra['size']['MODELNAMES'].append(d2d.get(
                    {'model.name'},
                    i+1,
                    d2d_instances[session['uid']]['DSET'], False, 'list'
                )['model.name'])
        except:
            pass
        try:
            extra['DSET'] = d2d_instances[session['uid']]['DSET']
            extra['size']['DSET'] = d2d.eval(
                "length(ar.model(" +
                str(d2d_instances[session['uid']]['MODEL']) +
                ").plot)", 1)
            extra['size'].update(d2d.get(
                    {'model.plot.name'},
                    d2d_instances[session['uid']]['MODEL'],
                    d2d_instances[session['uid']]['DSET'], False, 'list'
                ))
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

    if 'chi2' in options:
        d2d.simu()
        d2d.eval("arChi2")

    if 'update_graphs' in options or 'create_graphs' in options:
        d2d.simu('false', 'true', 'false', options)
        data = select_data(d2d_instances[session['uid']], options)

        data = create_dygraphs_data(data)

        d2d_instances[session['uid']]['data'] = data

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
                                   d2d_instance['ROUND'], 'list')

def create_dygraphs_data(data):
    """Prepares the selected data for use in dygraphs.js.
        Converts float('nan') to dygraph compatible values:
        "'NaN'" creates an actual gap in the data, "None" still allows to
        connect the seperated points.
        Errors need to be set to 0 if not available.
        Sometimes uFineSimu contains float('inf') (javascript: infinity)
        which does not work in dygraphs - set to None instead.

        Finally combibnes all plot data into a dygraphs compatible data set.
        """

    # make sure there are no two data entries with same timestamp, and if so
    # add epsilon to it, since dygraphs might get confused
    # for j in range(len(data['model.data.tExp'])):
    #     for i, key in enumerate(data['model.data.tExp'][j][0]):
    #         while (data['model.data.tExp'][j][0][i] in
    #                data['model.data.tFine'][j][0]):
    #             print(data['model.data.tExp'][j][0][i])
    #             data['model.data.tExp'][j][0][i] =\
    #              data['model.data.tExp'][j][0][i] + sys.float_info.epsilon

    data['plots'] = {}
    data['plots']['observables'] = []

    # for j in range(len(data['model.data.yNames'][0])):
    #     temp_fine = []
    #     temp_exp = []
    #     temp.append(data['model.data.tFine'][0][0])
    #
    #     for i in range(len(data['model.data.yExp'])):
    #         temp_exp.append(list(zip(data['model.data.yFineSimu'][i][j], data['model.data.ystdFineSimu'][i][j])))
    #         temp.append(list(zip(data['model.data.yFineSimu'][i][j], data['model.data.ystdFineSimu'][i][j])))
    #
    #     data['plots']['observables'].append(list(zip(*temp)))
    #
    #
    #
    #     for j in range(data['model.data.tFine']):
    #         zip(data['model.data.yExp'][j][0], data['model.yExpStd'][j][0])
    #         zip(data['model.data.yFineSimu'], data['model.yExpStd'])

    for i in range(len(data['model.data.yNames'][0])):
        data['plots']['observables'].append(
            data['model.data.tFine'][0][0].copy() +
                data['model.data.tExp'][0][0].copy())
        for j in range(len(data['model.data.tFine'][0][0])):
            data['plots']['observables'][i][j] = [data['plots']['observables'][i][j]]
            for k in range(len(data['model.data.yFineSimu'])):
                data['plots']['observables'][i][j].append([
                    data['model.data.yFineSimu'][k][i][j],
                    data['model.data.ystdFineSimu'][k][i][j]])
                data['plots']['observables'][i][j].append([None, 0])

        for j in range(len(data['model.data.tExp'][0][0])):
            c = j + len(data['model.data.tFine'][0][0])
            data['plots']['observables'][i][c] = [data['plots']['observables'][i][c]]
            for k in range(len(data['model.data.yExp'])):
                data['plots']['observables'][i][c].append([None, None])
                data['plots']['observables'][i][c].append(
                    [data['model.data.yExp'][k][i][j],
                     data['model.data.yExpStd'][k][i][j]])

        data['plots']['observables'][i].sort(key=lambda x: x[0])

    if len(data['model.z']) > 0:

        data['model.xNames'] = data['model.xNames'] + data['model.z']

        xzFine = []

        for i in range(len(data['model.condition.xFineSimu'])):

            xzFine.append(([data['model.condition.xFineSimu'][i] +
                                 data['model.condition.zFineSimu'][i]])[0])

        data['model.condition.xFineSimu'] = xzFine

    data['plots']['variables'] = []

    for i in range(len(data['model.xNames'])):
        data['plots']['variables'].append(data['model.condition.tFine'][0][0].copy())
        for j in range(len(data['model.condition.tFine'][0][0])):
            data['plots']['variables'][i][j] = [data['plots']['variables'][i][j]]
            for k in range(len(data['model.condition.xFineSimu'])):
                data['plots']['variables'][i][j].append(
                    data['model.condition.xFineSimu'][k][i][j])

    data['plots']['inputs'] = []

    if len(data['model.u']) > 0:
        if isinstance(data['model.u'][0], str):
            data['model.condition.uFineSimu'] = [data['model.condition.uFineSimu']]

        for i in range(len(data['model.u'])):
            data['plots']['inputs'].append(data['model.condition.tFine'][0][0].copy())
            for j in range(len(data['model.condition.tFine'][0][0])):
                data['plots']['inputs'][i][j] = [data['plots']['inputs'][i][j]]
                data['plots']['inputs'][i][j].append(
                    data['model.condition.uFineSimu'][0][0][i][j])

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

def select_data_old(d2d, options):

    ar = d2d.dic2obj()

    data = {}

    data['model'] = {}
    data['data'] = {}
    data['condition'] = {}
    data['config'] = {}
    data['fit'] = {}
    data['plot'] = {}
    data['plots'] = {}

    try:
        data['pLabel'] = ar.pLabel
        data['p'] = ar.p
        data['lb'] = ar.lb
        data['ub'] = ar.ub
        data['chi2fit'] = ar.chi2fit
    except AttributeError:
        pass

    try:
        data['config']['MaxIter'] = ar.config[0].optim[0].MaxIter
        data['config']['nFinePoints'] = ar.config[0].nFinePoints
    except AttributeError:
        pass

    if 'fit' in options:
        try:
            data['fit']['iter'] = ar.fit[0].iter
            data['fit']['iter_count'] = ar.fit[0].iter_count
            data['fit']['improve'] = ar.fit[0].improve
            data['fit']['p_hist'] = ar.fit[0].p_hist
            data['fit']['chi2_hist'] = ar.fit[0].chi2_hist
            if type(data['fit']['chi2_hist']) == float:
                data['fit']['chi2_hist'] = [[data['fit']['chi2_hist'], 'NaN']]
        except:
            pass

    try:
        data['model']['u'] = ar.model[M].u  # Input labels
        data['model']['pu'] = ar.model[M].pu
        data['model']['py'] = ar.model[M].py
        data['model']['pystd'] = ar.model[M].pystd
        data['model']['pv'] = ar.model[M].pv
        data['model']['pcond'] = ar.model[M].pcond
        data['model']['fu'] = ar.model[M].fu
        data['model']['name'] = ar.model[M].name
        data['model']['xNames'] = ar.model[M].xNames
        data['model']['z'] = ar.model[M].z  # derived variables labels
        data['model']['description'] = ar.model[M].description
    except AttributeError:
        pass

    try:
        data['condition']['tFine'] = ar.model[M].condition[0].tFine
        data['condition']['uFineSimu'] = \
            ar.model[M].condition[C].uFineSimu  # Inputs
    except AttributeError:
        pass

    try:
        data['data']['yNames'] = ar.model[M].data[D].yNames
        data['data']['tFine'] = ar.model[M].data[0].tFine
        data['data']['tExp'] = ar.model[M].data[0].tExp
    except AttributeError:
        pass

    data['data']['yFineSimu'] = []
    data['data']['ystdFineSimu'] = []
    data['data']['yExp'] = []
    data['data']['yExpStd'] = []
    data['condition']['xFineSimu'] = []
    data['condition']['zFineSimu'] = []

    if isinstance(ar.model[M].plot[D].dLink, (float, int)):
        data['plot']['dLink'] = [[ar.model[M].plot[D].dLink]]
    else:
        data['plot']['dLink'] = ar.model[M].plot[D].dLink


    for h, i in enumerate(data['plot']['dLink'][0]):
        i = int(i)
        data['data']['yFineSimu'].append(ar.model[M].data[i-1].yFineSimu)
        data['data']['ystdFineSimu'].append(ar.model[M].data[i-1].ystdFineSimu)
        data['data']['yExp'].append(ar.model[M].data[i-1].yExp)
        data['data']['yExpStd'].append(ar.model[M].data[i-1].yExpStd)
        data['condition']['xFineSimu'].append(ar.model[M].condition[h].xFineSimu)
        data['condition']['zFineSimu'].append(ar.model[M].condition[h].zFineSimu)

    return data

def create_dygraphs_data_old(data):
    """Prepares the selected data for use in dygraphs.js.
        Converts float('nan') to dygraph compatible values:
        "'NaN'" creates an actual gap in the data, "None" still allows to
        connect the seperated points.
        Errors need to be set to 0 if not available.
        Sometimes uFineSimu contains float('inf') (javascript: infinity)
        which does not work in dygraphs - set to None instead.

        Finally combibnes all plot data into a dygraphs compatible data set.
        """


    # make sure there are no two data entries with same timestamp, and if so
    # add epsilon to it, since dygraphs might get confused
    for i, key in enumerate(data['data']['tExp']):
        while data['data']['tExp'][i] in data['data']['tFine']:
            data['data']['tExp'][i][0] = numpy.nextafter(
                data['data']['tExp'][i][0], abs(
                    2*data['data']['tExp'][i][0] + 1))

    data['plots']['observables'] = []

    for i in range(len(data['data']['yNames'][0])):
        data['plots']['observables'].append(
            data['data']['tFine'] + data['data']['tExp'])
        for j in range(len(data['data']['tFine'])):
            for k in range(len(data['plot']['dLink'][0])):
                data['plots']['observables'][i][j] = \
                    data['plots']['observables'][i][j] + [[
                        data['data']['yFineSimu'][k][j][i],
                        data['data']['ystdFineSimu'][k][j][i]], [None, 0]]

        for j in range(len(data['data']['tExp'])):
            for k in range(len(data['plot']['dLink'][0])):
                c = j + len(data['data']['tFine'])
                data['plots']['observables'][i][c] = \
                    data['plots']['observables'][i][c] + [
                    [None, None], [data['data']['yExp'][k][j][i],
                                   data['data']['yExpStd'][k][j][i]]]

        data['plots']['observables'][i].sort(key=lambda x: x[0])


    if len(data['model']['z']) > 0:

        data['model']['xNames'] = [data['model']['xNames'][0] +
                                   data['model']['z'][0]]

        xzFine = []

        for i in range(len(data['plot']['dLink'][0])):
            xzFine.append([])
            for j in range(len(data['condition']['tFine'])):
                xzFine[i].append(data['condition']['xFineSimu'][i][j] +
                                 data['condition']['zFineSimu'][i][j])

        data['condition']['xFineSimu'] = xzFine

    data['plots']['variables'] = []

    for i in range(len(data['model']['xNames'][0])):
        data['plots']['variables'].append(data['condition']['tFine'].copy())
        for j in range(len(data['condition']['tFine'])):
            for k in range(len(data['plot']['dLink'][0])):
                data['plots']['variables'][i][j] = \
                    data['plots']['variables'][i][j] + [data['condition']['xFineSimu'][k][j][i]]

    data['plots']['inputs'] = []

    if len(data['model']['u']) > 0:
        if isinstance(data['model']['u'][0], str):
            data['condition']['uFineSimu'] = [data['condition']['uFineSimu']]
        else:
            data['condition']['uFineSimu'] = data['condition']['uFineSimu']

        if len(data['model']['u']) > 0:
            for i in range(len(data['model']['u'][0])):
                data['plots']['inputs'].append(data['condition']['tFine'].copy())
                for j in range(len(data['condition']['tFine'])):
                    data['plots']['inputs'][i][j] = \
                        data['plots']['inputs'][i][j] + \
                        [data['condition']['uFineSimu'][j][i]]

    # remove unecessary data to lower the traffic
    data['data'].pop('tFine', None)
    data['data'].pop('yFineSimu', None)
    data['data'].pop('ystdFineSimu', None)
    data['data'].pop('tExp', None)
    data['data'].pop('yExp', None)
    data['data'].pop('yExpStd', None)
    data['condition'].pop('tFine', None)
    data['condition'].pop('uFineSimu', None)
    data['condition'].pop('xFineSimu', None)
    data['condition'].pop('zFineSimu', None)
    data['condition'].pop('z', None)

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
                if (os.path.isdir(fn) and re.match('^.*(Compiled)$', fn)
                    is None):
                    child = create_filetree(fn, depth + 1, max_depth)
                    if child != None:
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


class JSONEncoder_ignore(JSONEncoder):

    """ set the encoder to ignore nan values, which will transfere all 'NaN',
    infinity etc. to "null" - somehing dygraph understands. """

    def __init__(self, **kwargs):
        """Leaves JSONEncoder as it is, just switches "ingore_nan" on."""
        kwargs['ignore_nan'] = True
        super(JSONEncoder_ignore, self).__init__(**kwargs)

if __name__ == "__main__":
    app.secret_key = os.urandom(24)
    app.permanent_session_lifetime = SESSION_LIFETIME
    app.debug = DEBUG
    app.json_encoder = JSONEncoder_ignore
    app.run(threaded=True,
            host="127.0.0.1",
            port=int("5000")
            )
