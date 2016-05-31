"""The D2D Presenter main application.

   Includes the flask web server and the the modifications to run with
   dygraphs.js.

   The options are set in this file after the python imports:

   SESSION_LIFETIME: Positive integer
                     The timespan in seconds, after which the client session
                     shall be terminated when no user input happened. Get's
                     resetted everytime the server recieves a request.

   COPY:             False|True
                     Decides if the working directory per client session should
                     be a temporary copy, which will be deleted after the
                     session ends, or if the clients will work with the
                     original model directories.

   HIDE_CONSOLE:     False|True
                     Decides if the input console to execute Matlab
                     commands shall be accesible. For public use this might not
                     be covered by the Matlab licence.

   DEBUG:            True|False
                     Sets the mode of the Flask webserver. Should be set to
                     "False" when used in a productive environment.

   HOST:             IP Adress
                     Sets the IP adress of the Flask webserver. For external
                     access set to "0.0.0.0", for internal usage on the same
                     server use "127.0.0.1".

   PORT:             Positive integer
                     Sets the port, on which the Flask webserver should listen.

   COMPRESS_LEVEL:   Integer between 0 and 9
                     If the Flask extension Flaks-compress (flask.ext.compress)
                     is available, sets the level of compression, where 0 is
                     no compression and 9 is maximum compression (highest CPU
                     usage). Highly recommended since it decreases the amount
                     of traffic immensiveley. Default is 5.

"""
import os
import re
import uuid
import time
import copy
import stat
from distutils.dir_util import copy_tree, remove_tree
from threading import Thread
from flask import Flask, jsonify, render_template, request, session, Markup
from flask import redirect
from simplejson import JSONEncoder
import pyd2d

# Configuration #
SESSION_LIFETIME = 600
DEBUG = True
COPY = False
HIDE_CONSOLE = False
HOST = '127.0.0.1'
PORT = '5000'
COMPRESS_LEVEL = 5

__author__ = 'Clemens Blank'
app = Flask(__name__)

try:
    from flask.ext.compress import Compress
    app.config['COMPRESS_LEVEL'] = COMPRESS_LEVEL
    Compress(app)
except:
    pass

FIELDS = {
    'pLabel', 'p', 'lb', 'ub', 'chi2fit', 'config.optim.MaxIter',
    'config.nFinePoints', 'model.u', 'model.pu', 'model.py', 'model.pc',
    'model.pystd', 'model.pv', 'model.pcond', 'model.fu', 'model.data.pystd',
    'model.data.pcond', 'model.data.pu', 'model.data.py',
    'model.name', 'model.xNames', 'model.z', 'model.description',
    'model.condition.tFine', 'model.condition.uFineSimu',
    'model.data.yNames', 'model.data.tFine', 'model.data.tExp',
    'model.data.yFineSimu', 'model.data.ystdFineSimu', 'model.data.yExp',
    'model.data.yExpStd', 'model.condition.xFineSimu',
    'model.condition.zFineSimu'
}

FIELDS_FIT = {
    'fit.iter', 'fit.iter_count', 'fit.improve', 'fit.chi2_hist',
    'fit.p_hist'
}

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

    return render_template('d2d_presenter.html', hide_console=HIDE_CONSOLE)


@app.route('/_start', methods=['GET'])
def start():
    # delete already existing matlab engine for this user
    try:
        del(d2d_instances[session['uid']])
    except:
        pass

    ROUND = int(request.args.get('round'))
    if ROUND is 0:
        ROUND = False

    status = {}
    status['nFinePoints_min'] = False
    status['arSimu'] = False

    nFinePoints = int(request.args.get('nFinePoints'))
    do_compile = str(request.args.get('compile'))
    # Set the directory you want to copy from
    rootDir = os.path.dirname(os.path.join(
       os.getcwd(),
       request.args.get('model')))

    # save the origin of the copied model
    originDir = os.path.join(os.getcwd(), request.args.get('model'))

    if COPY is True:
        temp_dir = os.path.join(os.getcwd(), 'temp', str(session['uid']))
        try:
            for root, dirs, files in os.walk(temp_dir):
                for fname in files:
                    full_path = os.path.join(root, fname)
                    os.chmod(full_path, stat.S_IWRITE)
            remove_tree(temp_dir)
        except:
            pass
        copy_tree(rootDir, temp_dir)
        rootDir = os.path.join(temp_dir,
                               os.path.basename(request.args.get('model')))
    else:
        rootDir = originDir

    d2d_instances.update(
        {
            session['uid']: {'d2d': pyd2d.d2d(), 'alive': 1,
                             'model': rootDir,
                             'origin': originDir,
                             'nFinePoints': nFinePoints,
                             'ROUND': ROUND,
                             'MODEL': 1,
                             'DSET': 1
                             }
        })

    load = None

    if not do_compile.endswith('on'):
        try:
            results = os.listdir(os.path.join(
                os.path.dirname(d2d_instances[session['uid']]['model']),
                'Results'))

            for savename in results:
                if savename.endswith('_d2d_presenter'):
                    load = savename
                    print("A result has been found and will be loaded.")
                    break
            if load is None:
                print("No saved results found, compiling " +
                      "from scratch and save the result. This might take " +
                      "some time. Keep in mind that arSave will only be " +
                      "persistent if COPY is False.")
        except:
            pass
    else:
        load = False
        print("The model will be compiled. This might take some time.")

    d2d_instances[session['uid']]['d2d'].load_model(
        d2d_instances[session['uid']]['model'], load=load,
        origin=d2d_instances[session['uid']]['origin'])

    try:
        nFinePoints_min = d2d_instances[session['uid']]['d2d'].get(
            {'d2d_presenter.nFinePoints_min'},
            1, 1, False, 'list')['d2d_presenter.nFinePoints_min'][0][0]
    except:
        nFinePoints_min = nFinePoints

    if nFinePoints_min > nFinePoints:
        nFinePoints = nFinePoints_min
        status['nFinePoints_min'] = nFinePoints

    d2d_instances[session['uid']]['d2d'].set(
        {'ar.config.nFinePoints': nFinePoints})
    d2d_instances[session['uid']]['d2d'].eval("arLink;")
    status['arSimu'] = d2d_instances[session['uid']]['d2d'].simu()
    d2d_instances[session['uid']]['d2d'].eval("arChi2;")

    # thread will shut down the matlab instance on inactivity
    t = Thread(target=d2d_close_instance, args=(session['uid'], ))
    t.start()

    if str(request.args.get('direct_access')).endswith('true'):
        return redirect("#main_page")
    else:
        return jsonify(status=status)


@app.route('/_update', methods=['GET'])
def update():

    if not session['uid'] in d2d_instances:
        return jsonify(session_expired=True)

    data = {}
    extra = {}
    status = {}

    extra['HIDE_CONSOLE'] = HIDE_CONSOLE

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

        if HIDE_CONSOLE is False:
            command = request.args.get('command')
            d2d.eval(command)

        extra['console'] = d2d.output_total

    if 'setup' in options:
        d2d.load_model(d2d_instances[session['uid']]['model'], load=False)

    if 'change_mdc' in options:

        if request.args.get('name') == 'MODEL':
            d2d_instances[session['uid']]['MODEL'] =\
                int(request.args.get('value'))
        elif request.args.get('name') == 'DSET':
            d2d_instances[session['uid']]['DSET'] =\
                int(request.args.get('value'))

    if 'simu_data' in options:
        status['arSimu'] = d2d.simu('false', 'false', 'false')
        status['arSimuData'] = d2d.eval(
            'arSimuData(' +
            str(d2d_instances[session['uid']]['MODEL']) + ',' +
            str(d2d_instances[session['uid']]['DSET']) + ');'
        )
        d2d.eval('arChi2;')

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
        extra['editor_data'] = editor(d2d.path, 'read', extra['filename'])

    if 'write' in options:
        extra['editor_data'] = editor(d2d.path, 'write', extra['filename'],
                                      request.args.get('content'))

    if 'update_graphs' in options:  # set new parameters
        d2d.set_pars_from_dict(request.args.to_dict())

    if 'chi2' in options:
        status['arSimu'] = d2d.simu()
        d2d.eval("arChi2")

    if 'update_graphs' in options or 'create_graphs' in options:
        status['arSimu'] = d2d.simu('false', 'true', 'false')
        data = select_data(d2d_instances[session['uid']], options)

        data = create_dygraphs_data(data)

        d2d_instances[session['uid']]['data'] = data

    return jsonify(ar=data, extra=extra, status=status)


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
    """Reads out the specified ar fields.
    MODEL set's the model from which to select the data.
    DSET set's the dataset from which to select the data.
    ROUND specifies how many significant digit should be respected (in order
    to decrease the amount of data sent to the client), however it comes at
    the cost of some ms (e.g. 10ms in Raia_CancerResearch2011 on medium
    hardware).
    'list' can be used to get python lists, 'numpy' to get numpy array (2ms
    faster).
    If option 'fit' is set, also the fields in FIELDS_FIT will be respected.
    """
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
    """Converts the data for use in dygraphs.js.

    Errors need to be set to 0 if not available.
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

    if data['model.data.yNames'][0]:

        tObs = data['model.data.tFine'][0][0].copy()

        if data['model.data.tExp'][0]:
            for i in range(len(data['model.data.tExp'])):
                tObs = tObs + data['model.data.tExp'][i][0].copy()
        tObs = [[x] for x in tObs]

        for i in range(len(data['model.data.yNames'][0])):
            data['plots']['observables'].append(copy.deepcopy(tObs))
            for k in range(len(data['model.data.tFine'])):
                for j in range(len(data['model.data.tFine'][k][0])):
                    data['plots']['observables'][i][j].append([
                        data['model.data.yFineSimu'][k][i][j],
                        data['model.data.ystdFineSimu'][k][i][j]])
                    data['plots']['observables'][i][j].append(
                        [float('nan'), 0])

            if data['model.data.tExp'][0]:
                c = len(data['model.data.tFine'][0][0])
                for k in range(len(data['model.data.tExp'])):
                    for j in range(len(data['model.data.tExp'][k][0])):
                        for l in range(len(data['model.data.tExp'])):
                            data['plots']['observables'][i][c].append(
                                [float('nan'), float('nan')])
                            if l is k:
                                data['plots']['observables'][i][c].append(
                                    [data['model.data.yExp'][l][i][j],
                                        data['model.data.yExpStd'][l][i][j]])
                            else:
                                data['plots']['observables'][i][c].append(
                                    [float('nan'), float('nan')])
                        c = c + 1

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
        data['plots']['variables'].append(
            data['model.condition.tFine'][0][0].copy())
        for j in range(len(data['model.condition.tFine'][0][0])):
            data['plots']['variables'][i][j] =\
                [data['plots']['variables'][i][j]]
            for k in range(len(data['model.condition.xFineSimu'])):
                data['plots']['variables'][i][j].append(
                    data['model.condition.xFineSimu'][k][i][j])

    data['plots']['inputs'] = []

    # if len(data['model.u']) > 0:
    #     if isinstance(data['model.u'][0], str):
    #         data['model.condition.uFineSimu'] =\
    #             [data['model.condition.uFineSimu']]

    for i in range(len(data['model.u'])):
        data['plots']['inputs'].append(
            data['model.condition.tFine'][0][0].copy())
        for j in range(len(data['model.condition.tFine'][0][0])):
            data['plots']['inputs'][i][j] = [data['plots']['inputs'][i][j]]
            for k in range(len(data['model.condition.uFineSimu'])):
                data['plots']['inputs'][i][j].append(
                    data['model.condition.uFineSimu'][k][i][j])

    # remove unecessary data to decrease the traffic
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
    data.pop('model.z', None)

    # remove spaces (spaces wont work in css/html ids)
    try:
        for i, key in enumerate(data['model.xNames']):
            data['model.xNames'][i] = key.replace(
                ' ', '_').replace('(', '').replace(')', '')
    except:
        pass
    try:
        for i, key in enumerate(data['model.u']):
            data['model.u'][i] = key.replace(
                ' ', '_').replace('(', '').replace(')', '')
    except:
        pass
    try:
        for i, key in enumerate(data['model.data.yNames'][0]):
            data['model.data.yNames'][0][i] = key.replace(
                ' ', '_').replace('(', '').replace(')', '')
    except:
        pass

    return data


def editor(path, option, filename, content=None):
    """Reads and writes the model files for the Editor tab.

       Checks if the file is indeed in the working directory and
       denies access outside of it.
    """

    newpath = os.path.join(os.getcwd(), os.path.dirname(filename))

    if newpath.startswith(path):
        file_content = {}

        print("jojojo")

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
    else:
        print("Illegal filename/directory.")
        return 0


def create_filetree(path=None, depth=0, max_depth=0):

    tree = None

    if max_depth == 0 or depth < max_depth:
        if path is None:
            path = os.getcwd()

        tree = dict(name=os.path.basename(path), children=[])

        try:
            lst = sorted(os.listdir(path))
        except OSError:
            pass  # ignore errors
        else:
            for name in lst:
                fn = os.path.join(path, name)
                if (
                    os.path.isdir(fn) and re.match('^.*(Compiled|Results)$',
                                                   fn)
                    is None
                ):
                    child = create_filetree(fn, depth + 1, max_depth)
                    if child is not None:
                        tree['children'].append(child)
                elif re.match('^.*\.(m|def|txt|csv)$', fn) is not None:
                    tree['children'].append(dict(name=fn.replace(
                        os.getcwd() + os.path.sep, "")))

    return tree


def d2d_close_instance(uid):
    """Shut's down the d2d instance and thread."""
    while True:
        try:
            if d2d_instances[uid]['alive'] == 0:
                path = os.path.dirname(d2d_instances[uid]['model'])

                # cleans up the global dict, deletes the instance
                del(d2d_instances[uid])

                # remove the copied temporary files
                if (COPY is True) and (str(uid) in path):
                    # set file acess to writeable before deletion
                    # (maybe necessary for Windows only)
                    for root, dirs, files in os.walk(path):
                        for fname in files:
                            full_path = os.path.join(root, fname)
                            os.chmod(full_path, stat.S_IWRITE)
                    remove_tree(path)
                break

            d2d_instances[uid]['alive'] = 0
            time.sleep(SESSION_LIFETIME)
        except Exception as e:
            print('Unable to shutdown thread ' + str(uid))
            print(e)
            break


class JSONEncoder_ignore(JSONEncoder):
    """Set the encoder to ignore nan values, which will transfere all 'NaN',
    infinity etc. to "null" - somehing dygraph understands."""
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
            host=HOST,
            port=int(PORT)
            )
