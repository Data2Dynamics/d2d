from math import log10, floor
import math

fields = {'pLabel', 'p', 'lb', 'ub', 'chi2fit', 'config.optim.MaxIter','config.nFinePoints',
'fit.iter', 'fit.iter_count', 'fit.improve', 'fit.chi2_hist', 'fit.p_hist',
'model.u', 'model.pu', 'model.py', 'model.pystd', 'model.pv', 'model.pcond', 'model.fu',
'model.name', 'model.xNames','model.z', 'model.description',
'model.condition.tFine', 'model.condition.uFineSimu', 'model.data.yNames',
'model.data.tFine', 'model.data.tExp', 'model.plot.dLink', 'model.data.yFineSimu',
'model.data.ystdFineSimu', 'model.data.yExp', 'model.data.yExpStd',
'model.condition.xFineSimu', 'model.condition.zFineSimu'}
fields2 = {'model.data.tFine', 'model.data.tExp', 'model.plot.dLink', 'model.data.yFineSimu', 'model.data.ystdFineSimu',
             'model.data.yExp', 'model.data.yExpStd', 'model.condition.xFineSimu', 'model.condition.zFineSimu'}

def loopit(data):
    """Prepares the selected data for use in dygraphs.js.
        Converts float('nan') to dygraph compatible values:
        "'NaN'" creates an actual gap in the data, "None" still allows to
        connect the seperated points.
        Errors need to be set to 0 if not available.
        Sometimes uFineSimu contains float('inf') (javascript: infinity)
        which does not work in dygraphs - set to None instead.

        Finally combibnes all plot data into a dygraphs compatible data set.
        """

    for j in range(len(data['data']['yFineSimu'])):

        data['data']['yFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                        data['data']['yFineSimu'][j]]
        data['data']['ystdFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                           data['data']['ystdFineSimu'][j]]
        data['data']['yExp'][j] = [[convert_value(i) for i in x] for x in data['data']['yExp'][j]]
        data['data']['yExpStd'][j] = [[convert_value(i) for i in x] for x in data['data']['yExpStd'][j]]

        data['condition']['uFineSimu'] = [[convert_value(i) for i in x] for x in
                                          data['condition']['uFineSimu']]

    return data

def loopit2(data):
    """Prepares the selected data for use in dygraphs.js.
        Converts float('nan') to dygraph compatible values:
        "'NaN'" creates an actual gap in the data, "None" still allows to
        connect the seperated points.
        Errors need to be set to 0 if not available.
        Sometimes uFineSimu contains float('inf') (javascript: infinity)
        which does not work in dygraphs - set to None instead.

        Finally combibnes all plot data into a dygraphs compatible data set.
        """

    for j in range(len(data['model.data.yFineSimu'])):

        data['model.data.yFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                        data['model.data.yFineSimu'][j]]
        data['model.data.ystdFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                           data['model.data.ystdFineSimu'][j]]
        data['model.data.yExp'][j] = [[convert_value(i) for i in x] for x in data['model.data.yExp'][j]]
        data['model.data.yExpStd'][j] = [[convert_value(i) for i in x] for x in data['model.data.yExpStd'][j]]

        data['model.condition.uFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                          data['model.condition.uFineSimu'][j]]
        data['model.condition.xFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                          data['model.condition.xFineSimu'][j]]
        data['model.condition.zFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                          data['model.condition.zFineSimu'][j]]


    return data

def loopit3(data):  # 10ms

    for j in range(len(data['data']['yFineSimu'])):

        data['data']['yFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                        data['data']['yFineSimu'][j]]
        data['data']['ystdFineSimu'][j] = [[convert_value(i) for i in x] for x in
                                           data['data']['ystdFineSimu'][j]]
        data['data']['yExp'][j] = [[convert_value(i) for i in x]
                                   for x in data['data']['yExp'][j]]
        data['data']['yExpStd'][j] = [[convert_value(i) for i in x]
                                      for x in data['data']['yExpStd'][j]]

        data['condition']['uFineSimu'] = [[convert_value(i) for i in x] for x in
                                          data['condition']['uFineSimu']]
        data['condition']['zFineSimu'] = [[convert_value(i) for i in x] for x in
                                          data['condition']['zFineSimu']]

    return data


def convert_value(x):

    if (x != x) or (math.isinf(x)):
        return None
    else:
        return x

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def create_dygraphs_data2(data):
    """Prepares the selected data for use in dygraphs.js.
        Converts float('nan') to dygraph compatible values:
        "'NaN'" creates an actual gap in the data, "None" still allows to
        connect the seperated points.
        Errors need to be set to 0 if not available.
        Sometimes uFineSimu contains float('inf') (javascript: infinity)
        which does not work in dygraphs - set to None instead.

        Finally combibnes all plot data into a dygraphs compatible data set.
        """

    for j in range(len(data['data']['yFineSimu'])):

        data['data']['yFineSimu'][j] = [[None if i != i else
                                        i for i in x] for x in
                                        data['data']['yFineSimu'][j]]
        data['data']['ystdFineSimu'][j] = [[0 if i != i else
                                           i for i in x] for x in
                                           data['data']['ystdFineSimu'][j]]
        data['data']['yExp'][j] = [[None if i != i else i for i in x]
                                   for x in data['data']['yExp'][j]]
        data['data']['yExpStd'][j] = [[0 if i != i else i for i in x]
                                      for x in data['data']['yExpStd'][j]]

        data['condition']['uFineSimu'] = [[None if i == float('inf') else
                                          i for i in x] for x in
                                          data['condition']['uFineSimu']]

    try:
        data['fit']['p_hist'] = [['NaN' if i != i else i for i in x]
                                 for x in data['fit']['p_hist']]
        data['fit']['chi2_hist'] = [['NaN' if i != i else i for i in x]
                                    for x in data['fit']['chi2_hist']]
    except:
        pass

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

    data['numpy'] = {}
    data['numpy']['test'] = numpy.asarray(data['condition']['tFine'])


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
