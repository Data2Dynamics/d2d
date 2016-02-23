import io
import os
import re
import matlab.engine
import numpy
from math import log10, floor
import math
import copy
import test

__author__ = 'Clemens Blank'

def convert3(x, ROUND=0):

    conv = numpy.asarray([])

    for _ in range(x.size[0]):
        lst = numpy.asarray(x._data[_::x.size[0]].tolist())


        if ROUND is not 0:
            lst = numpy.around(lst, ROUND)

        conv = numpy.append(conv, lst)
    return conv

class mat():
    def __init__(self):
        self.eng = matlab.engine.start_matlab()
        self.cd(os.getcwd())
        self.output_total = ''
        self.convert_value_vec = numpy.vectorize(self.convert_value, otypes=[numpy.object])
        self.round_sig_vec = numpy.vectorize(self.round_sig)

    def cd(self, arg):
        self.eng.cd(arg, nargout=0)

    def eval(self, exp, args=0):

        out = io.StringIO()
        err = io.StringIO()
        try:
            result = self.eng.eval(exp, stdout=out, stderr=err,
                                   nargout=args)
            self.output_total += out.getvalue()
            print(out.getvalue(), end='')
        except:
            result = None
            self.output_total += err.getvalue()
            print(err.getvalue(), end='')

        out.close()
        err.close()

        return result

    def set(self, dict, silent=True):
        if silent is True:
            silent = ';'
        else:
            silent = ','

        command = ""

        for key in dict:
            if type(dict[key]) == str:
                dict[key] = "'" + dict[key] + "'"

            command += key + "=" + str(dict[key]) + silent

        self.eval(command)

    def convert(self, x, ROUND=0):

        conv = []

        for _ in range(x.size[0]):
            lst = x._data[_::x.size[0]].tolist()
            conv.append(lst)

        return conv

    def round_sig(self, x, sig=2):

        return round(x, sig-int(floor(log10(abs(x))))-1)

    def convert_value(self, x):

        if (x != x) or (math.isinf(x)):
            return None
        else:
            return x

    def convert2(self, x, key=None):
        """ removes the matlab double type"""
        #return self.convert(x)
        # return [x._data[i:i+x.size[0]] for i in range(0,len(x._data), x._size[0])] # 33ms mit arrays
        #return list(zip(*[iter(x._data.tolist())]*x.size[0]))  # (35ms) mit tuplen
        #return [x._data[i:i+x.size[0]] for i in range(0,len(x._data.tolist()), x._size[0])] # 36ms mit arrays
        #return [x._data.tolist()[x.size[0]*i: x.size[0]*(i+1)] for i in range(x.size[1])] # --> way to slow (120ms to 32ms)

        return numpy.reshape(numpy.asarray(x._data), (x.size[1], x.size[0])).tolist() # (75ms)

        #num = numpy.reshape(numpy.asarray(x._data), (x.size[1], x.size[0])).tolist()

        # to speed it up, only the areas with graph data gets converted
        #if key.startswith('model.condition.') or key.startswith('model.data.'):
        #    num = self.convert_value_vec(num).tolist()
        #return num  # (40ms)

    def convert_struct2(self, dic):
        """ Converts all matlab.double types in structures into python lists.
        """

        try:
            dic = dic.copy()

            if isinstance(dic, dict):
                for key in dic:
                    dic[key] = self.convert_struct2(dic[key])
            elif isinstance(dic, list):
                try:
                    dic[0][0]  # check if already converted matlab.double
                    if isinstance(dic[0], matlab.double):
                        for i in range(len(dic)):
                            dic[i] = self.convert_struct2(dic[i])
                except:
                    for i in range(len(dic)):
                        dic[i] = self.convert_struct2(dic[i])

        except:
            if isinstance(dic, (str, float, bool, int, matlab.logical)):
                pass
            elif isinstance(dic, matlab.double):
                dic = self.convert2(dic)

        return dic

    def convert_struct(self, dic, ROUND=0):
        """ Converts all matlab.double types in structures into python lists.
        """

        try:
            dic = dic.copy()

            if isinstance(dic, dict):
                for key in dic:
                    dic[key] = self.convert_struct(dic[key], ROUND)
            elif isinstance(dic, list):
                try:
                    dic[0][0]  # check if already converted matlab.double
                    if isinstance(dic[0], matlab.double):
                        for i in range(len(dic)):
                            dic[i] = self.convert_struct(dic[i], ROUND)
                except:
                    for i in range(len(dic)):
                        dic[i] = self.convert_struct(dic[i], ROUND)

        except:
            if isinstance(dic, (str, float, bool, int, matlab.logical)):
                pass
            elif isinstance(dic, matlab.double):
                dic = self.convert(dic, ROUND)

        return dic



    def benchmark(self, method, args, runs=200):
        """ Can be used to measure the execution time of methods """

        import timeit
        a = timeit.default_timer()

        for i in range(runs):
            method(*args)

        b = timeit.default_timer() - a

        print(str(method) + " took an average time of " + str(b / runs))
        print(str(method) + " took for " + str(runs) + " runs " + str(b))
        return b

    def stop(self):
        self.__del__()

    def __del__(self):
        self.eng.quit()
        print("Beendet.")


class d2d(mat):

    """Enhances the mat class with d2d specific methods."""

    def load_model(self, path, load=None):
        self.path = os.path.dirname(path)
        self.filename = os.path.basename(path)
        self.cd(self.path)

        if load != None:
            self.eval("arLoad(" + str(load) + ")")
        else:
            self.eval(os.path.splitext(os.path.basename(path))[0])

        self.eval("arLink, arSimu(true, true, true)")
        self.ar = self.eval('arToPython;', 1)[0]

    def fastload(self):
        self.load_model(os.path.join(os.getcwd(), 'models/Raia_CancerResearch2011/Setup.m'))
        import d2d_presenter
        self.data = d2d_presenter.select_data(self, '')
        self.datac2 = self.convert_struct2(self.data)
        self.datac = self.convert_struct(self.data)
        print("datac2")
        print(self.datac2['data']['yExp'])
        print("datac")
        print(self.datac['data']['yExp'])

    def update(self):
        self.ar = self.eval("arToPython;", 1)[0]

    def dic2obj(self, d=None, convert=None):
        """Converts the Python ar structure to a more comfortable object. E.g.
            "ar['model'][0]['data'][1]['yExp']" becomes
            "ar.model[0].data[1].yExp".
            Provides the ability to convert all matlab double fields to
            Python lists, but keep in mind that this might increase
            calculation time by faktor 15(!) on common ar structs.
            It is reccomended to select the fields of interest first and
            convert only these via convert or convert_struct."""

        if d is None:
            d = self.ar

        top = type('new', (object,), d)
        seqs = tuple, list, set, frozenset
        for i, j in d.items():
            if convert is not None:
                if type(j) is matlab.double:
                    j = self.convert(j)
            if isinstance(j, dict):
                setattr(top, i, self.dic2obj(j, convert))
            elif isinstance(j, seqs):
                setattr(top, i, type(j)(self.dic2obj(sj, convert) if
                isinstance(sj, dict) else sj for sj in j))
            else:
                setattr(top, i, j)
        return top

    def set_pars_from_dict(self, dict):

        try:
            dict = {k: dict[k] for k in self.ar['pLabel'][0] if k in dict}

            self.eval("arSetPars(" +
                  str(dict.keys()).replace('dict_keys([', '{').replace(
                      '])', '}') + ', ' +
                  str(dict.values()).replace('dict_values(', '').replace(
                      ')', '').replace("'", "") + ");")
        except:
            print("Set parameters failed.")

    def get(self, fields, m=1, dset=1, sig=False, convert=False):
        """ Get's a set of fields in the ar struct with the possibility to
        chop the values after sig digits. """
        if sig is False:
            sig = 'false'

        try:
            data = dict(zip(fields, self.eng.eval("arGet(" + str(fields) +
                ", " + str(m) + ", " + str(dset) + ", " + str(sig) + ");", nargout=1)))
        except:
            print("Get failed. Please make sure all fields exist.")
            return None

        # Converts matlab.double arrays to python lists
        if convert is True:
            for key in data:
                if type(data[key]) == matlab.double:
                    data[key] = self.convert2(data[key], key)
                elif type(data[key]) == list:
                    for i in range(len(data[key])):
                        if type(data[key][i]) == matlab.double:
                            data[key][i] = self.convert2(data[key][i], key)

        return data

    def simu(self, sensi='true', fine='true', dynamics='true'):
        self.eval('arSimu(' + sensi + ', ' + fine + ', ' + dynamics + ')')
