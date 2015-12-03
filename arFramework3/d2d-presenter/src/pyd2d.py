import io
import os
import re
import matlab.engine
from math import log10, floor

__author__ = 'Clemens Blank'

class mat():
    def __init__(self):
        self.eng = matlab.engine.start_matlab()
        self.cd(os.getcwd())
        self.output_total = ''

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

            if ROUND is not 0:
                lst = [self.round_sig(elem, ROUND) if elem != 0 and
                       elem == elem else elem for elem in lst]

            conv.append(lst)

        return conv

    def round_sig(self, x, sig=2):
        return round(x, sig-int(floor(log10(abs(x))))-1)

    # testing the methods speed
    def benchmark(self, method, args, runs=200):
        import timeit
        a = timeit.default_timer()

        for i in range(runs):
            x = method(*args)

        b = timeit.default_timer() - a

        print(str(method) + " braucht durchschnittlich " + str(b / runs))
        print(str(method) + " brauchte insgesamt " + str(b))
        return b

    def stop(self):
        self.__del__()

    def __del__(self):
        self.eng.quit()
        print("Beendet.")


class d2d(mat):

    def load_model(self, path):
        self.path = os.path.dirname(path)
        self.filename = os.path.basename(path)
        self.cd(self.path)
        self.eval(os.path.splitext(os.path.basename(path))[0])
        self.ar = self.eval('arToPython;', 1)[0]

    def update(self):
        self.ar = self.eval("arToPython;", 1)[0]

    def dic2obj(self, d=None):
        if d is None:
            d = self.ar

        top = type('new', (object,), d)
        seqs = tuple, list, set, frozenset
        for i, j in d.items():
            if isinstance(j, dict):
                setattr(top, i, self.dic2obj(j))
            elif isinstance(j, seqs):
                setattr(top, i, type(j)(self.dic2obj(sj) if isinstance(
                    sj, dict) else sj for sj in j))
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

    def convert_dic(self, dic, ROUND=0):

        dic = dic.copy()

        for key in dic:
            if type(dic[key]) is dict:
                dic[key] = self.convert_dic(dic[key], ROUND)
            elif type(dic[key]) is matlab.double:
                dic[key] = self.convert(dic[key], ROUND)

        return dic

    def get(self, fields, m=1, d=1, c=1, convert=False):

        try:
            data = dict(zip(fields, self.eng.eval("arGet(" +
                        str(fields).replace('[', '{').replace(']', '}') +
                ", " + str(m) + ", " + str(d) +
                ", " + str(c) + ");", nargout=1)))
        except:
            print("get failed. Please make sure all fields exist.")
            return None

        # Converts matlab.double arrays to python lists
        if convert is True:
            for key in data:
                if type(data[key]) == matlab.double:
                    data[key] = self.convert(data[key])

        return data

    def simu(self, sensi='true', fine='true', dynamics='true'):
        self.eval('arSimu(' + sensi + ', ' + fine + ', ' + dynamics + ')')
