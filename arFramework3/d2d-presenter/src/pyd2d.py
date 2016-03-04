import io
import os
import matlab.engine
import numpy
from math import log10, floor

__author__ = 'Clemens Blank'


class mat():
    def __init__(self):
        self.output_total = ''
        self.eng = matlab.engine.start_matlab()
        self.cd(os.getcwd())

    def cd(self, arg):
        self.eval("cd '" + str(arg) + "'", 0)

    def eval(self, exp, args=0):
        print(str(exp) + "  - " + str(args))
        out = io.StringIO()
        err = io.StringIO()
        try:
            result = self.eng.eval(exp, stdout=out, stderr=err,
                                   nargout=args)
            self.output_total += out.getvalue()
            print(out.getvalue(), end='')
        except:
            result = False
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

    def round_sig(self, x, sig=2):
        """ round x to sig significant digits"""
        return round(x, sig-int(floor(log10(abs(x))))-1)

    def convert(self, x, key=None, convert='numpy'):
        """ removes the matlab double type"""
        if convert is 'list':
            return numpy.reshape(numpy.asarray(x._data),
                                              (x.size[1], x.size[0])).tolist()
        else:
            return numpy.reshape(numpy.asarray(x._data),
                                              (x.size[1], x.size[0]))

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

    def load_model(self, path, load=False):
        self.path = os.path.dirname(path)
        self.filename = os.path.basename(path)
        self.cd(self.path)

        if load != False:
            self.eval("arLoad('" + str(load) + "')")
        else:
            self.eval(os.path.splitext(os.path.basename(path))[0])
        self.update()

    def fastload(self):
        self.load_model(os.path.join(os.getcwd(), 'models/Raia_CancerResearch2011/Setup.m'))

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

            self.eval("arSetPars(" + str(
                dict.keys()).replace('dict_keys([', '{').replace(
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
            data = dict(zip(fields, self.eng.eval("arGet(" +
                        str(fields).replace('[', '{').replace(']', '}') +
                ", " + str(m) + ", " + str(dset) + ", " + str(sig) + ");", nargout=1)))
        except:
            print("Get failed. Please make sure all fields exist.")
            return None

        # Converts matlab.double arrays to python lists
        if convert is not False:
            for key in data:
                if type(data[key]) == matlab.double:
                    data[key] = self.convert(data[key], key, convert)
                elif type(data[key]) == list:
                    for i in range(len(data[key])):
                        if type(data[key][i]) == matlab.double:
                            data[key][i] = self.convert(
                                                        data[key][i],
                                                        key,
                                                        convert)

        return data

    def simu(self, sensi='true', fine='true', dynamics='true'):
        """Executes arSimu.

           Tries to reload save parameters when arSimu fails (e.g.
           CV_CONV_FAILURE) when using a bad parameter setting.
        """
        status = self.eval(
            'arSimu(' + sensi + ', ' + fine + ', ' + dynamics + ');')

        if status is False:

            params = self.get({'p', 'd2d_presenter.p', 'config.savepath'})

            if len(params['p']) == len(params['d2d_presenter.p']):
                self.eval("ar.p = ar.d2d_presenter.p;")
                status = 'p_reset'
                return status
            elif params['config.savepath']:
                self.eval("arLoad('" + params['config.savepath'] + "')")
                status = 'arLoad'
                return status
            else:
                status = 'broken'
                return status
        else:
            return status
