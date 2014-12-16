#/usr/bin/env python
"""Generate a lookup table in C (.c and .h file) for use with
   the BEST parametrization."""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
import numpy as np
import textwrap

# Classes
class Parametrization:
    """ Manages the data"""
    def __init__(self, name, comment = ''):
        self.Values = []
        self.Values_c = []
        self.Name = name
        self.Comment = comment

    def AddList(self, list):
        self.Values = list

    def AddByFunction(self, func, deriv, start, stop, step):
        for i in np.arange(start, stop, step):
            self.Values.append(func(i))
            self.Values_c.append(deriv(i))

    def ToString(self):
        Header = "/*This file is auto generated and contains the LUT  and the function for\n"
        Header +=" the BEST parametrization\n*/"
        ret = Header + '\n/*' + self.Comment + '/*\n'
        ret += '\n#include \"fpCommon.h\"\n'

        values = str()
        values_c = str()
        for i in range(len(self.Values)):
            values += str(self.Values[i]) + ', '
            values_c += str(self.Values_c[i]) + ', '

        ret += '\nconst static data[] = {\n'
        ret += fill(values) + '}\n'

        ret += '\nconst static data_c[] = {\n'
        ret += fill(values_c) + '}\n'

        ret += 'DLL_EXPORT double ' + self.name + '(short *err, double *param){\n'
        ret += '  return data[param[BATTERY_SOC]*' + str(len(self.Values)) + '];\n}\n'


        ret += 'DLL_EXPORT double ' + self.name + '_c(short *err, double *param){\n'
        ret += '  return data_c[param[BATTERY_SOC]*' + str(len(self.Values_c)) + '];\n}'

        return ret

    def SaveToFile(self,Filename):
        f = open(Filename, 'w')
        f.write(self.ToString())
        f.close()
