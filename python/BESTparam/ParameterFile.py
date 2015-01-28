#/usr/bin/env python
"""Generate a lookup table in C (.c and .h file) for use with
   the BEST parametrization."""

__author__     = "Nils Wenzler"
__copyright__  = "Nils Wenzler"
__maintainer__ = "Nils Wenzler"
__email__      = "wenzlern@ethz.ch"

# Imports
import numpy as np
from textwrap import fill

# Classes
class Parametrization:
    """ Manages the data"""
    def __init__(self, name, comment = '', timestep=0.01):
        self.Values = []
        self.Values_c = []
        self.Name = name
        self.Comment = comment
        self.timestep = timestep

    def AddByList(self, list, derivlist):
        self.Values = list
        self.Values_c = derivlist

    def AddValues(self, list):
        self.Values = list

    def AddValues_c(self, list):
        self.Values_c = list

    def AddByFunction(self, func, deriv, start, stop, step):
        for i in np.linspace(start, stop, (1/step) +1):
            self.Values.append(func(i))
            self.Values_c.append(deriv(i))

    def AddFunc(self, func, start, stop, step):
        for i in np.linspace(start, stop, (1/step) +1):
            self.Values.append(func(i))

    def AddFunc_c(self, func, start, stop, step):
        for i in np.linspace(start, stop, (1/step) +1):
            self.Values_c.append(func(i))

    def ToString(self):
        Header = "/*This file is auto generated and contains the LUT  and the function for\n"
        Header +=" the BEST parametrization\n*/"
        ret = Header + '\n/*' + self.Comment + '*/\n'
        ret += '\n#include \"fpCommon.h\"\n'

        values = str()
        values_c = str()
        print len(self.Values), len(self.Values_c)
        for i in range(len(self.Values)):
            values += str(self.Values[i]) + ', '
            values_c += str(self.Values_c[i]) + ', '

        # Fill wraps the text to 70 characters, [:-1] leaves away the last ','
        ret += '\nconst static double data[] = {\n'
        ret += fill(values)[:-1] + '};\n'

        ret += '\nconst static double data_c[] = {\n'
        ret += fill(values_c)[:-1]+ '};\n'

        ret += '\nDLL_EXPORT double ' + self.Name + '(short *err, double *param){\n'
        ret += '  int index = (int)(param[BATTERY_SOC]*' + str(len(self.Values)-1) + ');\n'
        ret += '  float rest  = param[BATTERY_SOC]*' + str(len(self.Values)-1) + '- index;\n'
        ret += '  return data[index] + rest*(data[index+1]-data[index]);\n}\n'


        ret += '\nDLL_EXPORT double ' + self.Name + '_c(short *err, double *param){\n'
        ret += '  int index = (int)(param[BATTERY_SOC]*' + str(len(self.Values_c)-1) + ');\n'
        ret += '  float rest  = param[BATTERY_SOC]*' + str(len(self.Values_c)-1) + '- index;\n'
        ret += '  return data_c[index] + rest*(data_c[index+1]-data_c[index]);\n}\n'

        return ret

    def SaveToFile(self,Filename):
        f = open(Filename, 'w')
        f.write(self.ToString())
        f.close()
