#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 13:50:26 2016

@author: aholaj
"""

###### read ice

from netCDF4 import Dataset
import numpy as np
import sys
import os

filename = sys.argv[1]
var=sys.argv[2]
fileH = Dataset(filename,mode='r')

data = fileH.variables[var][:]

fileH.close()
head, tail = os.path.split(filename)
print ' '
print 'file: ' + tail
print 'variable: ' + var
print 'min: '+ str(np.min(data)) + ' max: ' + str(np.max(data))