# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 19:31:21 2019

@author: Lihuan Yuan
"""
import numpy as np
import pandas

Number_Of_Cases = 1
'''
epsilon of sGC, 6C-sGC-NO and 5C-sGC-NO in M^-1cm^-1
6C: 6C sGC-NO
5C: 5C sGC-NO

Absorbance A = epislon * concentration * l
l: length in cm
See https://en.wikipedia.org/wiki/Molar_attenuation_coefficient

epsilon_sGC = 124e-3
epsilon_6C = 66e-3
epsilon_5C = 46e-3
'''
epsilons = np.array([124e3, 66e3, 46e3])
length = 100# cm


for i in range(0,Number_Of_Cases):
    filename = 'results{0}.csv'.format(i)
    pd = pandas.read_csv(filename)
    print(pd[1,:])