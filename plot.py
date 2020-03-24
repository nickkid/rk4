# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 19:31:21 2019

@author: Lihuan Yuan
"""
import numpy as np
import matplotlib.pyplot as plt
import csv

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
min_record_difference = 1e-3

plt.title('Time and Absorbance')
plt.xlabel('Time')
plt.ylabel('Absorbance')
plt.xscale('log')
plt.legend(('[NO]=0.47uM', '[NO]=1.53uM', '[NO]=6.6uM', '[NO]=500uM'), loc='upper right')
m = ('o','s','*','+')
Item_name = True
First_row = True

for i in range(0,Number_Of_Cases):
    filename = 'results{0}.csv'.format(i)
    with open(filename) as csvfile:
        plots = csv.reader(csvfile, delimiter=',')
        raw_data = np.empty((0,2))
        for row in plots:
            if Item_name:
                Item_name = False
                continue
            if First_row:
                First_row = False
                t0 = float(row[0])
                A0 = float(row[5])
                continue
            t = float(row[0])
            A = float(row[5])
            temp = np.array([t, A])
            raw_data = np.vstack((raw_data, temp))
        
        x = raw_data[:,0]
        A = raw_data[:,1]
        plt.plot(x, A, marker=m[i])
        
plt.show()
            

