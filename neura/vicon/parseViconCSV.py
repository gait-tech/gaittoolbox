# -*- coding: utf-8 -*-
"""
Convert Vicon csv to MATLAB mat

Note that I made a python code because I find it a lot easier to code this here

@author: Luke Sy z5151460
"""

import numpy as np
import scipy.io as sio
import os

def parseViconCSV(fname, output='viconcsv.mat'):
    with open(fname, 'r') as f:
        buf0 = f.readline().strip()
        buf1 = f.readline().strip()
        buf2 = f.readline().strip()
               
        retval = {}
        while buf2:
            if buf1[:5]=='Frame':
                var_names = filter(lambda x: x.strip() != '', buf0.split(','))
                var_names = list(map(lambda x: x.split(':')[1], var_names))
                # var_desc = filter(lambda x: x != '', buf1.split(','))
                # var_units = filter(lambda x: x != '', buf2.split(','))
                
                data = []
                
                buf2 = f.readline().strip()
                while buf2:
                    buf3 = list(map(lambda x: float(x.strip()) if x.strip()!='' else np.nan, buf2.split(',')))
                    buf3.extend([np.nan]*(3*len(var_names)+2-len(buf3)))
                    
                    data.append(buf3)
                    buf2 = f.readline().strip()
                
                data = np.array(data)
                
                print(var_names)
                print(data.shape)
                
                for i,j in enumerate(var_names):
                    retval[j] = data[:,i*3+2:i*3+5]
                buf1 = ''
            
            buf0 = buf1
            buf1 = buf2
            buf2 = f.readline().strip()
            
        sio.savemat(output, retval)

file_list = os.listdir('.')
csvs = filter(lambda x: x[-4:] == '.csv', file_list)
mats = filter(lambda x: x[-4:] == '.mat', file_list)

for i in csvs:
    if i[:-4]+'.mat' not in mats:
        print(i[:-4])
        parseViconCSV(i, output='{}.mat'.format(i[:-4]))
