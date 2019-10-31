# -*- coding: utf-8 -*-
"""
Convert Vicon csv to MATLAB mat

Note that I made a python code because I find it a lot easier to code this here

@author: Luke Sy z5151460
"""

import numpy as np
import pandas as pd
import scipy.io as sio
import os
import argparse
from StringIO import StringIO

def parseViconCSV(fname, output='viconcsv.mat'):
    retval = {}
    with open(fname, 'r') as f:
        buf0 = []
        
        for line in f:
            if line == '\r\n':
                df = pd.read_table(StringIO('\n'.join(buf0[2:])), sep=',', lineterminator='\n')
                df = df.iloc[:, :-1] # remove last column
                
                # fix column names
                newcolumnname = ['']
                for i in df.columns:
                    n = i.split(':')
                    if n[0] == 'Unnamed':
                        newcolumnname.append(newcolumnname[-1])
                    else:
                        n2 = n[1].split('.')
                        newcolumnname.append(n2[0])                
                
                for i, x in enumerate(df.loc[0,:]):
                    if x in ["X''", "Y''", "Z''"]:
                        newcolumnname[i+1] += "A"
                    elif x in ["X'", "Y'", "Z'"]:
                        newcolumnname[i+1] += "V"
                        
                df.columns = newcolumnname[1:]
                
                # get unique columns
                for i in set(df.columns):
                    if i != '':
                        retval[i] = df.loc[2:, i].values.astype(float)
                        
                # print(df.head())
                buf0 = []
            else:
                buf0.append(line.strip())
                
        sio.savemat(output, retval)

parser = argparse.ArgumentParser(description='parse exported vicon .csv to .mat')
parser.add_argument('-f', dest='dir', action='store', required=True,
                    help='html file name', default='.')
args = parser.parse_args()

file_list = os.listdir(args.dir)
csvs = filter(lambda x: x[-4:] == '.csv', file_list)
mats = filter(lambda x: x[-4:] == '.mat', file_list)

for i in csvs:
    if i[:-4]+'.mat' not in mats:
        n = os.path.join(args.dir, i)
        print(n)
        parseViconCSV(n, output='{}.mat'.format(n[:-4]))
