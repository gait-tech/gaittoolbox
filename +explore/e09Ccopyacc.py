import os
import csv
import shutil

dest = '/mnt/c/Users/z5151460/workspace/gaitrecon/neura-sparse01'
dest_imu = os.path.join(dest, 'imu/')
dest_calib = os.path.join(dest, 'calib/')

convert = {'00B40C49': '00B40C44', '00B40C4A': '00B40C47'}
# 00B40C49 to 00B40C44
# 00B40C4A to 00B40C47
for s in ['S01', 'S02']:
    # calibration
    for u, v in convert.items():
        srcfname = os.path.join(dest_calib, '{0}-Calib-SensorW2V-000_{1}.txt'.format(s, u))
        dstfname = os.path.join(dest_calib, '{0}-Calib-SensorW2V-000_{1}.txt'.format(s, v))
        print(srcfname+" "+dstfname)
        shutil.copy(srcfname, dstfname)
    
    # imu
    fnames = os.listdir(dest_imu)
    fnames = filter(lambda x: x[-4:] == '.txt' and x[0:3] == s, fnames)
    fnames = set(map(lambda x: x.split("-000_")[0], fnames))
    
    for f in fnames:
        for u, v in convert.items():
            srcfname = os.path.join(dest_imu, '{0}-000_{1}.txt'.format(f, u))
            dstfname = os.path.join(dest_imu, '{0}-000_{1}.txt'.format(f, v))
            print(srcfname+" "+dstfname)
            shutil.copy(srcfname, dstfname)