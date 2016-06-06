from netCDF4 import Dataset, num2date
import numpy as np
from datetime import datetime
import glob
from StringIO import StringIO
import dateutil.parser as dparser
import file_saver

"""
    parseRadiometrics.py
    Author: Greg Blumberg (OU/CIMMS)
    Email: wblumberg@ou.edu

    This script will parse the radiometrics text file of the Level 1 data (calibrated Tb variables)
    and transform the data into a daily netCDF file that can be read by MWRoe.  This data
    contains all available brightness temperature data (even off-zenith data) as well as surface
    pressure, surface temperature, and surface RH data from the MET tower provided with the radiometer.

    The lat/lon/alt of the radiometer is hard coded and should be changed by the person running this code
    to ensure that the metadata of the file is correct.

    This code was originally constructed to transform the University of Manitoba MWR data from PECAN
    into data that could be read by MWRoe.

"""

def parseRecords(lines):
    # Parse out the header
    header = lines[2].strip().replace("Ch","").split(",")
    header = header[:len(header)-1]

    # Parse out the frequency values
    freq = np.asarray(header[6:], dtype=float)
    data = lines[4:]
   
    # Parse out the lines that correspond to the surface met data 
    atmos_data = StringIO('\n'.join(data[::4]))
    atmos_data = np.genfromtxt(atmos_data, dtype='object', delimiter=',', unpack=True, converters={1: dparser.parse})
    
    # Pull out the brightness temperature observations
    all_idx = np.arange(0, len(data), 1)
    wrong_idx = np.arange(0, len(data), 4)
    bt_idx = np.setxor1d(all_idx, wrong_idx)
    bt_data = StringIO('\n'.join(np.asarray(data)[bt_idx]))
    bt_data = np.genfromtxt(bt_data, dtype='object', delimiter=',', unpack=True, converters={1: dparser.parse}).T
    bts = []
   
    # Every 3rd line is a new "observation"
    # Parse out each line.
    for i in np.arange(0, len(bt_data),3):
        sample = bt_data[i:i+3, 6:]
        if len(sample) != 3: # If this sample is incomplete, just ditch the data and break this loop
            break
        sample = np.where(sample != '', sample, -9999)
        sample = np.ma.masked_where(sample == -9999, sample)
        sample = np.ma.asarray(sample, dtype=float)
        z_bt = sample[1]
        elevs = np.asarray(bt_data[np.arange(i, i+3, 1),4], dtype=float)
        if not np.isclose(elevs[2], 180-elevs[0], rtol=1e-1):
            continue
        # Average the off-zenith values
        oz_bt = np.ma.mean(sample[[0,2],:], axis=0)
        z_bt = z_bt[:len(z_bt)-1]
        oz_bt = oz_bt[:len(oz_bt)-1]
        bts.append(np.ma.vstack((z_bt, oz_bt)))
    atmos_data = atmos_data.T
    bts = np.ma.asarray(bts)
    elev = np.unique(np.asarray(bt_data[:,4], dtype=float))
    elev = elev[:len(elev)/2 + 1]
    dts = atmos_data[:,1] 
    rainflag = np.asarray(atmos_data[:bts.shape[0],-1], dtype=int)
    #irTemp = np.asarray(atmos_data[:len(bts),-2], dtype=float)
    pres = np.asarray(atmos_data[:bts.shape[0],-3], dtype=float)
    temp = np.asarray(atmos_data[:bts.shape[0],-5], dtype=float)
    rh = np.asarray(atmos_data[:bts.shape[0],-4], dtype=float)
    # Variables that I have now:
    # freq, bts, pres, temp, rainflag, rh, elev
    dd = {}
    dd['elevation_angle'] = elev[::-1]
    dd['frequencies'] = freq
    dd['brightness_temperature'] = bts
    dd['air_temperature'] = temp
    dd['air_pressure'] = pres
    dd['flag'] = rainflag
    dd['dts'] = dts
    dd['lon'] = -99.5739974975586
    dd['lat'] = 38.95800018310547
    dd['alt'] = 237.4
    dd['radiometer_name'] = 'Ellis FP Uni. Manitoba MWR'
    dd['radiometer_system'] = 'Radiometrics'
    return dd

files = np.sort(glob.glob("raw/0*/"))
for f in files:
    files_for_date = np.sort(glob.glob(f + '*lv1.csv'))
    dd = None
    for fi in files_for_date:
        print fi
        data = open(fi, 'r')
        try:
            di = parseRecords(data.readlines()) 
        except Exception,e:
            print e
        if dd == None:
            dd = di
        else:
            for key in ['dts','air_pressure', 'air_temperature', 'flag', 'brightness_temperature']:
                print key
                dd[key] = np.ma.concatenate((dd[key],di[key]))
    file_saver.save(dd, 'ellis')
print files


