from datetime import datetime
from netCDF4 import Dataset, date2num
import numpy as np

def save(dd, pisa_name):
    fn = pisa_name + 'radiomwrC1.b1.' + datetime.strftime(dd['dts'][0], '%Y%m%d.%H%M%S.cdf')
    epoch_units = 'seconds since 1970-01-01 00:00:00+00:00'
    epoch_times = date2num(dd['dts'], epoch_units)

    d = Dataset(fn, 'w', format='NETCDF3_CLASSIC')
    time_dim = d.createDimension('time', None)
    freq_dim = d.createDimension('freq', dd['brightness_temperature'].shape[2])
    elev_dim = d.createDimension('angle', dd['elevation_angle'].shape[0])

    d.system = dd['radiometer_name'] + ' ' + dd['radiometer_system']
    d.institution = 'Data processed by Greg Blumberg at OU/CIMMS; contact: wblumberg@ou.edu'
    
    data = d.createVariable('latitude', 'f4')
    data.units = "degrees_north"
    data.long_name = "Latitude of microwave radiometer"
    data[:] = dd['lat']
    
    data = d.createVariable('longitude', 'f4')
    data.units = "degrees_east"
    data.long_name = "Longitude of microwave radiometer"
    data[:] = dd['lon']
    
    data = d.createVariable('altitude', 'f4')
    data.units = "m"
    data.long_name = "Height of microwave radiometer above mean sea level"
    data[:] = dd['alt']
    
    data = d.createVariable('time', 'i4', 'time')
    data[:] = np.arange(0, len(epoch_times))

    data = d.createVariable('time_since_19700101', 'i4', ('time',))
    data.units = epoch_units
    data.long_name = 'Time UTC'
    data[:] = epoch_times

    data = d.createVariable('frequencies', 'f4', ('freq',))
    data.units = "GHz"
    data.long_name = "Frequencies of brightness temperatures"
    data[:] = dd['frequencies']
    
    data = d.createVariable('elevation_angle', 'f4', ('angle',))
    data.units = "GHz"
    data.long_name = "Frequencies of brightness temperatures"
    data[:] = dd['elevation_angle']

    data = d.createVariable('brightness_temperature', 'f4', ('time','angle','freq',))
    data.units = "K"
    data.long_name = "brightness_temperature"
    data[:] = dd['brightness_temperature']
    
    data = d.createVariable('flag', 'i4', ('time',))
    data.long_name = 'Rain flag'
    data.flag_meanings = '0 - clear, 1 - rain'
    data[:] = dd['flag']

    data = d.createVariable('air_temperature', 'f4', ('time',))
    data.long_name = "enviromental temperature"
    data.units = 'K'
    data[:] = dd['air_temperature']

    data = d.createVariable("air_pressure", 'f4', ('time',))
    data.long_name = "environmental pressure"
    data.units = 'mb'
    data[:] = dd['air_pressure']

    d.close()

