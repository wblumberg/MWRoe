from netCDF4 import Dataset, date2num
import os
import numpy as np
import time as tm
from dateutil import tz
from datetime import datetime
import helper

def constructOutputFN(dts, config_dict):
    '''
        constructOutputFN

        This function constructs a string that is the output filename 
        for the file that contains the retrieval output.

        Parameters
        ----------
        dts: a list of datetime objects for each retrieval observation that
             will be retrieved.
        config_dict : a dictionary that contains the configuration variables for the
                      the retrieval.

        Returns
        -------
        a string that is the filename where the retireval output will be saved
    '''

    utc = tz.gettz('UTC')
    dt = dts[0].replace(tzinfo=utc)
    dt_string = datetime.strftime(dt, '.%Y%m%d.%H%M%S.cdf')

    return config_dict['output_path'] + "/" + config_dict['output_rootname'] + dt_string


def save_retrieval(out_filename, output, config_dict, prior_info, input_info):
    '''
        save_retrieval

        This function writes the retrieval output to the netCDF file.
        If the netCDF file already exists, this function appends the retrieval
        output to the file.

        Parameters
        ----------
        out_filename: a string containing the output filename for the retrieval output file
        output : a dictionary containing the retrieval output to be saved.
        config_dict : a dictionary containing the configuration variables
        prior_info :
        input_info : a dictionary containing the retrieval input variables
    
        Returns
        ------- 
        None
    '''
    Rd = 287
    K_C = 273.15

    # Make note of the error flag being used.
    error_flag = -9999

    # Set the RMSe threshold to specify whether or not a sample might be of suspect quality.
    rms_thres = 2.0 # K

    # This value is needed for profile extraction purposes.
    tot_num = len(output['temperature']) * 2

    # Create the current retrieval time's full filename accordingly.
    # need to pass this function the first datetime object of the file     
    if not os.path.isfile(out_filename):
        # Initialize the CDF file.
        rt_prf_grp = Dataset(out_filename, 'w', format='NETCDF3_CLASSIC')

        # Set the needed dimensions.
        time_dim = rt_prf_grp.createDimension('time', None)
        height_dim = rt_prf_grp.createDimension('height', len(output['height']))
        freq_dim = rt_prf_grp.createDimension('freq', len(output['Y']))
        dfs_dim = rt_prf_grp.createDimension('dfs', len(output['dfs']))
        arb_dim = rt_prf_grp.createDimension('arb', len(output['x_c']))
        exec_time_dim = rt_prf_grp.createDimension('exec_time',3)

        # Create the needed variables.
        prof_date = rt_prf_grp.createVariable('prof_date','i4')
        prof_time = rt_prf_grp.createVariable('prof_time','i4',('time',))
        base_time = rt_prf_grp.createVariable('base_time','i4')
        time_offset = rt_prf_grp.createVariable('time_offset','f8',('time',))
        hour = rt_prf_grp.createVariable('hour','f8','time')
        qc_flag = rt_prf_grp.createVariable('qc_flag','i2',('time',))
        rain_flag = rt_prf_grp.createVariable('rain_flag','i2',('time',))
        height = rt_prf_grp.createVariable('height','f4',('height',))
        temperature = rt_prf_grp.createVariable('temperature','f4',('time','height',))
        waterVapor = rt_prf_grp.createVariable('waterVapor','f4',('time','height',))
        lwp = rt_prf_grp.createVariable('lwp','f4',('time',))
        cbh = rt_prf_grp.createVariable('cbh','f4',('time',))
        cbh_flag = rt_prf_grp.createVariable('cbh_flag','i2',('time',))
        sigma_temperature = rt_prf_grp.createVariable('sigma_temperature','f4',('time','height',))
        sigma_waterVapor = rt_prf_grp.createVariable('sigma_waterVapor','f4',('time','height',))
        sigma_lwp = rt_prf_grp.createVariable('sigma_lwp','f4',('time',))
        converged_flag = rt_prf_grp.createVariable('converged_flag','i2',('time',))
        n_iter = rt_prf_grp.createVariable('n_iter','i2',('time',))
        exec_time = rt_prf_grp.createVariable('exec_time','i4',('time','exec_time',))
        rms = rt_prf_grp.createVariable('rms','f4',('time',))
        dfs = rt_prf_grp.createVariable('dfs','f4',('time','dfs',))
        sic = rt_prf_grp.createVariable('sic','f4',('time',))
        vres_temperature = rt_prf_grp.createVariable('vres_temperature','f4',('time','height',))
        vres_waterVapor = rt_prf_grp.createVariable('vres_waterVapor','f4',('time','height',))
        pressure = rt_prf_grp.createVariable('pressure','f4',('time','height',))
        freq = rt_prf_grp.createVariable('freq','f4',('freq',))
        freq_off = rt_prf_grp.createVariable('freq_offsets','f4',('freq',))
        rad_off = rt_prf_grp.createVariable('radiance_offsets','f4',('freq',))
        elev = rt_prf_grp.createVariable('elev','f4',('freq',))
        radiance_obs = rt_prf_grp.createVariable('radiance_obs','f4',('time','freq',))
        radiance_obs_uncertainty = rt_prf_grp.createVariable('radiance_obs_uncertainty','f4',('time','freq',))
        radiance_calculation = rt_prf_grp.createVariable('radiance_calculation','f4',('time','freq',))
        arb = rt_prf_grp.createVariable('arb','i2',('arb',))
        Xop = rt_prf_grp.createVariable('Xop','f4',('time','arb',))
        Sop = rt_prf_grp.createVariable('Sop','f8',('time','arb','arb',))
        Akernal = rt_prf_grp.createVariable('Akernal','f4',('time','arb','arb',))
        Xa = rt_prf_grp.createVariable('Xa','f4',('arb',))
        Sa = rt_prf_grp.createVariable('Sa','f8',('arb','arb',))
        lat = rt_prf_grp.createVariable('lat','f4')
        lon = rt_prf_grp.createVariable('lon','f4')
        alt = rt_prf_grp.createVariable('alt','f4')

        # Set the variable attributes.
        prof_date.long_name = 'The date of the retrieved profile(s)/MWR data collection'
        prof_date.format = 'YYYYMMDD'
        prof_time.long_name = 'The approximate time(s) of the retrieved profile(s)/MWR measurement.'
        prof_time.units = 'Hours,Minutes UTC'
        base_time.long_name = 'Epoch time'
        base_time.units = 'seconds since 1970/01/01 00:00:00+00:00 UTC'
        time_offset.long_name = 'Time offset from base-time'
        time_offset.units = 's'
        hour.long_name = 'Time'
        hour.units = 'Hours from 00:00 UTC'
        qc_flag.long_name = 'Manual QC flag'
        qc_flag.units = 'unitless'
        qc_flag.comment = 'A value of 0 implies quality is ok; non-zero values indicate that the sample has suspect quality.'
        qc_flag.value_descriptions = '0 = retrieval ok; 1 = no convergence; 2 = converged, but high rms value; 3 = converged, reasonable rms value, but rain detected'
        rain_flag.long_name = 'Rain flag (1 -> raining, 0 -> no rain detected)'
        rain_flag.units = 'unitless'
        rain_flag.comment = 'In general, a retrieval will be attempted even if rain is detected; but if rain is detected, note that the retrieval sample might not be completely accurate.'
        height.long_name = 'Height'
        height.units = 'km AGL'
        temperature.long_name = 'Temperature'
        temperature.units = 'C'
        waterVapor.long_name = 'Water Vapor Mixing Ratio'
        waterVapor.units = 'g/kg'
        lwp.long_name = 'Liquid Water Path'
        lwp.units = 'g/m^2'
        cbh.long_name = 'cloud base height'
        cbh.units = 'km AGL'
        cbh_flag.long_name = "Flag indicating the source of the cbh"
        cbh_flag.units = "unitless"
        cbh_flag.comment_0 = "Value 0 implies Clear Sky"
        cbh_flag.comment_1 = "Value 1 implies Inner Window"
        cbh_flag.comment_2 = "Value 2 implies Outer Window"
        cbh_flag.comment_3 = "Value 3 implies Default CBH"
        sigma_temperature.long_name = '1-sigma uncertainty in temperature'
        sigma_temperature.units = 'C'
        sigma_waterVapor.long_name = '1-sigma uncertainty in water vapor mixing ratio'
        sigma_waterVapor.units = 'g/kg'
        sigma_lwp.long_name = '1-sigma uncertainty in liquid water path'
        sigma_lwp.units = 'g/m^2'
        converged_flag.long_name = 'Convergence flag'
        converged_flag.units = 'unitless'
        converged_flag.comment = '1 indicates successful convergence, 0 implies no convergence'
        n_iter.long_name = 'Number of iterations performed'
        n_iter.units = 'unitless'
        exec_time.long_name = 'The execution time(s) of the retrieval sample(s)'
        exec_time.units = 'seconds'
        exec_time_comment = 'The execution time is the time it took for the time sample to be retrieved.'
        rms.long_name = 'Root mean square error between observed and computed spectrum'
        rms.units = 'K'
        dfs.long_name = 'Degrees of freedom of signal'
        dfs.units = 'unitless'
        dfs.comment = 'Total DFS, then DFS for each of temperature, waterVapor, and LWP.'
        sic.long_name = 'Shannon Information Content'
        sic.units = 'unitless'
        vres_temperature.long_name = 'Vertical resolution of the temperature profile'
        vres_temperature.units = 'km'
        vres_waterVapor.long_name = 'Vertical resolution of the water vapor profile'
        vres_waterVapor.units = 'km'
        pressure.long_name = 'Derived pressure'
        pressure.units = 'mb'
        pressure.comment = 'Derived from MWR surface pressure observation and the hypsometric equation using the thermodynamic profiles.'
        freq.long_name = 'Frequencies'
        freq.units = 'GHz'
        freq_off.long_name = 'Frequency offsets'
        freq_off.units = 'GHz'
        rad_off.long_name = 'MWR-observed brightness temperature offset'
        rad_off.long_name = 'K'
        elev.long_name = 'Elevations'
        elev.units = 'Deg.'
        radiance_obs.long_name = 'MWR-observed radiance'
        radiance_obs.units = 'K'
        radiance_obs_uncertainty.long_name = '1-sigma uncertainty in MWR-observed radiance'
        radiance_obs_uncertainty.units = 'K'
        radiance_calculation.long_name = 'Computed radiance (from the optimal state vector Xop)'
        radiance_calculation.units = 'K'
        arb.long_name = 'Arbitrary dimension'
        arb.units = 'mixed units'
        arb.comment = 'Contains: temperature profile (1), water vapor profile (2), and liquid cloud path (3).'
        Xop.long_name = 'Retrieved optimal solution/state vector'
        Xop.units = 'mixed units'
        Xop.comment = 'The first ' + str(len(output['temperature'])) + ' elements are temperature, the next ' + str(len(output['temperature'])) + ' elements are water vapor mixing ratio, and the last element is LWP.'
        Sop.long_name = 'Covariance matrix of the solution'
        Sop.units = 'mixed units'
        Akernal.long_name = 'Averaging kernal'
        Akernal.units = 'mixed units'
        Xa.long_name = 'Prior mean state vector'
        Xa.units = 'mixed units'
        Sa.long_name = 'Prior covariance matrix'
        Sa.units = 'mixed units'
        lat.long_name = 'MWR latitude'
        lat.units = 'degrees north'
        lon.long_name = 'MWR longitude'
        lon.units = 'degrees east'
        alt.long_name = 'MWR altitude'
        alt.units = 'm above MSL'

    else:
        # Re-initialize the CDF file.
        rt_prf_grp = Dataset(out_filename, 'r+', format='NETCDF3_CLASSIC')

        # Recover the existing variables.
        prof_date = rt_prf_grp.variables['prof_date']
        prof_time = rt_prf_grp.variables['prof_time']
        base_time = rt_prf_grp.variables['base_time']
        time_offset = rt_prf_grp.variables['time_offset']
        hour = rt_prf_grp.variables['hour']
        qc_flag = rt_prf_grp.variables['qc_flag']
        rain_flag = rt_prf_grp.variables['rain_flag']
        height = rt_prf_grp.variables['height']
        temperature = rt_prf_grp.variables['temperature']
        waterVapor = rt_prf_grp.variables['waterVapor']
        lwp = rt_prf_grp.variables['lwp']
        cbh = rt_prf_grp.variables['cbh']
        cbh_flag = rt_prf_grp.variables['cbh_flag']
        sigma_temperature = rt_prf_grp.variables['sigma_temperature']
        sigma_waterVapor = rt_prf_grp.variables['sigma_waterVapor']
        sigma_lwp = rt_prf_grp.variables['sigma_lwp']
        converged_flag = rt_prf_grp.variables['converged_flag']
        n_iter = rt_prf_grp.variables['n_iter']
        exec_time = rt_prf_grp.variables['exec_time']
        rms = rt_prf_grp.variables['rms']
        dfs = rt_prf_grp.variables['dfs']
        sic = rt_prf_grp.variables['sic']
        vres_temperature = rt_prf_grp.variables['vres_temperature']
        vres_waterVapor = rt_prf_grp.variables['vres_waterVapor']
        pressure = rt_prf_grp.variables['pressure']
        freq = rt_prf_grp.variables['freq']
        freq_off = rt_prf_grp.variables['freq_offsets']
        rad_off = rt_prf_grp.variables['radiance_offsets']
        elev = rt_prf_grp.variables['elev']
        radiance_obs = rt_prf_grp.variables['radiance_obs']
        radiance_obs_uncertainty = rt_prf_grp.variables['radiance_obs_uncertainty']
        radiance_calculation = rt_prf_grp.variables['radiance_calculation']
        arb = rt_prf_grp.variables['arb']
        Xop = rt_prf_grp.variables['Xop']
        Sop = rt_prf_grp.variables['Sop']
        Akernal = rt_prf_grp.variables['Akernal']
        Xa = rt_prf_grp.variables['Xa']
        Sa = rt_prf_grp.variables['Sa']
        lat = rt_prf_grp.variables['lat']
        lon = rt_prf_grp.variables['lon']
        alt = rt_prf_grp.variables['alt']

    #####################################################################################################################
    # Set variables in the netCDF file that describes the file, algorithm version, and any VIP settings set by the user #
    #####################################################################################################################

    rt_prf_grp.Algorithm_code = config_dict['alg_name']
    rt_prf_grp.Algorithm_author = config_dict['alg_authors']
    rt_prf_grp.Algorithm_version = "MWRoe v2" ##
    rt_prf_grp.Algorithm_reference = config_dict['alg_ref']

    rt_prf_grp.Datafile_created_on_date = tm.ctime(tm.time())
    rt_prf_grp.Datafile_created_on_machine = config_dict['mac_id_s']

    rt_prf_grp.Site = config_dict['mwr_site']
    rt_prf_grp.Instrument = config_dict['mwr_inst']
    rt_prf_grp.Dataset_contact = config_dict['alg_contacts']


    rt_prf_grp.Prior_dataset_comment = prior_info['prior_comment']
    rt_prf_grp.Prior_dataset_filename = prior_info['prior_filename']
    rt_prf_grp.Prior_dataset_number_profiles = prior_info['prior_numprofs']

    rt_prf_grp.Total_clock_execution_time_in_s = "1" #str(sum_runtime_s)
    rt_prf_grp.General_comment = "If for a given time sample, all solution/related values stored are shown as " + str(error_flag) + ", then the retrieval for that time sample failed as it was not able to converge properly."
    rt_prf_grp.Radiance_comment = "Note that all radiance values are given as brightness temperatures with units of degrees Kelvin (K)."

    rt_prf_grp.CF_mwr_type = config_dict['mwr_type']
    rt_prf_grp.CF_monortm_home = config_dict['monortm_exec_path']
    rt_prf_grp.CF_jac_option = str(config_dict['jac_type_flag'])
    rt_prf_grp.CF_jac_max_ht = str(round(config_dict['jac_max_alt'],2))
    rt_prf_grp.CF_first_guess = config_dict['use_prior_init']

    #################################################################################################################
    # Do any final calculations/processing of the input parameters before the final values are written to the file. #
    #################################################################################################################

    # Calculate needed additional time fields.
    utc_zone = tz.gettz('UTC')

    # Get base (epoch) time.
    base_tm_ob = input_info['dt_times'][0].replace(tzinfo=utc_zone)
    epoch_time = long(date2num(base_tm_ob,'seconds since 1970-01-01 00:00:00+00:00'))

    # Get time offset from base time.
    toff_tm_ob = input_info['dt_times'][output['sample_index']].replace(tzinfo=utc_zone)
    toff_time = long(date2num(toff_tm_ob,'seconds since 1970-01-01 00:00:00+00:00'))
    ep_time_offset = abs(toff_time - epoch_time)

    # Convert the altitude grid to km for storage consistency.
    alta_km = output['height']

    # Set the current QC flag.
    if output['converged_flag'] == 0:
        qc_flag_val = 1
    elif output['rms'] >= rms_thres: # If the RMS output was greater than the subjectively set "bad" threshold
        qc_flag_val = 2
    elif input_info['rainflags'][output['sample_index']] == 1:  # If it was raining
        qc_flag_val = 3
    else: # If it's none of these, then the retrieval is OK.
        qc_flag_val = 0

    # Make the variable markers for each element in the X_c matrix
    arbitrary = np.empty(len(output['x_c']))
    arbitrary[:tot_num/2] = 1
    arbitrary[tot_num/2:tot_num] = 2
    arbitrary[-1] = 3

    # Recover the needed uncertainty vectors and assign them to variables
    err_vars_op = np.diag(output['Sop'])
    num_alts = tot_num / 2
    err_stds_op = np.sqrt(err_vars_op)
    Terr_op = err_stds_op[0:num_alts]
    Qerr_op = err_stds_op[num_alts:tot_num]
    LWPerr_op = err_stds_op[tot_num]

    # Develop the frequency offset stuff (THIS IS PROBABLY THE LAST THING I NEED TO FIX.)
    freq_offsets = np.zeros(len(output["Y"]))
    tb_offsets = np.zeros(len(output["Y"]))
    freq_offsets[:len(input_info['z_freqs'])] = config_dict['freq_offsets']
    tb_offsets[:len(input_info['z_freqs'])] = config_dict['tb_offsets']

    # Separate out the diagonal of the Se matrix for storage
    Y_err_vars = np.diag(output['Se'])
    Y_err = np.sqrt(Y_err_vars)

    # Retrieve the current number of time stamps in the current netCDF file.
    ntimes = len(rt_prf_grp.dimensions['time'])

    # Save the retrieval inputs and outputs to the file.
    prof_date[:] = long(datetime.strftime(toff_tm_ob, '%Y%m%d'))

    ###################################################
    # Set values to the variables in the netCDF file. #
    ###################################################
    
    dec_hour = helper.convert_time(datetime.strftime(toff_tm_ob, '%H%M'), 'h.hf')
    
    prof_time[ntimes:ntimes+1] = datetime.strftime(toff_tm_ob, '%H%M')
    base_time[:] = epoch_time
    time_offset[ntimes:ntimes+1] = ep_time_offset
    time_offset[len(time_offset) - 1] = ep_time_offset
    hour[ntimes:ntimes+1] = dec_hour
    hour[len(hour) - 1] = dec_hour
    qc_flag[ntimes:ntimes+1] = qc_flag_val
    qc_flag[len(hour) - 1] = qc_flag_val
    rain_flag[ntimes:ntimes+1] = input_info['rainflags'][output['sample_index']]
    rain_flag[len(hour) - 1] = input_info['rainflags'][output['sample_index']]
    height[:] = alta_km
    temperature[ntimes:ntimes+1,:] = output['temperature']
    waterVapor[ntimes:ntimes+1,:] = output['waterVapor'] # already solved for
    lwp[ntimes:ntimes+1] = output['LWP'] # already solved for
    lwp[len(lwp) - 1] = float(output['LWP'])
    cbh[ntimes:ntimes+1] = output['cbh']
    cbh[len(cbh) - 1] = output['cbh']
    cbh_flag[ntimes:ntimes+1] = int(output['cbh_flag'])
    cbh_flag[len(cbh_flag) - 1] = int(output['cbh_flag'])
    sigma_temperature[ntimes:ntimes+1,:] = Terr_op # already solved for
    sigma_waterVapor[ntimes:ntimes+1,:] = Qerr_op # already solved for
    sigma_lwp[ntimes:ntimes+1] = LWPerr_op   # already solved for
    sigma_lwp[len(sigma_lwp) - 1] = LWPerr_op
    converged_flag[ntimes:ntimes+1] = output['converged_flag']
    converged_flag[len(converged_flag) - 1] = output['converged_flag']
    n_iter[ntimes:ntimes+1] = output['iter_count']
    n_iter[len(n_iter) - 1] = output['iter_count']
    exec_time[ntimes:ntimes+1,:] = output['retr_time']
    exec_time[len(exec_time) - 1, :] = output['retr_time']
    rms[ntimes:ntimes+1] = output['rms']
    rms[len(rms) - 1] = output['rms']
    dfs[ntimes:ntimes+1,:] = output['dfs']
    sic[ntimes:ntimes+1] = output['sic']
    sic[len(sic) - 1] = output['sic']
    vres_temperature[ntimes:ntimes+1,:] = output['T_vres'] / 1000.
    vres_waterVapor[ntimes:ntimes+1,:] = output['Q_vres'] / 1000.
    pressure[ntimes:ntimes+1,:] = output['pressure'] # already solved for
    freq[:] = input_info['all_freqs'] 
    freq_off[:] = freq_offsets 
    rad_off[:] = tb_offsets
    elev[:] = input_info['all_elevs'] 
    radiance_obs[ntimes:ntimes+1,:] = output['Y'] # already extracted
    radiance_obs_uncertainty[ntimes:ntimes+1,:] = Y_err # already solved for
    radiance_calculation[ntimes:ntimes+1,:] = output['F_xop'] # already extracted
    arb[:] = arbitrary 
    Xop[ntimes:ntimes+1,:] = output['x_c'] # already extracted
    Sop[ntimes:ntimes+1,:,:] = output['Sop'] # already extracted
    Akernal[ntimes:ntimes+1,:] = output['Akernal_op'] # already extracted
    Xa[:] = output['Xa'] # already extracted
    Sa[:,:] = output['Sa'] # already squeezed
    lat[:] = input_info['lat']
    lon[:] = input_info['lon']
    alt[:] = input_info['alt']

    # Close the netCDF file.
    rt_prf_grp.close()



def writeMonoRTMFreqs(oe_input, config, offsets):
    '''
        writeMonoRTMFreqs

        This function writes the frequency file used to run the 
        MonoRTM using the frequencies specified by the retrieval.

        Parameters
        ----------
        oe_input : a dictionary containing the retrieval input variables
        config : a dictionary containing all of the retrieval configurationv variables

        Returns
        -------
        None
    '''
    # Write the ZENITH frequencies we want on the frequency file.
    freq_zenith_file = config['working_dir'] + "/" + config['monortm_freqs_fname_base'] + '.zen'
    num_zenith_freqs = str(len(oe_input['z_freqs']))
    zenith_freqs = [str(i + o) for o,i in zip(offsets['z_freq_offsets'], oe_input['z_freqs'])]

    freq_file = open(freq_zenith_file,'w')
    freq_file.write('\n')
    freq_file.write('%s\n' % num_zenith_freqs)

    for i in range(int(num_zenith_freqs)):
        #print '%s\n' % zenith_freqs[i]
        freq_file.write('%s\n' % zenith_freqs[i])

    freq_file.close()
    monoRTM_freq_files = [freq_zenith_file]

    # Write the OFF-ZENITH frequencies we want to a separate frequency file.
    num_off_zenith_freqs = str(len(oe_input['oz_freqs']))
    if int(num_off_zenith_freqs) != 0:
        freq_oz_file = config['working_dir'] + "/" + config['monortm_freqs_fname_base'] + '.ozen'
        oz_freqs = oe_input['oz_freqs']
        off_zenith_freqs = [str(i + o) for o,i in zip(offsets['oz_freq_offsets'],oz_freqs)]
        freq_file = open(freq_oz_file,'w')
        freq_file.write('\n')
        freq_file.write('%s\n' % num_off_zenith_freqs) # Write the number of frequencies for the monoRTM to the freq file

        # Write the off-zenith frequencies to the monortm_freqs file
        for i in range(int(num_off_zenith_freqs)):
            freq_file.write('%s\n' % off_zenith_freqs[i])
        freq_file.close()

        monoRTM_freq_files.append(freq_oz_file)

    return monoRTM_freq_files



def writeMonoRTMConfig(alt, config_dict):
    """
        writeMonoRTMConfig

        This function accepts a height grid (km), along with with the config dictionary.
        Using the height grid, this function will write the height grid to the 
        monortm_config file and set it as an evironmental variable.

        Parameters
        ----------
        alt : the height grid (km)
        config_dict : the VIP config dictionary made by reader.readVip()

        Returns
        -------
        Nothing

    """
    # Write the MonoRTM config. file with the altitude grid from the prior file.
    if np.max(alt) < 70:
        alt_km = np.hstack((alt, np.arange(np.max(alt)+5, 75, 5)))

    num_alts = len(alt_km)
    n_s = str(len(alt_km))
    alt_s = [str(i) for i in alt_km] # Make the array containing the heights for the height grid from alta_km

    monortm_configs_fname = config_dict['working_dir'] + '/' + config_dict['monortm_configs_fname'] + '_' + n_s + 'L.txt'
    mcf_fncp = monortm_configs_fname
    print mcf_fncp
    mcf_file = open(mcf_fncp,'w')
    mcf_file.write(config_dict['monortm_path'] + '\n')
    mcf_file.write(config_dict['monortm_spectral_dat'] + '\n')
    mcf_file.write('0       Verbose level (0 is quiet, 1 is DEBUG)\n')
    mcf_file.write('6       Default atmos (1->trop, 2->MLS, 3->MLW, 4->SAS, 5->SAW, 6->USStd)\n')
    mcf_file.write('1       Ispd (0->use all lines and go slow, 1->use subset and go faster)\n')
    mcf_file.write('1.000       Self-broadened WV continuum multiplier\n')
    mcf_file.write('1.000       Foreign-broadened WV continuum multiplier\n')
    mcf_file.write('1.000       CO2 continuum multiplier\n')
    mcf_file.write('1.000       O3 continuum multiplier\n')
    mcf_file.write('1.000       O2 continuum multiplier\n')
    mcf_file.write('1.000       N2 continuum multiplier\n')
    mcf_file.write('15.0        Minimum altitude [km MSL] that sonde must achieve to be used, Number of model levels (next line), Heights in km above minimum input sonde level (all lines following), and a blank line at the end\n')

    mcf_file.write('%s\n'%n_s) # Write the number of levels for the height grid
    for i in range(num_alts):
         mcf_file.write('%s\n'%alt_s[i])
    mcf_file.write('\n')
    mcf_file.close()

    # Set the monortm_config enviromental variable
    os.system ("chmod 744" + " " + mcf_fncp)
    os.environ['monortm_config'] = mcf_fncp

    return mcf_fncp


def makeMonoRTMCDF(sonde_file,alta,T_z,P_z,RH_z):
    '''
        makeMonoRTMCDF
    
        This function creates a netCDF file containing the thermodynamic profile
        that will be used to run the MonoRTM.  The netCDF file is formatted off of
        the ARM radiosonde format, as Dave Turner's (NSSL) MonoRTM wrapper script
        that is called by the MWRoe program requires a netCDF file in this format.

        Parameters
        ----------
        sonde_file : a string that is the sonde file name to be written to.
        alta : an array containing the height grid [meters]
        T_z : an array containing the temperature profile [C]
        P_z : an array containing the pressure profile [mb]
        RH_z : an array contianing the relative humidity profile [%]

        Returns
        -------
        None
    '''

    # Initialize the CDF file.
    sf_tmp_grp = Dataset(sonde_file, 'w', format='NETCDF3_CLASSIC')

    # Set the needed dimensions.
    time_dim = sf_tmp_grp.createDimension('time', None)

    # Create the needed variables.
    base_time = sf_tmp_grp.createVariable('base_time','i4')
    time_offset = sf_tmp_grp.createVariable('time_offset','f8',('time',))
    pres = sf_tmp_grp.createVariable('pres','f4',('time',))
    tdry = sf_tmp_grp.createVariable('tdry','f4',('time',))
    rh = sf_tmp_grp.createVariable('rh','f4',('time',))
    dp = sf_tmp_grp.createVariable('dp','f4',('time',))

    wspd = sf_tmp_grp.createVariable('wspd','f4',('time',))
    deg = sf_tmp_grp.createVariable('deg','f4',('time',))
    lat = sf_tmp_grp.createVariable('lat','f4',('time',))
    lon = sf_tmp_grp.createVariable('lon','f4',('time',))
    alt = sf_tmp_grp.createVariable('alt','f4',('time',))

    # Set the needed attributes.
    sf_tmp_grp.Description = 'Temporary / artificial radiosonde file to generate MonoRTM calcs for the MWR T/Q retrievals'
    sf_tmp_grp.Primary_source = 'MWR thermodynamic retrieval primary function: compute_thermodynamic_mwr_retrievals'
    sf_tmp_grp.Secondary_source = 'MWR thermodynamic retrieval secondary function: create_temp_sonde_file'
    sf_tmp_grp.Time_last_created = tm.ctime(tm.time())
    base_time.long_name = 'Time since 1 Jan 1970 00:00:00 UTC'
    base_time.units = 'seconds'
    time_offset.long_name = 'Time offset from base_time'
    time_offset.units = 'seconds'
    pres.long_name = 'Pressure'
    pres.units = 'hPa'
    tdry.long_name = 'Dry bulb temperature'
    tdry.units = 'C'
    rh.long_name = 'Relative humidity'
    rh.units = '%'
    dp.long_name = 'Dewpoint temperature'
    dp.units = 'C'
    wspd.long_name = 'Wind speed'
    wspd.units = 'm/s'
    deg.long_name = 'Wind direction'
    deg.units = 'Degrees'
    lat.long_name = 'Latitude'
    lat.units = 'degrees North'
    lon.long_name = 'Longitude'
    lon.units = 'degrees East'
    alt.long_name = 'Altitude above mean sea level'
    alt.units = 'm'

    alt_nums = np.array(range(len(alta)))

    # Set the variables using the input profiles.
    base_time[:] = int(tm.time())
    time_offset[:] = 2 * alt_nums
    pres[:] = P_z
    tdry[:] = T_z
    rh[:] = RH_z
    dp[:] = 1.0 * alt_nums
    wspd[:] = 1.0 * alt_nums
    deg[:] = 1.0 * alt_nums
    lat[:] = 1.0 * alt_nums
    lon[:] = 1.0 * alt_nums
    alt[:] = alta

    # Close the CDF file.
    sf_tmp_grp.close()

    # Set the permissions for the CDF file.
    os.system ("chmod 744" + " " + sonde_file)


def writeMonoRTMFreqs_FM(z_freqs,oz_freqs, config):
    '''
        writeMonoRTMFreqs

        This function writes the frequency file used to run the 
        MonoRTM using the frequencies specified by the retrieval.

        Parameters
        ----------
        oe_input : a dictionary containing the retrieval input variables
        config : a dictionary containing all of the retrieval configurationv variables

        Returns
        -------
        None
    '''
    # Write the ZENITH frequencies we want on the frequency file.
    freq_zenith_file = config['working_dir'] + "/" + config['monortm_freqs_fname_base'] + '_FM.zen'
    num_zenith_freqs = str(len(z_freqs))
    zenith_freqs = [str(i) for i in z_freqs]

    freq_file = open(freq_zenith_file,'w')
    freq_file.write('\n')
    freq_file.write('%s\n' % num_zenith_freqs)

    for i in range(int(num_zenith_freqs)):
        #print '%s\n' % zenith_freqs[i]
        freq_file.write('%s\n' % zenith_freqs[i])

    freq_file.close()
    monoRTM_freq_files = [freq_zenith_file]

    # Write the OFF-ZENITH frequencies we want to a separate frequency file.
    num_off_zenith_freqs = str(len(oz_freqs))
    if int(num_off_zenith_freqs) != 0:
        freq_oz_file = config['working_dir'] + "/" + config['monortm_freqs_fname_base'] + '_FM.ozen'
        off_zenith_freqs = [str(i) for i in oz_freqs]
        freq_file = open(freq_oz_file,'w')
        freq_file.write('\n')
        freq_file.write('%s\n' % num_off_zenith_freqs) # Write the number of frequencies for the monoRTM to the freq file

        # Write the off-zenith frequencies to the monortm_freqs file
        for i in range(int(num_off_zenith_freqs)):
            freq_file.write('%s\n' % off_zenith_freqs[i])
        freq_file.close()

        monoRTM_freq_files.append(freq_oz_file)

    return monoRTM_freq_files


