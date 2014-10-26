from netCDF4 import Dataset, num2date, date2num
import sys
import glob
import numpy as np
from datetime import datetime
import os

def constructPrior(prior_filename, config_dict):
    """
        constructPrior

        This function is called by the main script and does the following things:
            1.) Opens the prior file specified by the user.
            2.) Loads in the height grid from the prior file.
            3.) Sets the Xa numpy array that is passed to the retrieval (thermo profile + LWP_mean)
            4.) Sets the Sa matrix that is passed to the retrieval (thermo profile + LWP_std)

        Parameters
        ----------
        prior_filename : the path to the prior .cdf file.
        config_dict : the dictionary loaded in from readVIP

        Returns
        -------
        Xa : a numpy array of the mean state vector X
        Sa : a numpy matrix of the covariance around the state vector X
        alta : the height grid for the retrieval (meters) of length (len(Xa)-1)/2
        prior_info : a dictionary/structure containing information about the prior file.
    """

    # Extract/form the needed (month dependent) climatology/prior parameters x_a and Sa. Also extract the altitude grid.
    # x_a is a vector of size n x 1, Sa is a matrix of size n x n.
    # Load in the prior (Xa, Sa):
    prior = Dataset(prior_filename)
    alta = prior.variables['height'][:]*1000.
    Xa = prior.variables['mean_prior'][:]
    Xa = np.concatenate((Xa, [config_dict['mean_prior_LWP']]))

    # Make the Se matrix
    Sa_pre = prior.variables['covariance_prior'][:]
    Sa = np.zeros((len(Xa),len(Xa)))
    Sa[:len(Xa)-1,:len(Xa)-1] = Sa_pre
    Sa[-1,-1] = config_dict['std_prior_LWP']**2
    Sa = np.matrix(Sa)

    # Load additional info from the prior file.
    prior_comment = str(getattr(prior,'Comment', "No comment specified"))
    num_profsc = str(getattr(prior,'Nsonde', "Number of profiles used to generate the prior was not specified"))
    prior_numprofs = [s for s in num_profsc.split() if s.isdigit()][0]
    if len(prior_numprofs) == 0:
        prior_numprofs = -9999.
    
    prior_info = {}
    prior_info['prior_filename'] = prior_filename
    prior_info['prior_comment'] = prior_comment
    prior_info['prior_numprofs'] = prior_numprofs
                
    # Close the prior file.
    prior.close()

    return Xa, Sa, alta, prior_info

def findVIPVariable(variable, vip_string):
    """
        findVIPVariable

        This function will search for the argument "variable" contained
        within the VIP file, which is passed to the code via the vip_string
        argument.

        Parameters
        ----------
        variable : the variable name to search for in the VIP file
        vip_string : the already parsed VIP file in string format

        Returns
        -------
        either a string for a float of the value assigned to the variable

        Returns None if the variable name is not found.
    """    
    #This function searches for the value associatied with the key
    #(key is variable) within the VIP file and returns it either as a
    #float or a string

    ini=vip_string.find(variable)+(len(variable)+1)
    rest=vip_string[ini:]
    search_enter=rest.find('\n')
    data = (rest[:search_enter])
    if len(data) == 0:
        #print "VIP variable: " + variable + " not found in VIP file."
        return None
    datas = data.split(";")
    string = datas[0].replace('=', '')
    try: 
        return float(string)
    except:
        return string

def readVIP(vip_fn):
    """
        readVIP

        This function will load the configuration data/variables from the VIP text file provided
        by the user of the retrieval.  This VIP file is very similar to the ones used in AERIoe.

        In this code, several different things are loaded from the VIP file and machine that
        are used to run the program.  In the following code, the different categories
        of things that are read in are commented below.

        Parameters
        ----------
        vip_fn : the filename to the VIP file

        Returns
        -------
        config_dict : a dictionary/structure of VIP file, MWRoe code, and machine information
    """
    
    # Read the VIP file.
    vip_file = open(vip_fn,'r')
    vip_string = vip_file.read()
    vip_file.close()
    
    vip_string = vip_string.split('\n')
    vip_string[:] = [i for i in vip_string if i != '']
    vip_string[:] = [i for i in vip_string if not i.startswith('#')] # should now have length == 55
    vip_string = '\n'.join(vip_string)

    config_dict = {}    
    
    ############################################################################
    #
    #   VIP information about zenith obs to be used in the retrieval
    #
    ############################################################################
    var = findVIPVariable('zenith_freqs', vip_string).strip()

    # 1) Load in the zenith frequencies specified in the VIP file.
    try:
        # Try to read in the zenith frequencies
        var = findVIPVariable('zenith_freqs', vip_string).strip()
    except:
        # If there's an error, stop the program and tell the user.
        print "Error reading the zenith frequencies."
        sys.exit()

    # Set the zenith frequencies
    config_dict['zenith_freqs'] = var.split(',')
    config_dict['zenith_freqs'] = np.asarray(config_dict['zenith_freqs'], dtype=float)

    # 2) Load in the Tb uncertainties specified in the VIP file for the zenith Tb observations.
    try:
        # Try to read in the zenith frequency uncertantities
        var = findVIPVariable('zenith_noise', vip_string).strip()
    except:
        # If there's an error, stop the program and tell the user.
        print "Error reading the zenith frequency uncertanities."
        sys.exit()
    
    # Check to see if there are an equal number of elements between the zenith frequencies and the zenith uncertaintites in the VIP file.
    var = np.asarray(var.split(','), dtype=float)
    if len(config_dict['zenith_freqs']) != len(var):
        # This means that the user hasn't specified the right number of zenith uncertanities
        print "Number of zenith channel uncertanities listed in VIP file (" + str(len(var)) + ") does not match the number of specified zenith channel frequencies (" + str(len(config_dict['zenith_freqs'])) + "). Aborting."
        sys.exit()
    else:
        config_dict['zenith_uncert'] = var

    ############################################################################
    #
    #   VIP information about Off-zenith obs to be used in the retrieval
    #
    ############################################################################
    
    # 1.) Load in the off-zenith frequencies and the elevation angles used for the retrieval
    try:
        var = findVIPVariable('off_zenith_freq', vip_string)
        var2 = findVIPVariable('off_zenith_elev', vip_string)
    except:
        print "Error reading the off-zenith frequenices from VIP file."
        sys.exit()
    
    # 2.) Check to see if the VIP fields actually means don't use the off-zenith frequencies in the retrieval.
    if var == 0 or var2 == 0:
        config_dict['off_zenith_freqs'] = []
        config_dict['off_zenith_elevs'] = []
    else:
        config_dict['off_zenith_freqs'] = var.strip().split(',')
        config_dict['off_zenith_freqs'] = np.asarray(config_dict['off_zenith_freqs'], dtype=float)
        config_dict['off_zenith_elevs'] = var2.strip().split(',')
        config_dict['off_zenith_elevs'] = np.sort(np.asarray(config_dict['off_zenith_elevs'], dtype=float))

    # 3.) Load in the off-zenith frequency uncertainities specified in the VIP file.
    if var != 0 or var2 != 0:
        # Check to see if there are off-zenith frequency uncertanities listed (only happens if we have off zenith frequencie listed)
        try:
            # Try to read in the off_zenith frequency uncertantities
            var = findVIPVariable('off_zenith_noise', vip_string).strip()
        except:
            # If there's an error, stop the program and tell the user.
            print "Error reading the off-zenith frequency uncertanities."
            sys.exit()
      
        # Gather the off-zenith channel uncertanities since they exist and put them into an array
        var = np.asarray(var.split(','), dtype=float)
        if len(config_dict['off_zenith_freqs']) != len(var):
            # This means that the user hasn't specified the right number of zenith uncertanities
            print "Number of off-zenith channel uncertanities listed in VIP file does not match the number of specified off-zenith channel frequencies. Aborting."
            sys.exit()
        else:
            config_dict['off_zenith_uncert'] = var
    else:
        config_dict['off_zenith_uncert'] = []

   
    ############################################################################
    #
    #   VIP information about frequency offsets to account for lack of band-pass
    #
    ############################################################################

    # 1.) Check for frequency offsets (ONLY FOR ZENITH OBS) listed in the VIP file.
    try: 
        var = findVIPVariable('freq_offset', vip_string)
    except:
        print "Error reading frequency offsets from VIP file."
        sys.exit()
 
    # 2.) Load them in if they exist.
    if var == 0:
        # Means that the user does not want any frequency offsets
        config_dict['freq_offsets'] = np.zeros(config_dict['zenith_freqs'].shape[0])
        print "Using no offsets for the zenith frequencies."
    elif type(var) == str and len(var.strip().split(',')) != len(config_dict['zenith_freqs']):
        # This part checks that the user listed the same number of zenith frequencies and frequency offsets in the VIP file.
        print "Number of frequency offsets listed does not match the number of zenith frequencies listed.  Aborting."
        sys.exit()
    else:
        config_dict['freq_offsets'] = np.asarray(var.split(','), dtype=float)
        print "Using offsets for the zenith frequencies."


    ############################################################################
    #
    #   VIP information about frequency offsets to account for Tb Biases
    #
    ############################################################################
   
    # 1.) Check for brightness temperature offsets (ONLY FOR ZENITH OBS) listed in the VIP file.
    try:
        var = findVIPVariable('tb_offset', vip_string)
    except:
        print "Error reading observed brightness temperature offsets from VIP file."
        sys.exit()
    if type(var) == str:
        var = var.strip()

    # 2.) Load them in if they exist.
    if var == 0:
        config_dict['tb_offsets'] = np.zeros(config_dict['zenith_freqs'].shape[0])
        print "Using no offsets for the observed zenith brightness temperature calculations."
    elif type(var) == str and len(var.split(',')) != len(config_dict['zenith_freqs']):
        # This means that the user hasn't given us the correct number of brightness temperature offsets
        print "Number of TB offsets listed does not match the number of zenith frequencies listed.  Aborting."
        sys.exit()
    else:
        config_dict['tb_offsets'] = np.asarray(var.split(','), dtype=float)
        print "Using offsets for the observed zenith brightness temperature calculations."

    # Compute the total expected length of each Y vector.
    config_dict['num_Y_elems'] = len(config_dict['zenith_freqs']) + (len(config_dict['off_zenith_freqs']) * len(config_dict['off_zenith_elevs']))

    ############################################################################
    #
    #   VIP information about the Jacobian calculations
    #
    ############################################################################   
    
    #Retrieved variable (T,Q,LWP) perturbation settings for finite differencing Jacobian
    config_dict['jac_Tpert'] = 1.0
    config_dict['jac_Qpert'] = .99
    config_dict['jac_LWPpert'] = 10.0
    
    #Jacobian settings
    config_dict['jac_type_flag'] = int(findVIPVariable("jac_option", vip_string)) # Not working yet.
    config_dict['jac_max_alt'] = float(findVIPVariable('jac_max_ht', vip_string))

    ############################################################################
    #
    #   VIP information about the Prior (Sa, Xa)
    #
    ############################################################################   

    # Load in the prior data mean and standard deviation from the VIP file
    try:
        config_dict['mean_prior_LWP'] = float(findVIPVariable('prior_lwp_mn', vip_string)) #float(vip_string[30])
    except:
        config_dict['mean_prior_LWP'] = 0
    try:    
        config_dict['std_prior_LWP'] = float(findVIPVariable('prior_lwp_sd', vip_string)) #float(vip_string[31])
    except:
        config_dict['std_prior_LWP'] = 10.0

    config_dict['use_prior_init'] = int(findVIPVariable('first_guess', vip_string))

    ############################################################################
    #
    #   VIP information regarding the microwave radiometer properties and any 
    #   corrections that should be made.
    #
    ############################################################################   

    # MWR Data corrections
    config_dict['mwr_type'] = int(findVIPVariable('mwr_type', vip_string))
    config_dict['mwr_lat'] = findVIPVariable('mwr_lat', vip_string)
    config_dict['mwr_lon'] = findVIPVariable('mwr_lon', vip_string)
    config_dict['mwr_alt'] = findVIPVariable('mwr_alt', vip_string)
    config_dict['mwr_calib_pres'] = np.asarray(findVIPVariable('mwr_calib_pres',vip_string).strip().split(','), dtype=float)
    config_dict['mwr_site'] = findVIPVariable('globatt_Site', vip_string).strip()
    config_dict['mwr_inst'] = findVIPVariable('globatt_Instrument', vip_string).strip()

    ############################################################################
    #
    #   VIP information about how the retrieval should handle clouds and ceilometer data.
    #
    ############################################################################   

    # Ceilometer information
    config_dict['cbh_type'] = int(findVIPVariable('cbh_type', vip_string))
    config_dict['cbh_path'] = findVIPVariable('cbh_path', vip_string).strip()
    config_dict['cbh_window_in'] = int(findVIPVariable('cbh_window_in', vip_string))
    config_dict['cbh_window_out'] = int(findVIPVariable('cbh_window_out', vip_string))
    config_dict['cbh_default_ht'] = float(findVIPVariable('cbh_default_ht', vip_string))

 
    
    ############################################################################
    #
    #   VIP information about the MonoRTM executable and files
    #
    ############################################################################   

    #MonoRTM settings
    config_dict['monortm_freqs_fname_base'] = 'mwr_retr_monortm_freqs'
    config_dict['monortm_configs_fname'] = 'mwr_retr_monortm_config'
    config_dict['num_monortm_calc_header_lines'] = 11 #int(vip_string[37])
    config_dict['monortm_exec_path'] = findVIPVariable('monortm_exec', vip_string).strip()
    
    #Working directory for MonoRTM
    config_dict['working_dir'] = findVIPVariable('working_dir', vip_string).strip()
    config_dict['monortm_path'] = findVIPVariable('monortm_path', vip_string).strip()
    config_dict['monortm_spectral_dat'] = findVIPVariable('monortm_spectral_dat', vip_string).strip()

    ############################################################################
    #
    #   VIP information about input and output paths for the retrieval
    #
    ############################################################################  

    #Paths to MWR data 
    config_dict['mwr_dir'] = findVIPVariable('mwr_path', vip_string).strip()
    
    #Path to retrieved data
    config_dict['output_path'] = findVIPVariable('output_path', vip_string).strip()
    config_dict['output_clobber'] = int(findVIPVariable('output_clobber', vip_string))
    config_dict['output_rootname'] = findVIPVariable('output_rootname', vip_string).strip()

    ############################################################################
    #
    #   Information about the MWRoe Code to be embedded in the retrieval output.
    #
    ############################################################################ 

    #Algorithm information
    config_dict['alg_authors'] = 'Greg Blumberg, Stephen Castleberry, Dave Turner'
    config_dict['alg_contacts'] = 'wblumberg@ou.edu'
    config_dict['alg_name'] = 'MWRoe Version 2'
    config_dict['alg_ref'] = 'Blumberg et. al. 2015 JAMC'
    config_dict['dataset_contact'] = findVIPVariable('globatt_Dataset_contact', vip_string).strip()
    # Get the identifier of the machine on which the algorithm is being run.  
    config_dict['mac_id_s'] = '-'.join(os.uname())

    ############################################################################
    #
    #   VIP information on the max number of iterations the retrieval should make.
    #
    ############################################################################ 

    config_dict['max_iterations'] = int(findVIPVariable('max_iterations', vip_string))


    # Return the dictionary (or structure) containing all of the configuration information needed to run the retrieval.        
    return config_dict

def read_vceil_data(config, date, times):
    """
        read_vceil_data

        This function performs the following functions:
            a.) Reads in the ceilometer data (if it exists)
            b.) Splits the data into portions specified by the times array
            c.) Provides best estimates of the cloud base height via the innerWindow and outerWindow specified by the VIP file.
            d.) Returns an array

        Parameters
        ----------
        config : the config_dict returned by readVIP
        date : the date of the retrieval time in YYYYMMDD format.
        times : an array of the times MWRoe will be retrieving in epoch time.

        Returns
        -------
        cbh : the cloud base height array (same length as times) (km)
        cbh_flag : an array the same length as cbh with flags indicating how the cloud base height was picked.
    """

    vceil_files = glob.glob(config['cbh_path'] + '/*' + str(date[0]) + '*.cdf')
    inner_window = config['cbh_window_in']
    outer_window = config['cbh_window_out']

    if config['cbh_type'] == 1 and len(vceil_files) != 0: # ARM Vceil
        vceil_cdf = Dataset(vceil_files[0])
        cloudHeight = vceil_cdf.variables['first_cbh'][:]/1000.
        no_cloud_idx = np.where(cloudHeight < 0)[0]
    
        epoch_times = vceil_cdf.variables['base_time'][:] + vceil_cdf.variables['time_offset'][:]
    elif config['cbh_type'] == 2 and len(vceil_files) != 0: # Greg's Ceilometer files
        vceil_cdf = Dataset(vceil_files[0])
        cloudHeight = vceil_cdf.variables['cloudHeight'][:] # km AGL

        epoch_times = vceil_cdf.variables['base_time'][:] + vceil_cdf.variables['time_offset'][:]
    else:
        print "VCEIL data not found.  Default CBH will be used."
        cbh = np.repeat(2.0, len(times))
        cbh_flag = np.repeat(3, len(times))
        return cbh, cbh_flag

    default_cbh = config['cbh_default_ht']

    cbh = np.repeat(default_cbh, len(times))
    cbh_flag = np.empty((len(times)))

    for i in range(len(times)):
        epoch = times[i]
        inner_idx = np.where((epoch_times > epoch - (inner_window * 60)) & (epoch_times < epoch + (inner_window * 60)))[0]
        outer_idx = np.where((epoch_times > epoch - (outer_window * 60)) & (epoch_times < epoch + (outer_window * 60)))[0]
        
        if len(np.where(cloudHeight[inner_idx] >= 0)[0]) > 0:
            # Found a cloud in the inner window
            minCBH = np.min(cloudHeight[inner_idx])
            cbh_flag[i] = 1            
        elif len(np.where(cloudHeight[outer_idx] >= 0)[0]) > 0:
            # Found a cloud in the outer window
            minCBH = np.min(cloudHeight[outer_idx])
            cbh_flag[i] = 2
        elif len(inner_idx) == 0 and len(outer_idx) == 0:
            # There aren't any observations from in either the inner or outer window.  Using default CBH.
            minCBH = default_cbh
            cbh_flag[i] = 3
        else:
            # There are observations, but they all say no cloud.
            minCBH = default_cbh
            cbh_flag[i] = 0
            continue

        cbh[i] = minCBH
    print "Reading in VCEIL data...Done."
    return cbh, cbh_flag

def read_mwr_data(config, date, btime, etime):
    """
        read_mwr_data

        This function performs the following functions:
            a.) Reads in the microwave radiometer data for a specific MWR type and date. (if it exists)
            b.) Stores the needed fields for retrieving the thermodynamic profile.
            c.) Uses the VIP config dictionary to develop the Y vector.
            d.) Saves the brightness temperature uncertainities.

        Parameters
        ----------
        config : the config_dict returned by read VIP
        date : the date of the retrieval time in YYYYMMDD format.
        btime : a string of the time that indicates the beginning sample to be retrieved.
        etime : a string of the time that indicates the end sample to be retrieved.

        Returns
        -------
        oe_input : a dictionary containing the variables from the microwave radiometer file 
                   needed for the OE equation to work.
    """

    # Generate a list of files that have the string "date" in the file
    mwr_files = glob.glob(config['mwr_dir'] + '/*' + date + '*.cdf')
    if len(mwr_files) >= 2:
        print "Retrieval found " + str(len(mwr_files)) + " files that have " + date + \
                " in its filename within the directory: " + config['mwr_dir']
        print "MWRoe cannot isolate one MWR data file to retrieve from.  Aborting the retrieval."
        sys.exit()
    elif len(mwr_files) == 0:
        print "Unable to find any microwave radiometer files to perform the retrieval on."
        sys.exit()
    else:
        print "Retrieval found: " + mwr_files[0]
        mwr_fn = mwr_files[0]

    if config['mwr_type'] == 1:
        print "MWRoe is expecting to find a HATPRO file..."
        oe_input = read_HATPRO(mwr_fn, config, date, btime, etime)
    else:
        print "No other type of microwave radiometer files are supported right now."
        sys.exit()

    return oe_input


def read_HATPRO(mwr_fn, config, date, btime, etime):
    """
        read_HATPRO

        This function reads HATPRO microwave radiometer files containing
        both zenith and off-zenith data.  It:
            a.) Reads in the microwave radiometer data for a specific date.(if it exists)
            b.) Stores the needed fields for retrieving the thermodynamic profile.
            c.) Uses the VIP config dictionary to develop the Y vector.
            d.) Saves the brightness temperature uncertainities.

        Parameters
        ----------
        mwr_fn : the path to the HATPRO file
        config : the config_dict returned by read VIP
        date : the date of the retrieval time in YYYYMMDD format.
        btime : a string of the time that indicates the beginning sample to be retrieved.
        etime : a string of the time that indicates the end sample to be retrieved.

        Returns
        -------
        oe_input : a dictionary containing the variables from the microwave radiometer file 
                   needed for the OE equation to work.
    """
    mwr_file = Dataset(mwr_fn, 'r')

    # Load in the time fields
    #times = mwr_file.variables['time'][:]
    epoch_times = mwr_file.variables['time_since_19700101'][:]
    dts = num2date(epoch_times, 'seconds since 1970-01-01 00:00:00+00:00')

    # Find the bounds for the time frame we want to retrieve from.
    start_dt = datetime.strptime(date + btime, '%Y%m%d%H%M')
    end_dt = datetime.strptime(date + etime, '%Y%m%d%H%M')
    start_dt = date2num(start_dt, 'seconds since 1970-01-01 00:00:00+00:00')
    end_dt = date2num(end_dt, 'seconds since 1970-01-01 00:00:00+00:00')
    idx = np.where((start_dt < epoch_times) & (end_dt > epoch_times))[0]

    # Try to ensure that a sample gets retrieved.
    if len(idx) == 0:
        end_dt = start_dt + (60*15)
        print "MWRoe was unable to find a spectra within the requested time frame to retrieve on."
        print "Setting the end time to the requested start time + 15 minutes to search for a sample..."
        idx = np.where((start_dt < epoch_times) & (end_dt > epoch_times))[0]
        if len(idx) == 0:
            print "Unable to find a sample to retrieve on."
            sys.exit()
        print str(len(idx)) + " sample found."

    dts = dts[idx]
    epoch_times = epoch_times[idx]

    # Load needed variables from current HATPRO MWR file.
    rainflags = mwr_file.variables['flag'][idx]
    p_sfcs = mwr_file.variables['air_pressure'][idx]
    elevs = mwr_file.variables['hatpro_elevation_angle'][:]
    freqs = mwr_file.variables['frequencies'][:]
    tbs = mwr_file.variables['brightness_temperature'][idx,:,:]

    # Apply the correction to the MWR pressure sensor.
    p_sfcs = config['mwr_calib_pres'][1] * p_sfcs + config['mwr_calib_pres'][0]

    # Apply any corrections to the MWR altitude, longitude, or latitude fields
    if config['mwr_lat'] == -1:
        lat = mwr_file.variables['latitude'][:]
    else:
        lat = np.asarray([config['mwr_lat']])

    if config['mwr_lon'] == -1:
        lon = mwr_file.variables['longitude'][:]
    else:
        lon = np.asarray([config['mwr_lon']])

    if config['mwr_alt'] == -1:
        alt = mwr_file.variables['altitude'][:]
    else:
        alt = np.asarray(config['mwr_alt'])

    # Isolate the indices related to the zenith observations
    zenith_idx = np.where(abs(90-elevs)<0.1)[0]
    tb_zenith = tbs[:,zenith_idx,:].squeeze()
    
    if len(tb_zenith.shape) == 1:
        tb_zenith = tb_zenith[np.newaxis,:]
    
    # Load in the zenith observations
    z_freqs_idxs = []
    for z_freq in config['zenith_freqs']:
        for file_freq_idx in range(len(freqs)):
            if abs(freqs[file_freq_idx] - z_freq)<0.01:
                z_freqs_idxs.append(file_freq_idx)

    z_freqs = freqs[z_freqs_idxs]
    filtered_tb_zenith = tb_zenith[:,z_freqs_idxs]
    tb_z_uncert = config['zenith_uncert'][z_freqs_idxs]
    elevations = 90 * np.ones(len(z_freqs))

    # Load in the off-zenith observations (need to loop through the specified elevations)
    elev_indexs = []
    for vip_e in config['off_zenith_elevs']:
        for file_e in range(len(elevs)):
            if abs(vip_e - elevs[file_e]) < 0.1:
                elev_indexs.append(file_e)

    all_freqs = z_freqs

    # If there are off-zenith observations, then let's add them to the Y vector.
    if len(elev_indexs) > 0:

        # Find the indices for the off-zenith frequencies
        oz_freqs_idxs = []
        for oz_freq in config['off_zenith_freqs']:
            for file_freq_idx in range(len(freqs)):
                if abs(freqs[file_freq_idx] - oz_freq)<0.01:
                    oz_freqs_idxs.append(file_freq_idx)

        oz_freqs = freqs[oz_freqs_idxs]
        all_freqs = np.hstack((all_freqs, np.tile(oz_freqs,len(elev_indexs)).flatten()))
        elevations = np.hstack((elevations, np.repeat(elevs[elev_indexs], len(oz_freqs))))

        tb_ozenith = tbs[:, elev_indexs,:]
        tb_ozenith = tb_ozenith[:,:, np.asarray(oz_freqs_idxs)].squeeze()

        if len(tb_ozenith.shape) == 2:
            tb_ozenith = tb_ozenith[np.newaxis,:,:]

        tb_ozenith = np.reshape(tb_ozenith, (tb_ozenith.shape[0], tb_ozenith.shape[1]*tb_ozenith.shape[2]))
        tb_oz_uncert = config['zenith_uncert'][oz_freqs_idxs]

        Y = np.hstack((filtered_tb_zenith, tb_ozenith))
        tb_uncert = tb_z_uncert

        for i in elev_indexs:
            tb_uncert = np.hstack((tb_uncert, tb_oz_uncert))
    else:
        Y = filtered_tb_zenith
        oz_freqs = []
        tb_uncert = tb_z_uncert
    
    oe_input = {}
    # Save the retrieval inputs to a dictionary
    oe_input["Y"] = Y
    oe_input["p"] = p_sfcs
    oe_input["rainflags"] = 0
    oe_input["epoch_times"] = epoch_times
    oe_input["dt_times"] = dts
    oe_input["z_freqs"] = z_freqs
    oe_input["oz_freqs"] = oz_freqs
    oe_input["all_freqs"] = 0
    oe_input["all_elevs"] = 0
    indexes = np.unique(elevations, return_index=True)[1]
    oe_input["elevations_unique"] = [elevations[index] for index in sorted(indexes)]
    oe_input["elevations"] = elevations
    oe_input["lat"] = lat
    oe_input["lon"] = lon
    oe_input["alt"] = alt
    oe_input["tb_uncert"] = tb_uncert
    oe_input["rainflags"] = rainflags

    return oe_input



