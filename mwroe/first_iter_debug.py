
from netCDF4 import Dataset
from pylab import *
import os
import subprocess
import math as mt
import time as tm

def makeMonoRTMCDF(sonde_file,alta,T_z,P_z,RH_z):

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

# <codecell>

def getProfilePresRH(alt, T_z, Q_z, sfc_pres):
    K_C = 273.15
    Rd = 287
    g = 9.81

    # Convert T_zp to Kelvin.
    T_z_K = T_z + K_C

    p = np.empty(T_z_K.shape[0])
    p[0] = sfc_pres
    # Compute P_zp.
    for i in np.arange(1, len(T_z_K), 1):
        avg = (T_z_K[i] + T_z_K[i-1])/2.
        p[i]= p[i-1] * np.exp((-g * (alt[i] - alt[i-1]) )/(Rd*avg))

    P_z = p
    # Compute RH_zp.
    RH_z = q2rh(Q_z ,P_z ,T_z)

    return P_z, RH_z

def q2rh(Q_z, P_z, T_z):
    K_C = 273.15

    # Set constant parameters.
    epsilon = 622.0 # empirical coeff. to obtain Q(z), dimmensionless
    T0 = 273.15 # base temperature in K
    Rv = 461.5 # specific water vapor gas constant in J/(kg*K)
    L = 2.5*(10**6) # latent heat of vaporization in J/kg
    es0 = 611 # measured saturation vapor pressure at temp. T0 in Pa

    # Make needed unit conversions.
    T_z_K = T_z + K_C # K
    P_z_pa = P_z * 100.0 # Pa

    # Compute saturation vapor pressure using the Clausius-Clapeyron equation.
    es_z = es0 * np.exp((L/Rv)*((1/T0)-(1/T_z_K)))
    print Q_z.shape, P_z_pa.shape
    # Compute the relative humidity.
    e_z = (Q_z * P_z_pa) / (epsilon + Q_z)
    RH_z = (e_z / es_z) * 100.0 # %

    return RH_z

def jacobian(freq_filenames, LWP_n, config_dict, directory_dict, elevations, F_x, X, sfc_pres, alt, cbh, cth):

    # Pre-allocate the Jacobian matrix - array Ka.
    Ka = np.zeros((len(F_x),len(X)), dtype=np.float64)
    truth = np.zeros((len(F_x),len(X)), dtype=np.float64)
    pert = np.zeros((len(F_x),len(X)), dtype=np.float64)

    # Create an array specifying what type of variable each index of X corresponds to
    # var_types is an array of variables that describes what each variable of the X vector is
    # 0 = temperature, 1 = RH, 2 = LWP
    x = np.asarray(X)
    var_types = np.zeros(x.shape[0])
    profile_size = (var_types.shape[0]-1)/2
    var_types[profile_size:profile_size*2] = 1
    var_types[-1] = 2

    # Find the indices corresponding to the maximum height of the Jacobian.
    stacked_alt = np.hstack((alt, alt, [0])) # Concatenates the height array with itself and adds a [0] element
    idx = np.where(stacked_alt/1000. <= 18)[0]

    # This for loop computes the Jacobian for each element of the X array.
    for level in idx:  # Looping through each element in the idx array

        # Take the X vector we are currently using and pull out the T & Q profiles 
        T_z = np.squeeze(x[:profile_size])
        Q_z = np.squeeze(x[profile_size:profile_size * 2])

        # Save the unperturbed value
        old_level = float(x[level])
        if var_types[level] == 0: # means that the index level corresponse to a T value in the profile
            #perturb the temperature by the value listed in the config_dict['jac_Tpert']
            pert_val = np.interp(stacked_alt[level], [0,6000], [1,3], left=1, right=3)
            x[level] = x[level] + pert_val
            denom_pert = x[level] - old_level
            LWP = LWP_n
        elif var_types[level] == 1: # means that the index level corresponds to a Q value in the profile
            #perturb the Q by multiplying it by the
            pert_val = np.interp(stacked_alt[level], [0,6000], [.95,.25], left=.95, right=.25)
            x[level] = x[level] * pert_val
            LWP = LWP_n
            denom_pert = x[level] - old_level
        elif var_types[level] == 2: # means that the index level corresponds to a LWP value in the profile
            #perturb the LWP
            denom_pert = 10
            LWP = LWP_n + 10

        # This just transforms newly perturbed X matrix object to 2 array objects
        T_zp = np.squeeze(x[:profile_size])
        Q_zp = np.squeeze(x[profile_size: profile_size*2])

        # Get the new RH and pressure profiles
        P_zp, RH_zp = getProfilePresRH(alt ,T_zp ,Q_zp ,sfc_pres)

        # Write the new profile to a sonde netCDF file
        makeMonoRTMCDF('jacobian.cdf', alt, T_zp, P_zp, RH_zp)
        sonde_file = 'jacobian.cdf'
        
        # Run the radiative transfer model
        F_xp = gen_Fx(sonde_file, freq_filenames, LWP, config_dict, directory_dict, elevations, cbh, cth)
        os.system('rm ' + sonde_file)

        # Save the unperturbed (truth) and perturbed (pert) spectra
        truth[:,level] = F_x
        pert[:, level] = F_xp

        # Save the Ka(:,level) to be equal to the Jacobian calculation. 
        Ka[:, level] = np.asarray((F_xp - F_x) / denom_pert, dtype=np.float64)

        # Return the unperturbed value to the X vector.
        x[level] = float(old_level)

    return Ka

# <codecell>

def gen_Fx(sonde_file, freq_filenames, LWP_n, config_dict, directory_dict, elevations, cbh, cth):
    # Set the temp. sonde data file location.
    #sonde_file = config_dict['asf_dir'] + "/" + directory_dict['asf_fname']

    # Set MonoRTM adaptable parameters.
    pwv_sf = "1.0" # PWV scale factor is 1
    LWP = str(LWP_n) # g/m^2
    cbh = str(cbh) # km
    cth = str(cth) # km

    Fx = None
    elevs = None
    freqs = None

    for i in range(len(elevations)): # Loop through all the elevations used in the retrieval
        elev = elevations[i]
        for f in freq_filenames: # Loop through the 2 files (.zen or .ozen)
            if elev == 90.0 and '.zen' in f: # If the elevation is zenith and the file is .zen
                # Set the monortm_freqs environment variable
                os.system ("chmod 744" + " " + f)
                os.environ['monortm_freqs'] = f
            elif elev != 90.0 and '.ozen' in f: # If it's an off zenith and the file is .ozen
                # Set the monortm_freqs env. variable
                os.system ("chmod 744" + " " + f)
                os.environ['monortm_freqs'] = f

            # Set the monortm_freqs environment variable (Dave's suggestion)
            #exp = ' '.join(['export', 'monortm_freqs=' + f])
            #os.system(exp)

            # Build the string to be used to run MonoRTM on the command line 
            monoStr = ['monortm_v4', sonde_file, pwv_sf, LWP, cbh, cth, str(90 - elev)]

            # Run the string and get the output
            output = subprocess.Popen(monoStr, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]

            # Try to parse out the brightness temperatures from the MonoRTM output
            try:
                output = output.split('ODliq')[1].strip().split()
            except:
                # Means the MonoRTM didn't run properly
                print "ERROR with MonoRTM", output
                sys.exit()

            # The MonoRTM ran successfully, parse out the brightness temperatures
            len_out = len(output)
            output = np.asarray(output, dtype=np.float64)
            output = np.reshape(output, (len_out/7, 7))
            monoRTM_output = output[:,1] # The brightness temperatures

            if Fx == None:
                # If the Fx variable hasn't been assigned to anything yet
                # then we need to assign it.  This is the zenith Tbs.
                Fx = monoRTM_output
            else:
                # If the Fx variable already is assigned to something, we need to append
                # This will happen if there are additional Tbs (off-zenith).
                Fx = np.hstack((Fx, monoRTM_output))

    return Fx

# Load in the observations into the Y vector:
data = Dataset('../../../data/theoretical_data_ddt/xxxhatpromwrX1.x1.20071001.003000.cdf')
index = 0
Y = np.matrix(data.variables['brightness_temperature'][index,0,:14]).T
freqs = data.variables['frequencies'][:]
sfc_pres = 963.5

# Load in the prior (Xa, Sa):
prior = Dataset('../../../prior/Xa_Sa_datafile.55_levels.all_months.fkb.cdf')
alta = prior.variables['height'][:]*1000.
x_a = prior.variables['mean_prior'][:]
x_a = np.concatenate((x_a, [0]))

# Make the Se matrix
Sa_pre = prior.variables['covariance_prior'][:]
Sa = np.zeros((111,111))
Sa[:110,:110] = Sa_pre
Sa[-1,-1] = 100**2
Sa = np.matrix(Sa)

# Make the Se Matrix
Se = np.empty((14,))
Se[:7] = .5
Se[7:11] = 1.
Se[11:] = .6
Se = np.power(np.matrix(np.diag(Se)),2)


# Set the x_c equal to the prior for the first iteration
x_c = x_a
x_a = np.matrix(x_a).T
for i in range(10):
    
    # Make the prior netCDF file for the forward model calculation    
    T_z = np.squeeze(x_c[:55])
    Q_z = np.squeeze(x_c[55:55+55])
    p, RH = getProfilePresRH(alta, T_z, Q_z, sfc_pres)
    makeMonoRTMCDF('prior_comp.cdf', alta, T_z, p, RH)
   
    # Make the F_x calculation
    config_dict = 'dummy'
    directory_dict = 'dummy'
    monortm_freqs_files = ['tmp/mwr_retr_monortm_freqs27251.zen']
    LWP_n = 0
    sonde_file = 'prior_comp.cdf'
    os.environ['monortm_config'] = 'tmp/mwr_retr_monortm_config_66L.txt'
    F_x = gen_Fx(sonde_file, monortm_freqs_files, LWP_n, config_dict, directory_dict, [90], 2, 2.3)

    # Compute the Jacobian K:
    K = jacobian(monortm_freqs_files, LWP_n, config_dict, directory_dict, [90], F_x, x_c, sfc_pres, alta, 2, 2.3)
    K = np.matrix(K)

    # Fix the dimensions of the numpy matrices so the OE equation doesn't go fubar
    F_x = np.matrix(F_x).T
    x_c = np.matrix(x_c).T
    
    # Solve the OE Equation:
    gamma = 1
    inv_Sop = (gamma*np.linalg.inv(Sa)) + (K.T * np.linalg.inv(Se) * K)
    Sop = np.linalg.inv(inv_Sop)
    A = Sop * K.T * np.linalg.inv(Se) * K
    x = x_a + ((Sop * K.T * np.linalg.inv(Se)) * ((Y - F_x) + (K * (x_c - x_a))))
    print "LWP: " + str(np.asarray(x)[-1][0]) + ' +/-' + str(np.sqrt(np.asarray(Sop[-1,-1])))
    conv_norm = np.asarray((x - x_c).T * inv_Sop * (x - x_c))[0] 
    stop
    if conv_norm[0] < ( len(x) / 10. ): 
        print "CONVERGED!"
        break
    else:
        x_c = x
        LWP_n = np.asarray(x)[-1][0]
        print LWP_n
        if LWP_n < 0:
            LWP_n = 0
        print "NEED TO ITERATE AGAIN."
        print np.asarray(x_c[:,0]).squeeze().shape
        x_c =  np.asarray(x_c[:,0]).squeeze()
