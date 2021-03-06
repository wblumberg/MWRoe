import subprocess
import numpy as np
import writer
import os
import helper
import sys

###
#
# This module allows for interactions between MWRoe and the MonoRTM.
#
# It contains functions to generate the F(X) vector and the Jacobian (K) matrix.
#
###

def jacobian_Kb(sonde_file, LWP_n, config_dict, oe_inputs, offsets, F_x, cbh, cth, delta_tb=0, debug=False):
    '''
        Jacobian Kb

        This function computes the Jacobian that describes the sensitivity of the forward model
        to perturbations in frequency offsets and brightness temperature offsets.


    '''
   
    tb_pert = 1 # K
    f_pert = 0.01 # GHz

    # Step 1 write new zenith and off-zenith freqency files by adding the f_pert to the known freqs.
    #writer.
    freq_filenames = writer.writeMonoRTMFreqs_FM(f_pert + np.asarray(oe_inputs['z_freqs']) + np.asarray(offsets['z_freq_offsets']), f_pert + np.asarray(oe_inputs['oz_freqs']) + np.asarray(offsets['oz_freq_offsets']), config_dict)

    # Step 2 Compute the F_xp using these new frequency files.  Subtract from F_x.
    F_xp = gen_Fx(sonde_file, freq_filenames, LWP_n, oe_inputs['elevations_unique'], cbh, cth, delta_tb, debug=debug)
    J1 = np.diag((F_xp - F_x)/f_pert)
    
    # Step 3 Compute the Brightness temperature portion of the Jacobian.  This will be linear,
    # so it'll just the 1.  Use the identity matrix here.
    J2 = np.identity(len(np.asarray(F_x).squeeze()))
        
    return np.vstack((J1, J2))

def jacobian_Ka(freq_filenames, LWP_n, config_dict, elevations, F_x, X, sfc_pres, alt, cbh, cth, delta_tb=0, debug=False):
    '''
        jacobian

        This function computes the Jacobian needed by the optimal estimation equation.
        It uses the gen_Fx function to do its dirty work.

        Parameters
        ----------
        freq_filenames : a list of the paths to the MonoRTM frequency files that MonoRTM needs.
        LWP_n : a float that is the current liquid water path value [g/m^2]
        config_dict : a dictionary that contains the configuration variables
        elevations : an array that contains the elevation levels of the radiometer.
        F_x : the precomputed F(X) vector that will be the base state for the Jacobian.
        X : the state vector containing the temperature [C], water vapor mixing ratio [g/kg], and LWP [g/m^2] values
        sfc_pres : a float of the surface pressure value [mb].
        alt : an array containing the height grid [meters]
        cbh : the current cloud base height [km]
        cth : the current cloud top height [km]

        Returns
        -------
        K : the optimal estimation Jacobian
    '''
    
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
    idx = np.where(stacked_alt/1000. <= config_dict['jac_max_alt'])[0]

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
        P_zp, RH_zp = helper.getProfilePresRH(alt ,T_zp ,Q_zp ,sfc_pres)

        # Write the new profile to a sonde netCDF file
        sonde_file = config_dict['working_dir'] + '/jacobian.cdf'
        writer.makeMonoRTMCDF(sonde_file, alt, T_zp, P_zp, RH_zp)
        
        # Run the radiative transfer model
        F_xp = gen_Fx(sonde_file, freq_filenames, LWP, elevations, cbh, cth, delta_tb, debug=debug)
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

def gen_Fx(sonde_file, freq_filenames, LWP_n, elevations, cbh, cth, delta_tb=0, debug=False):
    '''
        gen_Fx

        This function runs the MonoRTM given a frequency file, a config file, a sonde file,
        and any additional parameters the MonoRTM needs.  It outputs a vector of brightness
        temperature values for different frequencies.

        Parameters
        ----------
        sonde_file : a string that contains the path to the netCDF sonde file to be run by MonoRTM
        freq_filenames : an array containing the paths to the MonoRTM frequency files needed.
        LWP_n : a float containing the liquid water path value [g/m^2]
        elevations : an array containing the elevation values to run the MonoRTM at [deg]
        cbh : a float representing the cloud base height [km]
        cth : a float representing the cloud top height [km]

        Returns
        -------
        F_x : the F_x brightness temperature vector
    '''

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
            else:
                continue

            # Build the string to be used to run MonoRTM on the command line 
            monoStr = ['monortm_v4', sonde_file, pwv_sf, LWP, cbh, cth, str(90 - elev)]

            # Run the string and get the output
            output = subprocess.Popen(monoStr, stdout=subprocess.PIPE, stderr=subprocess.STDOUT).communicate()[0]
            
            if debug is True:
                print output
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

            if Fx is None:
                # If the Fx variable hasn't been assigned to anything yet
                # then we need to assign it.  This is the zenith Tbs.
                Fx = monoRTM_output
            else:
                # If the Fx variable already is assigned to something, we need to append
                # This will happen if there are additional Tbs (off-zenith).
                Fx = np.hstack((Fx, monoRTM_output))

    return delta_tb + Fx
