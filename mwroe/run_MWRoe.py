import reader
import writer
import forwardmodel
import numpy as np
import os
import sys
from netCDF4 import Dataset
from pylab import *
from datetime import datetime
import helper

###################################################################
#
#       MWRoe - Microwave Radiometer optimal estimation
#       Developed by Greg Blumberg (OU/CIMMS)
#                    Dave Turner (NSSL)
#                    Stephen Castleberry (OU)
#
#       Copyrighted 2014
#
###################################################################

try:
    date = sys.argv[1]
    vip = sys.argv[2]
    prior = sys.argv[3]
    btime = sys.argv[4]
    etime = sys.argv[5]
except:
    print "Incorrect command line arguments!  Need the following:\n" + \
          "python run_MWRoe.py YYYYMMDD <path to VIP> <path to prior> HHMM HHMM"
    sys.exit()

os.system('rm -R run_mon*')
# TODO:
# - create the ability to run in Append Mode (via command line)
#   Append Mode will not take a date, but will grab the current date/time until 2359, look for the input data file
#   and then find the existing output file and append to it.

#config = reader.readVIP('../../../vips/theoretical/mwroe.vip.theoretical-VKnoOZ.txt')
#config = reader.readVIP('../../../vips/theoretical/mwroe.vip.theoretical-VK.txt')
#prior_filename = "../../../prior/Xa_Sa_datafile.55_levels.all_months.fkb.cdf"
#date = '20071001'

print "\n<----------- Running MWRoe ------------>"

# Read in the VIP
print "Using VIP file: " + vip
config = reader.readVIP(vip)
prior_filename = prior

# Wipe the working directory clean
os.system('rm ' + config['working_dir'] + '/*')

# Read in the MWR data, the VCEIL data, and the prior
oe_inputs = reader.read_mwr_data(config, date, btime, etime)
cbh, cloud_flag = reader.read_vceil_data(config, date, oe_inputs['epoch_times'])
Xa, Sa, alt, prior_info = reader.constructPrior(prior_filename, config)

# Write the MonoRTM Config and Freqs files
monortm_config_file = writer.writeMonoRTMConfig(alt/1000., config)
monortm_freqs_files = writer.writeMonoRTMFreqs(oe_inputs, config)

# Get the output filename.
out_filename = writer.constructOutputFN(oe_inputs['dt_times'], config)

if os.path.exists(out_filename):
    if config['output_clobber'] == 0:
        print "File: " + out_filename + ' already exists.\n Exiting MWRoe.'
        sys.exit()
    elif config['output_clobber'] == 1:
        print "Deleting file: " + out_filename
        os.system('rm ' + out_filename)
    elif config['output_clobber'] == 2:
        # Enter the APPEND MODE (yet to be coded)
        sys.exit()

# Make the Se Matrix
Se = np.matrix(np.diag(np.power(oe_inputs['tb_uncert'],2)))

# Set the x_c equal to the prior for the first iteration
x_c = Xa
Xa = np.matrix(Xa).T

print "Retrieving MWR data for the date: " + datetime.strftime(oe_inputs['dt_times'][0], '%m/%d/%Y') + " for the time interval of " + \
datetime.strftime(oe_inputs['dt_times'][0], '%H:%M:%S') + " UTC to " + datetime.strftime(oe_inputs['dt_times'][-1], '%H:%M:%S') \
+ " UTC."

# Loop through the samples and retrieve each one.
for samp_idx in range(len(oe_inputs['p'])):
    print "Beginning retrieval for sample: " + str(samp_idx + 1) + ' out of ' + str(len(oe_inputs['Y'])) + \
          ", observation timestamp: " + datetime.strftime(oe_inputs['dt_times'][samp_idx], '%Y-%m-%d %H:%M:%S UTC') + \
          " [CBH = " + str(np.round(cbh[samp_idx],2)) + " km]"
    # Extract the spectra, surface pressure, and cloud base height to be used in for the following retrieval
    Y = np.matrix(oe_inputs['Y'][samp_idx]).T
    sfc_pres = oe_inputs['p'][samp_idx]
    cloud_base = cbh[samp_idx]
    cloud_top = cloud_base + .3

    i = 0

    # Build arrays used to save the results from each iteration.
    conv_norms = np.zeros(config['max_iterations'])
    x_cs = np.zeros((config['max_iterations'], len(Xa)))
    Sops = np.zeros((config['max_iterations'], len(Xa), len(Xa)))
    As = np.zeros((config['max_iterations'], len(Xa), len(Xa)))
    Fxs = np.zeros((config['max_iterations'], len(Y)))
    RMSs = np.zeros(config['max_iterations'])

    conv_norms[0] = 9.9999e5

    begin_dt = datetime.today()
    converged_flag = 0
    while i < config['max_iterations']:
        # Make the Xc (or Xa if i == 0) netCDF file for the forward model calculation    
        T_z = np.squeeze(x_c[:55])
        Q_z = np.squeeze(x_c[55:55+55])
        LWP_n = x_c[-1]
        p, RH = helper.getProfilePresRH(alt, T_z, Q_z, sfc_pres)

        # Make the F_x calculation
        sonde_file = config['working_dir'] + '/prior_comp.cdf'
        writer.makeMonoRTMCDF(sonde_file, alt, T_z, p, RH)
        os.environ['monortm_config'] = monortm_config_file
        F_x = forwardmodel.gen_Fx(sonde_file, monortm_freqs_files, LWP_n, oe_inputs['elevations_unique'], cloud_base, cloud_top)

        # Compute the RMS
        RMSs[i] = helper.rms(np.asarray(Y).squeeze(), F_x)

        fmt_conv = '%.4E' % conv_norms[i]
        print "    iter is " + str(i+1) + ", di2m is " + fmt_conv + ", and RMS is " + str(np.round(RMSs[i],2))
        
        # Compute the Jacobian K:
        K = forwardmodel.jacobian(monortm_freqs_files, LWP_n, config, oe_inputs['elevations_unique'], F_x, x_c, sfc_pres, alt, cloud_base, cloud_top)
        K = np.matrix(K)

        # Fix the dimensions of the numpy matrices so the OE equation doesn't go fubar
        F_x = np.matrix(F_x).T
        x_c = np.matrix(x_c).T
        
        # Solve the OE Equation:
        gamma = 1
        inv_Sop = (gamma*np.linalg.inv(Sa)) + (K.T * np.linalg.inv(Se) * K)
        Sop = np.linalg.inv(inv_Sop)
        A = Sop * K.T * np.linalg.inv(Se) * K
        x = Xa + ((Sop * K.T * np.linalg.inv(Se)) * ((Y - F_x) + (K * (x_c - Xa))))

        # Solve for the convergence limit.
        conv_norm = np.asarray((x - x_c).T * inv_Sop * (x - x_c))[0] 

        As[i] = np.asarray(A)
        Sops[i] = np.asarray(Sop)
        x_cs[i] = np.asarray(x_c).squeeze()
        Fxs[i] = np.asarray(F_x).squeeze()

        if conv_norm[0] < ( len(x) / 10. ): # Reached the strict converged criteria.
            i = 10 # This will end the loop.
            converged_flag = 1
            # Save the profile
            x_c = x
        else:
            # Didn't meet the strict converged criteria, keep iterating, but save the profiles.
            conv_norms[i+1] = conv_norm

            # Save the profile and arrays for later.
            x_c = x
            
            # If the LWP < 0, need to set it to 0 to prevent MonoRTM from crashing.
            LWP_n = np.asarray(x)[-1][0]
            if LWP_n < 0:
                LWP_n = 0

            # Reset the x_c so it's oriented properly for the OE equation.
            x_c =  np.asarray(x_c[:,0]).squeeze()

            # Increment the index and try the retrieval again.
            i = i + 1

    converged_idx = np.where(conv_norms < (len(x)))[0]
    if converged_flag == 1:
        # Meets strict converged criteria.
        print "Converged! ( di2m << nY )"
        idx = np.where( conv_norms != 0)[0]
        iter_count = len(idx)
        idx = idx[-1]

    elif converged_flag == 0 and len(converged_idx) > 0:
        # Didn't meet the strict converged criteria, but at least one profile met the converged criteria
        idx = np.argmin(conv_norms[converged_idx])
        iter_count = len(conv_norms)
        converged_flag = 1
        print "Converged! ( di2m < nY )"

    elif converged_flag == 0:
        # Didn't converge, save the profile with the lowest RMS.
        idx = np.argmin(RMSs)
        iter_count = config['max_iterations']
        print "Max Iterations Reached -- Retrieval did not converge!  Saving iteration (iter " + str(idx) + ") with the best RMS."

    # Extract the profiles/variables from the retrieval solution
    x_c = np.asarray(x_cs[idx]).squeeze()
    T_z = x_c[:len(alt)]
    Q_z = x_c[len(alt):len(alt)*2]
    LWP_n = x_c[-1]
    if LWP_n < 0:
        LWP_n = 0

    # Perform the final RMS calculation
    p, RH = helper.getProfilePresRH(alt, T_z, Q_z, sfc_pres)
    writer.makeMonoRTMCDF(sonde_file, alt, T_z, p, RH)
    os.environ['monortm_config'] = monortm_config_file
    F_x = forwardmodel.gen_Fx(sonde_file, monortm_freqs_files, LWP_n, oe_inputs['elevations_unique'], cloud_base, cloud_top)

    # Compute additional variables about the retrieval.
    rms = helper.rms(np.asarray(Y).squeeze(), F_x)
    sic = helper.sic(Sa, Sops[idx])
    dfs = helper.dfs(As[idx])
    T_vres, Q_vres = helper.vres(As[idx], alt)

    output = {}

    # Additional calculations
    output['temperature'] = T_z
    output['waterVapor'] = Q_z
    output['LWP'] = LWP_n
    output['pressure'] = p
    output['sic'] = sic
    output['dfs'] = dfs
    output['T_vres'] = T_vres
    output['Q_vres'] = Q_vres
    output['height'] = alt / 1000. # from m to km

    # Retrieval quality flags
    output['rms'] = rms
    output['converged_flag'] = converged_flag
    output['retr_time'] = (datetime.now() - begin_dt).seconds
    output['iter_count'] = iter_count + 1

    # Matrices used in the OE equation
    output['x_c'] = x_c
    output['Akernal_op'] = As[idx]
    output['Sop'] = Sops[idx]
    output['Sa'] = Sa
    output['Xa'] = np.asarray(Xa).squeeze().T
    output['Se'] = Se
    output['F_xop'] = F_x
    output['Y'] = Y
    output['cbh'] = cloud_base
    output['cbh_flag'] = cloud_flag[samp_idx]
    output['sample_index'] = samp_idx

    # NEED TO FIX THIS PART
    writer.save_retrieval(out_filename, output, config, prior_info, oe_inputs)

    # Reset x_c = Xa for the next sample.
    x_c = np.asarray(Xa).squeeze().T

print "MWRoe finished retrieving profiles."

