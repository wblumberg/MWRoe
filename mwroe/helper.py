import math
import numpy as np
import writer
import forwardmodel
from datetime import datetime, timedelta

def average_spectra(config, oe_inputs):
    tres = config['tres']
    td = tres*30
    epoch_times = oe_inputs['epoch_times']
    center = epoch_times[0] + td
    Ys = []
    Ps = []
    Dts = []
    Ep = []
    Rain = []
    while center <= epoch_times[-1]:
        idx = np.where((epoch_times >= center - td) & (epoch_times <= center + td))[0]
        new_Y = np.mean(oe_inputs['Y'][idx,:], axis=0)
        new_p = np.mean(oe_inputs['p'][idx])
        Ep.append(center)
        Dts.append(datetime.utcfromtimestamp(center))
        Ps.append(new_p)
        Ys.append(new_Y)
        if np.sum(oe_inputs['rainflags']) > 0:
            Rain.append(1)
        else:
            Rain.append(0)
        center = center + 2*td
    print "Spectral averaging reduced the number of spectra from: " + str(len(oe_inputs['Y'])) + ' to ' + str(len(Ys))
    oe_inputs['dt_times'] = np.asarray(Dts)
    oe_inputs['Y'] = np.asarray(Ys)
    oe_inputs['rainflags'] = np.asarray(Rain)
    oe_inputs['p'] = np.asarray(Ps)
    oe_inputs['epoch_times'] = np.asarray(Ep)    
    return oe_inputs

def getProfilePresRH(alt, T_z, Q_z, sfc_pres):
    '''
        getProfilePresRH

        This function converts surface pressure value and water vapor 
        mixing ratio profile into a profile of relative humidity
        and pressure.  It uses the q2rh function and the hypsometric
        equation (the virtual temperature correction is not applied here.)

        Parameters
        ----------
        alt : an array containing the height grid [meters]
        T_z : an array containing the temperature profile [C]
        Q_z : an array containing the water vapor mixing ratio profile [g/kg]
        sfc_pres : the surface pressure value [mb]

        Returns
        -------
        P_z : an array of the pressure profile [mb]
        RH_z : an array of the relative humidity profile [%]
    '''

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
        p[i]= p[i-1] * np.exp( (-g * (alt[i] - alt[i-1])) / (Rd*avg) ) 

    P_z = p
    # Compute RH_zp.
    RH_z = q2rh(Q_z ,P_z ,T_z)

    return P_z, RH_z

def q2rh(Q_z, P_z, T_z):
    '''
        q2rh

        This helper function converts a profile of water vapor mixing ratio
        into a profile of relative humidity.

        Parameters
        ----------
        Q_z : an array of water vapor mixing ratio values [g/kg]
        P_z : an array of the pressure profile [mb]
        T_z : an array of the temperature profile [C]

        Returns
        -------
        RH_z : an array of the relative humidity profile [%]
    '''
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

    # Compute the relative humidity.
    e_z = (Q_z * P_z_pa) / (epsilon + Q_z)
    RH_z = (e_z / es_z) * 100.0 # %

    return RH_z


def rms(Y, Fx):
    """
        rms

        Calculates the root mean squared error.

        Parameters
        ----------
        Y : observed spectra used in the retrieval (K)
        Fx : forward model calculation from the state vector (K)

        Returns
        -------
        rmse : the root mean squared error (K)
    """
    # Calculate the RMSe using Y and F_xop.
    rmse = np.sqrt(np.mean(np.power(Y-Fx, 2)))

    return rmse

def sic(Sa, Sop):
    """
        sic

        Calculates the Shannon Information Content.

        Parameters
        ----------
        Sa : the prior covariance matrix
        Sop : the retrieved profile (optimal) covariance matrix

        Returns
        -------
        sic : the Shannon Information Content
    """
    mp = Sa * np.linalg.inv(Sop)
    det_mp = np.linalg.det(mp)
    sic = 0.5 * np.log(det_mp)

    return sic

def dfs(A):
    """
        dfs

        Calculates the degrees of freedom of signal for portions of the retrieved profile.

        Parameters
        ----------
        A : the averaging kernal from the retrieval

        Returns
        -------
        DFS : an array with DFS for the entire retrieval, the temperature profile,
                the water vapor mixing ratio profile, and the LWP
    """
    tot_num = A.shape[0] - 1
    num_alts = tot_num / 2.
    dfs = np.empty(4)

    # Extract the diagonal of A.
    diag_A = np.diag(A,0)

    # Compute the total dfs.
    dfs[0] = np.sum(diag_A)

    # Compute T(z) dfs.
    dfs[1] = np.sum(diag_A[0:num_alts])

    # Compute Q(z) dfs.
    dfs[2] = np.sum(diag_A[num_alts:tot_num])

    # Compute LWP dfs.
    dfs[3] = diag_A[tot_num]

    return dfs

def vres(A, alt):
    """
        vres

        Calculates the vertical resolution of the retrieved profile via the method used in 
        Hewson (2007).

        Parameters
        ----------
        A : the averaging kernal from the retrieval
        alt : the height grid (km)

        Returns
        -------
        T_vres : the vertical resolution of the temperature profile
        Q_vres : the vertical resolution of the water vapor mixing ratio profile
    """
    tot_num = A.shape[0] - 1
    num_alts = tot_num / 2.
    zres = [(alt[1] - alt[0])/2.]
    for i in np.arange(1,num_alts - 1):
        zres.append((alt[i+1]-alt[i-1])/2.)
    zres.append((alt[num_alts-1]- alt[num_alts-2])*2.)

    A_diag = np.diag(A)
    A_diag = np.ma.masked_where(A_diag == 0, A_diag)
    T_vres = zres/A_diag[:num_alts]
    Q_vres = zres/A_diag[num_alts: num_alts*2]

    return T_vres, Q_vres

def convert_time(time, flag):
    '''
        convert_time

        This function allows for string time conversions in the format of
        'hhmm' and 'h.hf' using a datetime object.

        Parameters
        ----------
        time : the datetime object to be formatted
        flag : a flag indicating the format of the string to be returned.
        
        Returns
        -------
        time_s : a string of the time in the specified format
    '''

    if flag == 'hhmm':
        t_dec,t_int = math.modf(time)
        min_s = str(int(t_dec * 60))
        hr_s = str(int(t_int))
        if len(min_s) < 2:
            min_s = '0' + min_s
        if len(hr_s) < 2:
            hr_s = '0' + hr_s
        time_s = hr_s + min_s

    elif flag == 'h.hf':
        hr_s = time[0:2]
        min_s = time[2:4]
        min_s = str(float(min_s) / 60).replace('0.','.')
        time_s = str(float(hr_s + min_s))

    return time_s
