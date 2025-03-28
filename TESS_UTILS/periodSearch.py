"""
Module: periodSearch
This module provides functionality for period searching TESS data.
"""

import numpy as np

def BLS(t, y, dy, pmin=3/1440, pmax=None, qmin=2e-2, qmax=0.12, dlogq=0.1, freqs_to_remove=None):

    import cuvarbase.bls as bls
    from astropy.timeseries import BoxLeastSquares
    import astropy.units as u

    # tmean = np.mean(t)
    # t=t-tmean
    
    # Derive baseline from the data 
    baseline = max(t) - min(t)    

    # Set up search parameters
    search_params = dict(qmin=qmin, qmax=qmax,
                         dlogq=dlogq, # logarithmic spacing of q
                         noverlap=3 # number of overlapping phase bins
                         )


    df = search_params['qmin'] / baseline
    fmax = 1/pmin     
    if pmax == None:
        fmin = 4/baseline
    else:
        fmin = 1/pmax

    nf = int(np.ceil((fmax - fmin) / df))
    freqs = fmin + df * np.arange(nf)

    power = bls.eebls_gpu_fast(t, y, dy, freqs, **search_params)

    if freqs_to_remove is not None:

        for pair in freqs_to_remove:
            idx = np.where((freqs < pair[0]) | (freqs > pair[1]))[0]
            freqs = freqs[idx]
            power = power[idx]

    imax = np.argmax(power)
    freq_best = freqs[imax]
    period=1.0/freq_best

    # Get best duty cycle (q), epoch (phi0)
    durations = np.geomspace(qmin, qmax, 100) * period * u.day
    model = BoxLeastSquares(t*u.day, y)
    power_q = model.power(np.array([period])*u.day, durations)
    imax = np.argmax(power_q.power)
    q = power_q.duration[imax].value / period
    phi0 = (power_q.transit_time[imax].value % period) / period # mid-transit
    
    return period, q, phi0, freqs, power


def load_BLS_results(mydir):

    import os

    data = np.empty((0, 17))
    fnames = os.listdir(mydir)
    fnames.sort()
    for f in fnames:
        sector = int(f[4:6])
        cam = int(f[7])
        ccd = int(f[9])
        
        tmp_data = np.loadtxt(mydir+f, delimiter=',')

        if len(tmp_data.shape) == 1:
            tmp_data = np.expand_dims(tmp_data, 0)
        
        add_data = np.array([ [sector, cam, ccd] ])
        add_data = np.repeat(add_data, len(tmp_data), axis=0)
        tmp_data = np.append(tmp_data, add_data, axis=1)
        
        data = np.append(data, tmp_data, axis=0)
        
    cols = ['TICID', 'RA', 'DEC', 'PEAK_SIG', 'ECLP_SNR', 'PEAK_WIDTH',
            'PERIOD', 'PERIOD_MIN', 'Q', 'PHI0', 'EPO', 'R_P', 'N_TRANSIT',
            'DPHI', 'SECTOR', 'CAM', 'CCD']

    bls_results = {}
    for i in range(len(cols)):
        bls_results[cols[i]] = data[:,i]
        
    return bls_results
    
    
