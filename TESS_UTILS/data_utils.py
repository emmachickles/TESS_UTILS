import numpy as np

# def phase_fold(t, P, t0=0, mid_eclipse=False):
#     t = t - t0
#     phases = (t % P) / P

#     if mid_eclipse:
#         inds = np.nonzero( phases > 0.5 )
#         phases[inds] -= 1

#     return phases

# def bin(t, y, P, t0, nbins=200, cycles=3, mid_eclipse=False):

#     phases = phase_fold(t, P, t0)
#     bin_size = 1/nbins
#     binned_lightcurve = []

#     for i in range(bin_size):
#         phase_bin = bin_size*i        
#         phase_start = phase_bin - bin_size/2
#         phase_end = phase_bin + bin_size/2

#         inds_bin = np.nonzero( (phases>phase_start) * (phases<phase_end) )
        
#         binned_lightcurve.append( np.nanmean( y[inds_bin] ) )

        

