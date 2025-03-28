# from astroquery.mast import Observations
# # Search using TIC ID
# target = "TIC 1201247611"
# obs_table = Observations.query_object("TIC 1201247611", radius="0s")

from astropy.io import fits
from LIGHTCURVE_UTILS.fold_utils import binning
import numpy as np
from astropy.timeseries import LombScargle
import matplotlib.pyplot as plt
plt.ion()
plt.rcParams['font.family'] = 'serif'

# P = 0.034519
P = 49.7115504/1440
t0 = 57153.6890966

f = 'gaia14aae_20s/mastDownload/TESS/tess2021175071901-s0040-0000001201247611-0211-a_fast/tess2021175071901-s0040-0000001201247611-0211-a_fast-lc.fits'
hdul = fits.open(f)
t = hdul[1].data['TIME']
y = hdul[1].data['PDCSAP_FLUX']
dy = hdul[1].data['PDCSAP_FLUX_ERR']
quality = hdul[1].data['QUALITY']
good_inds = np.nonzero( quality != 0. )
t, y, dy = t[good_inds], y[good_inds], dy[good_inds]

good_inds = np.nonzero(~np.isnan(y))
t, y, dy = t[good_inds], y[good_inds], dy[good_inds]

# frequency, power = LombScargle(t, y, dy).autopower()
# good_inds = np.nonzero( frequency > 24 )
# frequency, power = frequency[good_inds], power[good_inds]
# P = 1/frequency[ np.argmax( power ) ]


binned_lc = binning(t, y, dy, P, N=50, cycles=3)
plt.figure(figsize=(8,2))
plt.title('Sector 40')
plt.errorbar( binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
              c='k')
plt.xlabel('Time [TJD]')



lc_dir = '/data/echickle/TESS_FFI/s0058/s0058-lc/'
ticid_arr = np.load(lc_dir+'id-4-3.npy')
flux_arr = np.load(lc_dir+'lc-4-3.npy')
t1 = np.load(lc_dir+'ts-4-3.npy')
ind = np.nonzero( ticid_arr == 1201247611 )[0][0]
y1 = flux_arr[ind]
dy1 = np.ones( y1.shape )*np.std( y1 )

good_inds = np.nonzero(~np.isnan(y1))
t1, y1, dy1 = t1[good_inds], y1[good_inds], dy1[good_inds]

frequency, power = LombScargle(t1, y1, dy1).autopower()
good_inds = np.nonzero( frequency > 24 )
frequency, power = frequency[good_inds], power[good_inds]
P = 1/frequency[ np.argmax( power ) ]
binned_lc = binning(t1, y1, dy1, P, N=50, cycles=3)

plt.figure(figsize=(8,2))
plt.title('Sector 59')
plt.errorbar( binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
              c='k')
plt.xlabel('Time [TJD]')

t = np.append(t, t1)
y = np.append(y, y1)
dy = np.append(dy, dy1)



# frequency, power = LombScargle(t, y, dy).autopower()
# good_inds = np.nonzero( frequency > 24 )
# frequency, power = frequency[good_inds], power[good_inds]
# P = 1/frequency[ np.argmax( power ) ]

P = 0.01725985871313891*2

binned_lc = binning(t, y, dy, P, N=50, cycles=3)
plt.figure(figsize=(8,2))
plt.title('Combined')
plt.errorbar( binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
              c='k')
plt.xlabel('Time [TJD]')
