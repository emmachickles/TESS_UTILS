import pandas as pd
import numpy as np
from TESS_UTILS.periodSearch import BLS
from LIGHTCURVE_UTILS.fold_utils import binning
import matplotlib.pyplot as plt
import os

ticid = 378898110
out_dir = 'TIC378898110/'
os.makedirs(out_dir, exist_ok=True)

topcat_file = '/home/echickle/data/cycle5_gaia.txt'
with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break
df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)


row = np.where(df['TICID'] == ticid)[0][0]

sector = int(df['SECTOR'][row])
cam = int(df['CAM'][row])
ccd = int(df['CCD'][row])
ticid = int(df['TICID'][row])
print(df['PEAK_SIG'][row])
print(df['RA_1'][row], df['DEC_1'][row])

lc_dir = '/data/echickle/TESS_FFI/s00{}/s00{}-lc/'.format(sector, sector)
ticid_arr = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
flux_arr = np.load(lc_dir+'lc-{}-{}.npy'.format(cam, ccd))
t = np.load(lc_dir+'ts-{}-{}.npy'.format(cam, ccd))

ind = np.nonzero( ticid_arr.astype('int') == ticid )[0][0]
y = flux_arr[ind]
dy = np.ones( y.shape ) * np.std( y )

dy /= np.abs(np.nanmedian(y))
y /= np.abs(np.nanmedian(y))
if np.nanmedian(y) < 0:
    y += 2

period, q, phi0, freqs, power = BLS(t, y, dy, pmin=400/86400)
print(np.max(power) / np.median(np.abs(power-np.median(power))))
period = df['PERIOD'][row]

binned_lc = binning(t, y, dy, period, N=100, cycles=2, t0= 3014.56006)
binned_lc[:,2] /= np.abs(np.nanmedian(binned_lc[:,1]))
binned_lc[:,1] /= np.abs(np.nanmedian(binned_lc[:,1]))
if np.nanmedian(binned_lc[:,1]) < 0:
    binned_lc[:,1]+=2

plt.figure(figsize=(10, 5))
plt.plot(t, y, '.k')
plt.xlabel('Time [BJD-2457000]')
plt.ylabel('Relative flux')
plt.savefig(out_dir+'lc-{}-{}-{}-{}.png'.format(sector, cam, ccd, ticid), dpi=300)

plt.figure(figsize=(10, 5))
plt.errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='')
plt.xlabel('Phase')
plt.ylabel('Relative flux') 
plt.savefig(out_dir+'phase-{}-{}-{}-{}.png'.format(sector, cam, ccd, ticid), dpi=300)

plt.figure(figsize=(10, 5))
plt.plot(freqs, power, '-k', lw=1)
plt.xlabel('Frequency [1/d]')
plt.ylabel('Power')
plt.savefig(out_dir+'power-{}-{}-{}-{}.png'.format(sector, cam, ccd, ticid), dpi=300)
plt.show()