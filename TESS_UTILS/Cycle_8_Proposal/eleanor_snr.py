import numpy as np

# import eleanor
# star = eleanor.Source(tic=942081746)
# data = eleanor.TargetData(star)
# q = data.quality == 0
# t = data.time[q]
# y = data.corr_flux[q]

# y = y / np.nanmedian(y)
# dy = np.ones(len(y)) * np.nanstd(y)
# np.savetxt('eleanor_TIC942081746.txt', np.array([t,y,dy]).T)

from TESS_UTILS.lightCurve import clean
t, y, dy = np.loadtxt('eleanor_TIC942081746.txt').T
t, y, dy = clean(t, y, dy)

from TESS_UTILS.periodSearch import BLS
period, q, phi0, freqs, power = BLS(t, y, dy, pmin=400/86400)

print(np.max(power) / np.median(np.abs(power-np.median(power))))

from LIGHTCURVE_UTILS.fold_utils import binning

period = 0.043939651502263144

binned_lc = binning(t, y, dy, period, N=25, cycles=2, t0 = 3014.56006)

import matplotlib.pyplot as plt

plt.ion()
plt.rcParams['font.family'] = 'serif'
fig, ax = plt.subplots(figsize=(4,2.5), ncols=2, nrows=2)
ax[0][0].errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
                          c='k')
inds = np.nonzero( (1/freqs)< 126/1440 )
ax[0][1].plot(1/freqs[inds]*1440, power[inds], '-k', lw=1)


plt.figure()
plt.plot(t, y, '.k')

plt.figure()
plt.errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='')

plt.figure()
plt.plot(freqs, power, '-k', lw=1)

import csv
topcat_file = '/home/echickle/data/cycle5_gaia.txt'
with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break
import pandas as pd
df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)

ticid = 942081746
row = np.where(df['TICID'] == ticid)[0][0]

sector = int(df['SECTOR'][row])
cam = int(df['CAM'][row])
ccd = int(df['CCD'][row])
ticid = int(df['TICID'][row])
print(df['PEAK_SIG'][row])

lc_dir = '/data/echickle/TESS_FFI/s00{}/s00{}-lc/'.format(sector, sector)
ticid_arr = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
flux_arr = np.load(lc_dir+'lc-{}-{}.npy'.format(cam, ccd))
t = np.load(lc_dir+'ts-{}-{}.npy'.format(cam, ccd))

ind = np.nonzero( ticid_arr.astype('int') == ticid )[0][0]
y = flux_arr[ind]
dy = np.ones( y.shape ) * np.std( y )


period, q, phi0, freqs, power = BLS(t, y, dy, pmin=400/86400)
print(np.max(power) / np.median(np.abs(power-np.median(power))))

binned_lc = binning(t, y, dy, period, N=25, cycles=2, t0= 3014.56006)
binned_lc[:,2] /= np.abs(np.median(binned_lc[:,1]))
binned_lc[:,1] /= np.abs(np.median(binned_lc[:,1]))
if np.median(binned_lc[:,1]) < 0:
    binned_lc[:,1]+=2

plt.figure()
plt.plot(t, y, '.k')

plt.figure()
plt.errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='')

plt.figure()
plt.plot(freqs, power, '-k', lw=1)

ax[1][0].errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
                          c='k')
inds = np.nonzero( (1/freqs)< 126/1440 )
ax[1][1].plot(1/freqs[inds]*1440, power[inds], '-k', lw=1)

# ax[0][0].set_xticks(
ax[0][0].set_xlim( ax[1][0].get_xlim() )
ax[0][1].set_xlim( ax[1][1].get_xlim() )
ax[0][0].set_xticklabels([])
ax[0][1].set_xticklabels([])
# ax[0][0].set_ylim( ax[1][0].get_ylim() )
ax[0][1].set_ylim( ax[1][1].get_ylim() )
ax[1][0].set_xlabel('Orbital Phase')
ax[1][1].set_xlabel('Period [min]')
ax[0][1].yaxis.tick_right()
ax[1][1].yaxis.tick_right()
fig.text(0.02, 0.5, 'Relative Flux', va='center', rotation='vertical', ha='left')
fig.text(0.98, 0.5, 'BLS Power', va='center', rotation=270, ha='right')
fig.subplots_adjust(wspace=0, hspace=0, bottom=0.16, left=0.16, top=0.99, right=0.82)
out_dir = '/home/echickle/out/250324/'

fig.savefig(out_dir+'eleanor.pdf')
