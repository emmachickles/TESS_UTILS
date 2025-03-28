import numpy as np
import os
import matplotlib.pyplot as plt
from LIGHTCURVE_UTILS.fold_utils import binning
from astropy.io import fits

plt.ion()
plt.rcParams['font.family'] = 'serif'

mydir = '/data/echickle/BLS_Results/TESS_Cycle_5/'
fnames = os.listdir(mydir)
fnames.sort()

# 0      1   2    3    4    5    6       7           8  9     10   11
# ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi
data = np.empty((0, 14))
for f in fnames:
    tmp_data = np.loadtxt(mydir+f, delimiter=',')
    data = np.append( data, tmp_data, axis = 0 )

plt.figure(figsize=(8,2))
plt.hist(data[:,3], bins=100, log=True)
plt.xlabel('Significance')
plt.tight_layout()

good_inds = np.nonzero( data[:,3] > 25 )
print('Cut: {} / {}'.format(len(good_inds[0]), len(data)))
plt.figure(figsize=(8,2))
# plt.figure(figsize=(16,4))
res = plt.hist(data[:,6][good_inds], bins=1000, log=True)
plt.xlabel('Period [d]')
plt.ylabel('Number of objects')
plt.tight_layout()
plt.savefig('hist_period.png', dpi=300)

good_inds = np.nonzero( (data[:,3] > 25) * (data[:,6] < 6/24) )
print('Cut: {} / {}'.format(len(good_inds[0]), len(data)))
plt.figure(figsize=(8,2))
# plt.figure(figsize=(16,4))
plt.hist(data[:,6][good_inds]*1440, bins=1000, log=True)
plt.xlabel('Period [min]')
plt.ylabel('Number of objects')
plt.tight_layout()
plt.savefig('hist_period<6hr.png', dpi=300)

lc_dir = '/data/echickle/TESS_FFI/s0061/s0061-lc/'
ticid_arr = np.load(lc_dir+'id-1-1.npy')
flux_arr = np.load(lc_dir+'lc-1-1.npy')
t = np.load(lc_dir+'ts-1-1.npy')
ind = np.nonzero( ticid_arr == str(803489769) )[0][0]
y = flux_arr[ind]

ind = np.nonzero( data[:,0] == 803489769 )[0]
P = data[:,6][ind]
t0 = data[:,10][ind]

dy = np.ones( y.shape )*np.std( y )
binned_lc = binning(t, y, dy, P, t0=t0, N=200, cycles=3)

plt.figure(figsize=(8,2))
plt.errorbar( binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
              c='k')
plt.xlabel('Time [TJD]')


# data = np.empty((0, 14))
# for f in fnames[:4]:
#     tmp_data = np.loadtxt(mydir+f, delimiter=',')
#     data = np.append( data, tmp_data, axis = 0 )

# hr_dir = '/data/echickle/TESS_FFI/s0056/s0056-hr/'
# hr_data = np.empty((2,0))
# for i in [1,2,3,4]:
#     tmp_data = np.loadtxt(hr_dir+'hr-1-{}.txt'.format(i))
#     hr_data = np.append( hr_data, tmp_data, axis=1)

# per_arr = data[:,6]
# bprp_arr = hr_data[0]
# absmag_arr = hr_data[1]
        
# fig, ax = plt.subplots(figsize=(4,2))
# sc = ax.scatter(bprp_arr, absmag_arr, c=np.log10(per_arr), cmap='viridis', s=1)
# ax.set_xlabel('Gaia BP-RP')
# ax.set_ylabel('Absolute Magnitude (Gaia G)')
# fig.colorbar(sc, ax=ax, label='Period [d]')
# ax.invert_yaxis()
# fig.savefig('hr.png', dpi=300)
    
# ind = np.nonzero(wd_cat[0].to_numpy().astype('int') == objid)[0][0]
# gid = np.int64(wd_cat.iloc[ind][3])
