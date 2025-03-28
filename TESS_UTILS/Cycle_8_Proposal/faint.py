import numpy as np
import csv
import matplotlib.pyplot as plt
from TESS_UTILS.periodSearch import load_BLS_results
from LIGHTCURVE_UTILS.fold_utils import binning

plt.ion()
plt.rcParams['font.family'] = 'serif'


bls_results_dir = '/data/echickle/BLS_Results/ZTF_Cross_Match/'
out_dir = '/home/echickle/out/250324/'

# TICID      RA               DEC               G
# 806302837  123.672490973881 -64.4478883060073 20.1591 61 4 1
# 2025680944 336.990200348119 15.0492521992705  20.1247 56 1 2
# ticid_list = [806302837]

# # Cross match this with TOPCAT
# bls_results = load_BLS_results(bls_results_dir)
# with open(out_dir+'ZTF_CM.txt', 'w', newline='') as f:
#     f.write('# ' + '\t'.join(bls_results.keys()) + '\n')    
#     writer = csv.writer(f, delimiter='\t')
#     for row in zip(*bls_results.values()):
#         writer.writerow(row)

import pandas as pd
topcat_file = '/home/echickle/data/ZTF_CM_Gaia.txt'

with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break

df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)

plt.figure()
plt.plot( df['Gmag'][df['RECOVERED']].to_numpy(), df['PERIOD'][df['RECOVERED']].to_numpy(), '.k' )

sorted_inds = np.argsort( df['Gmag'][df['RECOVERED']] )
row_numbers = df['Gmag'][df['RECOVERED']].index
row_numbers = row_numbers[sorted_inds]
# row_numbers = row_numbers[-24:-16]
# row_numbers = row_numbers[-16:-8]
# row_numbers = row_numbers[-8:]
row_numbers = [221, 58, 13, 377, 75, 244, 378, 337] # 3 not in WD catalog
row_numbers = [377, 214, 196, 9, 287, 75, 54, 244] # all not in WD catalog

fig, ax = plt.subplots(figsize=(8,3), ncols=4, nrows=2)

for i, row in enumerate(row_numbers):

    axrow = i // 4
    axcol = i - (i//4 * 4)
    
    sector = int(df['SECTOR'][row])
    cam = int(df['CAM'][row])
    ccd = int(df['CCD'][row])
    ticid = int(df['TICID'][row])
    period = df['col4'][row]
    t0 = df['EPO'][row]
    gmag = df['Gmag'][row]
    print(ticid)
    print(df['PEAK_SIG'][row])
    print(df['col1'][row])
    print()

    if i in [1, 3, 4]:
        t0 += period/2

    lc_dir = '/data/echickle/TESS_FFI/s00{}/s00{}-lc-ZTF/'.format(sector, sector)
    ticid_arr = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
    flux_arr = np.load(lc_dir+'lc-{}-{}.npy'.format(cam, ccd))
    time = np.load(lc_dir+'ts-{}-{}.npy'.format(cam, ccd))

    ind = np.nonzero( ticid_arr == ticid )[0][0]
    flux = flux_arr[ind]
    fluxe = np.ones( flux.shape ) * np.std( flux )    
    binned_lc = binning(time, flux, fluxe, period, t0=t0, N=20, cycles=2)
    binned_lc[:,2] /= np.abs(np.median(binned_lc[:,1]))
    binned_lc[:,1] /= np.abs(np.median(binned_lc[:,1]))
    if np.median(binned_lc[:,1]) < 0:
        binned_lc[:,1]+=2

    from scipy.interpolate import make_smoothing_spline
    spl = make_smoothing_spline(binned_lc[:,0], binned_lc[:,1], lam=0.0001)
    xnew = np.linspace(-1, 1, 1000)
    ax[axrow][axcol].plot(xnew, spl(xnew), '-k')
        
    ax[axrow][axcol].errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='',
                              c='gray')    
    
    if axrow == 0:
        # ax[axrow][axcol].set_ylim([0.942, 1.07])
        ax[axrow][axcol].set_ylim([0.89, 1.1])        
        if axcol == 0:
            # ax[axrow][axcol].set_yticks([0.95, 1, 1.05])
            # ax[axrow][axcol].set_yticklabels(['0.950', '1.000', '1.050'])
            ax[axrow][axcol].set_yticks([0.9, 1, 1.1])
            ax[axrow][axcol].set_yticklabels(['0.900', '1.000', '1.100'])       
        else:
            ax[axrow][axcol].set_yticks([])
            ax[axrow][axcol].set_yticklabels([])
    if axrow == 1:
        # ax[axrow][axcol].set_ylim([0.984, 1.008])
        ax[axrow][axcol].set_ylim([0.985, 1.008])        
        if axcol == 0:
            # ax[axrow][axcol].set_yticks([1, 0.985])
            # ax[axrow][axcol].set_yticklabels(['1.000', '0.985'])            
            ax[axrow][axcol].set_yticks([1, 0.99])
            ax[axrow][axcol].set_yticklabels(['1.000', '0.990'])
        else:
            ax[axrow][axcol].set_yticks([])            
            ax[axrow][axcol].set_yticklabels([])
    ax[axrow][axcol].set_xlim([-1., 1.])
    if axrow == 0:
        ax[axrow][axcol].set_xticklabels([])
    if axrow == 1:
        ax[axrow][axcol].set_xticks([-0.5, 0, 0.5])
        ax[axrow][axcol].set_xticklabels(['-0.5', '0.0', '0.5'])

    ax[axrow][axcol].text(0.75, 0.05, 'Gmag\n'+str(np.round(gmag, 2)),
                          transform=ax[axrow][axcol].transAxes,
                          ha='center', va='bottom')
    ax[axrow][axcol].text(0.25, 0.05, 'Period\n'+str(np.round(period*1440, 1))+' m',
                          transform=ax[axrow][axcol].transAxes,
                          ha='center', va='bottom')    

fig.text(0.55, 0.02, 'Orbital Phase', ha='center', va='bottom')
fig.text(0.02, 0.5, 'Relative Flux', va='center', rotation='vertical', ha='left')
# ax[0][0].set_ylabel('Relative flux')
# ax[1][0].set_ylabel('Relative flux')
plt.subplots_adjust(hspace=0, wspace=0, bottom=0.15, left=0.1, top=0.97, right=0.98)
# plt.savefig(out_dir+'faint.png', dpi=300)
plt.savefig(out_dir+'faint.pdf')
    
    



    
    
