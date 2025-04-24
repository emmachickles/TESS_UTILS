import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from LIGHTCURVE_UTILS.fold_utils import binning


# plt.ion()
plt.rcParams['font.family'] = 'serif'

out_dir = '/home/echickle/out/250324/'

topcat_file = '/home/echickle/data/cycle5_gaia_cut.txt'
with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break
df_cut = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)

inds = (df_cut['PEAK_SIG']>25) & (df_cut['DPHI']<0.001) & \
    (df_cut['PERIOD']<8/24) * (df_cut['PERIOD']>7/1440)

topcat_file = '/home/echickle/data/cycle5_gaia.txt'
with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break
df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)

# ticid_list = [938779482, 942081746, 860512401, 1052063842]
ticid_list = [942081746, 860512401, 938779482, 1052063842]

row_numbers = [np.where(df['TICID'] == ticid)[0][0] for ticid in ticid_list]


import matplotlib.gridspec as gridspec
fig =  plt.figure(figsize=(8,3))
gs = gridspec.GridSpec(2, 4, height_ratios=[1,1], wspace=0.0, hspace=0.6)
axs_top = [fig.add_subplot(gs[0, i]) for i in range(4)]

for i in range(4):

    row = row_numbers[i]

    sector = int(df['SECTOR'][row])
    cam = int(df['CAM'][row])
    ccd = int(df['CCD'][row])
    ticid = int(df['TICID'][row])
    period = df['PERIOD'][row]
    t0 = df['EPO'][row]
    gmag = df['phot_g_mean_mag'][row]

    N = 25
    if i == 0:
        axs_top[i].set_ylim([0.75, 1.1])
        axs_top[i].set_ylabel('Relative Flux')
    else:
        axs_top[i].set_yticks([])
        axs_top[i].set_yticklabels([])

    if i == 2:
        period *= 2
        N=50
        
    axs_top[i].set_xticks([-0.5,0,0.5])
    axs_top[i].set_xticklabels(['-0.5', '0.0', '0.5'])

    lc_dir = '/data/echickle/TESS_FFI/s00{}/s00{}-lc/'.format(sector, sector)
    ticid_arr = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
    flux_arr = np.load(lc_dir+'lc-{}-{}.npy'.format(cam, ccd))
    time = np.load(lc_dir+'ts-{}-{}.npy'.format(cam, ccd))

    ind = np.nonzero( ticid_arr.astype('int') == ticid )[0][0]
    flux = flux_arr[ind]
    fluxe = np.ones( flux.shape ) * np.std( flux )    
    binned_lc = binning(time, flux, fluxe, period, t0=t0, N=N, cycles=2)
    binned_lc[:,2] /= np.abs(np.median(binned_lc[:,1]))
    binned_lc[:,1] /= np.abs(np.median(binned_lc[:,1]))
    if np.median(binned_lc[:,1]) < 0:
        binned_lc[:,1]+=2

    from scipy.interpolate import make_smoothing_spline
    spl = make_smoothing_spline(binned_lc[:,0], binned_lc[:,1], lam=0.0001)
    xnew = np.linspace(-1, 1, 1000)
    axs_top[i].plot(xnew, spl(xnew), '-k')        

    axs_top[i].errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='', c='gray')
    # plt.close()

    axs_top[i].text(0.75, 0.05, 'Gmag\n'+str(np.round(gmag, 2)),
                          transform=axs_top[i].transAxes,
                          ha='center', va='bottom')
    axs_top[i].text(0.25, 0.05, 'Period\n'+str(np.round(period*1440, 1))+' m',
                          transform=axs_top[i].transAxes,
                          ha='center', va='bottom')        
    


ax_bottom = fig.add_subplot(gs[1, :])
ax_bottom.hist(df_cut['PERIOD'][inds]*1440., bins=200, log=True)
ax_bottom.set_xlabel('Period [min]')
ax_bottom.set_ylabel('N')

fig.text(0.5, 0.54, 'Orbital Phase', ha='center', va='bottom')
fig.subplots_adjust(hspace=0.02, bottom=0.16, top=0.99, right=0.97, left=0.08)
fig.savefig(out_dir+'period_hist.pdf')

# topcat_file = '/home/echickle/data/cycle5_gaia_cut.txt'
# with open(topcat_file) as f:
#     for line in f:
#         if line.startswith('#'):
#             header = line[1:].strip().split()
#             break
# df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)

# row_numbers = np.nonzero( (df['PERIOD'] > 0.007) * (df['PERIOD'] < 0.009) *
#                           (df['PEAK_SIG'] > 15) )[0]

row_numbers = np.nonzero( (inds * (df_cut['PERIOD']>6/1440) * (df_cut['PERIOD']<14/1440)).to_numpy() )[0]


ticid_list = [378898110]
# ticid_list = [2054836149]
row_numbers = [np.where(df_cut['TICID'] == ticid)[0][0] for ticid in ticid_list]
for i in range(len(row_numbers)):

    # row = np.random.choice(row_numbers)
    row = row_numbers[i]

    sector = int(df_cut['SECTOR'][row])
    cam = int(df_cut['CAM'][row])
    ccd = int(df_cut['CCD'][row])
    ticid = int(df_cut['TICID'][row])
    period = df_cut['PERIOD'][row]
    peak_sig = df_cut['PEAK_SIG'][row]
    peak_wid = df_cut['PEAK_WIDTH'][row]
    dphi = df_cut['DPHI'][row]
    t0 = df_cut['EPO'][row]

    print(ticid)
    print(period)
    print(peak_sig)
    print(peak_wid)
    print(dphi)
    print(sector)
    print(cam)
    print(ccd)
    print()

    lc_dir = '/data/echickle/TESS_FFI/s00{}/s00{}-lc/'.format(sector, sector)
    ticid_arr = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
    flux_arr = np.load(lc_dir+'lc-{}-{}.npy'.format(cam, ccd))
    time = np.load(lc_dir+'ts-{}-{}.npy'.format(cam, ccd))

    ind = np.nonzero( ticid_arr.astype('int') == ticid )[0][0]
    flux = flux_arr[ind]
    fluxe = np.ones( flux.shape ) * np.std( flux )    
    binned_lc = binning(time, flux, fluxe, period, t0=t0, N=50, cycles=2)
    binned_lc[:,2] /= np.abs(np.median(binned_lc[:,1]))
    binned_lc[:,1] /= np.abs(np.median(binned_lc[:,1]))
    if np.median(binned_lc[:,1]) < 0:
        binned_lc[:,1]+=2

    fig, ax = plt.subplots(figsize=(8,2))
    ax.set_title(str( np.round( period*1440, 2 ))+' min')
    ax.errorbar(binned_lc[:,0], binned_lc[:,1], binned_lc[:,2], ls='')
    fig.savefig(out_dir+str(period)+'_'+str(peak_sig)+'_'+str(ticid)+'.png')
    plt.show()
    # plt.close()
