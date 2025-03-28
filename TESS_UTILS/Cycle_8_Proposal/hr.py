import numpy as np
import csv
import matplotlib.pyplot as plt
from TESS_UTILS.periodSearch import load_BLS_results
from LIGHTCURVE_UTILS.fold_utils import binning

plt.ion()
plt.rcParams['font.family'] = 'serif'


bls_results_dir = '/data/echickle/BLS_Results/TESS_Cycle_5/'
out_dir = '/home/echickle/out/250324/'

# # Cross match this with TOPCAT
# bls_results = load_BLS_results(bls_results_dir)
# with open(out_dir+'cycle5.txt', 'w', newline='') as f:
#     f.write('# ' + '\t'.join(bls_results.keys()) + '\n')    
#     writer = csv.writer(f, delimiter='\t')
#     for row in zip(*bls_results.values()):
#         writer.writerow(row)

import pandas as pd
topcat_file = '/home/echickle/data/cycle5_gaia_cut.txt'

with open(topcat_file) as f:
    for line in f:
        if line.startswith('#'):
            header = line[1:].strip().split()
            break

df = pd.read_csv(topcat_file, comment='#', delim_whitespace=True, names=header)


color = df['bp_rp']
abs_mag = df['phot_g_mean_mag'] + 5*(np.log10(df['parallax']) - 2)
period = df['PERIOD']

print(len(period))
inds = np.nonzero( ((df['PEAK_SIG']>50) * (period > 7/1440)* (period<10) * \
                    (df['PEAK_WIDTH']>5)).to_numpy())[0]
print(len(inds))
fig, ax = plt.subplots(figsize=(4,4), nrows=2)
cmap = plt.cm.viridis  # or whatever colormap you're using
norm = plt.Normalize(vmin=period[inds].min(), vmax=period[inds].max())
# from matplotlib.colors import LogNorm
# norm = LogNorm(vmin=period[period > 0].min(), vmax=period.max())
sc = ax[0].scatter(color[inds], abs_mag[inds], c=period[inds], cmap=cmap, norm=norm, s=0.5, alpha=0.2)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax[0], label='Period [d]')
# ax[0].set_xlabel('Gaia BP-RP')
# ax[0].set_ylabel('Absolute Mag. (Gaia G)')
ax[0].set_xlabel(r'$\mathrm{G}_{\rm BP}-\mathrm{G}_{\rm RP}$')
ax[0].set_ylabel(r'$\mathrm{M}_{\rm G}$')
ax[0].set_xlim([-0.5, 1.75])
ax[0].set_ylim([15.5, 2.5])
ax[1].hist(period[inds], bins=500, log=True)
ax[1].set_xlabel('Period [d]')
ax[1].set_ylabel('N')
ax[1].set_xlim([0, 10])
plt.subplots_adjust(left=0.15, bottom=0.12, top=0.98, right=0.96, hspace=0.3)
# fig.savefig(out_dir+'hr.png', dpi=300)
fig.savefig(out_dir+'hr.pdf')


# ax2 = ax.twinx()
# M_sun_G = 4.67
# mag_min, mag_max = ax.get_ylim()
# logL_min = -0.4 * (mag_min - M_sun_G)
# logL_max = -0.4 * (mag_max - M_sun_G)
# ax2.set_ylim([logL_min, logL_max])
# ax2.set_ylabel(r'$\log_{10}(L / L_\odot)$')



# # Add top axis: log(Teff)
# ax_top = ax.twiny()
# color_ticks = np.array([0.0, 0.5, 1.0, 1.5, 2.0])  # example ticks in BP-RP
# logTeff_ticks = 3.962 - 0.305 * color_ticks + 0.064 * color_ticks**2
# ax_top.set_xlim(ax.get_xlim())
# ax_top.set_xticks(color_ticks)
# ax_top.set_xticklabels([f"{lt:.2f}" for lt in logTeff_ticks])
# ax_top.set_xlabel(r'$\log_{10}(T_{\rm eff}/\mathrm{K})$')


# plt.tight_layout()
#ax.invert_yaxis()

period_min = period*1440
inds = np.nonzero( ((period_min < 60) * (period_min>7/1440) * (df['PEAK_SIG']>50) * (df['PEAK_WIDTH']>5)).to_numpy())[0]
print(len(inds))
fig, ax = plt.subplots()
cmap = plt.cm.viridis  # or whatever colormap you're using
norm = plt.Normalize(vmin=period_min[inds].min(), vmax=period_min[inds].max())
# from matplotlib.colors import LogNorm
# norm = LogNorm(vmin=period[period > 0].min(), vmax=period.max())
sc = ax.scatter(color[inds], abs_mag[inds], c=period_min[inds], cmap=cmap, norm=norm, s=1, alpha=0.75)
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Period [m]')
ax.set_xlabel('Gaia BP-RP')
ax.set_ylabel('Absolute Magnitude (Gaia G)')
ax.set_xlim([-0.5, 1.75])
ax.set_ylim([15.5, -5])

# ax2 = ax.twinx()
# M_sun_G = 4.67
# mag_min, mag_max = ax.get_ylim()
# logL_min = -0.4 * (mag_min - M_sun_G)
# logL_max = -0.4 * (mag_max - M_sun_G)
# ax2.set_ylim([logL_min, logL_max])
# ax2.set_ylabel(r'$\log_{10}(L / L_\odot)$')

# # Add top axis: log(Teff)
# ax_top = ax.twiny()
# color_ticks = np.array([0.0, 0.5, 1.0, 1.5, 2.0])  # example ticks in BP-RP
# logTeff_ticks = 3.962 - 0.305 * color_ticks + 0.064 * color_ticks**2
# ax_top.set_xlim(ax.get_xlim())
# ax_top.set_xticks(color_ticks)
# ax_top.set_xticklabels([f"{lt:.2f}" for lt in logTeff_ticks])
# ax_top.set_xlabel(r'$\log_{10}(T_{\rm eff}/\mathrm{K})$')


#ax.invert_yaxis()
plt.tight_layout()
fig.savefig(out_dir+'hr_short.png', dpi=300)

