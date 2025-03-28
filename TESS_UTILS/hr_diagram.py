import numpy as np
import matplotlib.pyplot as plt
plt.ion()
from astropy.io import fits


wd_main='/data/GaiaEDR3_WD_main.fits'
rp_ext='/data/GaiaEDR3_WD_RPM_ext.fits'

source_id = np.empty(0)
bp_rp = np.empty(0)
parallax = np.empty(0)
gmag = np.empty(0)

maincat = fits.open(wd_main)
source_id = np.append(source_id, maincat[1].data['source_id'])
bp_rp = np.append(bp_rp, maincat[1].data['bp_rp'])
parallax = np.append(parallax, maincat[1].data['parallax'])
gmag = np.append(gmag, maincat[1].data['phot_g_mean_mag'])

rpmext = fits.open(rp_ext) 
source_id = np.append(source_id, rpmext[1].data['source_id'])
bp_rp = np.append(bp_rp, rpmext[1].data['bp_rp'])
parallax = np.append(parallax, rpmext[1].data['parallax'])
gmag = np.append(gmag, rpmext[1].data['phot_g_mean_mag'])

abs_mag = gmag+5*(np.log10(parallax)-2)

# fig, ax = plt.subplots(figsize=(4,2))
# # _ = ax.hist2d(bp_rp, abs_mag, bins=200, range=[[-0.6, 1.9], [4, 15.5]], density=True,
# #               cmin=0.03)
# ax.plot(bp_rp, abs_mag, '.k', ms=1, alpha=0.2)
# ax.set_xlabel('Gaia BP-RP')
# ax.set_ylabel('Absolute Magnitude (Gaia G)')
# ax.invert_yaxis()

import os
import pandas as pd
wd_tab='/home/echickle/data/WDs.txt'
wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
wd_tic = wd_cat[0].to_numpy()
wd_ra = wd_cat[1].to_numpy()
wd_dec = wd_cat[2].to_numpy()
wd_gid = wd_cat[3].to_numpy()


data_dir = '/data/echickle/TESS_FFI/'
for sector in range(56, 70):
    print('Sector {}'.format(sector))
    
    lc_dir = data_dir+'s00{}/s00{}-lc/'.format(sector, sector)
    out_dir = data_dir+'s00{}/s00{}-hr/'.format(sector, sector)

    os.makedirs(out_dir, exist_ok=True)

    fnames = os.listdir(lc_dir)
    fnames = [f for f in fnames if 'id-' in f]
    fnames.sort()

    if sector == 56:
        fnames = fnames[4:]
    
    for f in fnames:
        print(f[3:-4])
        lc_tic = np.load(lc_dir+f)

        bprp_arr = []
        absmag_arr = []
        per_arr = []
        
        for ticid in lc_tic:
            ind = np.nonzero( wd_tic == ticid )[0]
            gaiaid = wd_gid[ind]
            ra = np.float64(wd_ra[ind])
            dec = np.float64(wd_dec[ind])

            ind = np.nonzero(source_id == np.float64(gaiaid))[0]
            if len(ind) > 0:
                bprp_arr.append(bp_rp[ind][0])
                absmag_arr.append(abs_mag[ind][0])
            else:
                import astropy.units as u
                from astropy.coordinates import SkyCoord
                from astroquery.gaia import Gaia
                Gaia.ROW_LIMIT = 5
                Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
                coord = SkyCoord(ra=ra, dec=dec,
                                 unit=(u.degree, u.degree), frame='icrs')
                width, height = u.Quantity(2, u.arcsec), u.Quantity(2, u.arcsec)
                r = Gaia.query_object_async(coordinate=coord, width=width, height=height)

                # j = Gaia.cone_search_async(coord, radius=u.Quantity(3, u.arcsec))
                if len(r['phot_g_mean_mag']) > 0 and len(r['bp_rp']) > 0:

                # if len(j.get_results()['phot_g_mean_mag']) > 0 and \
                #    len(j.get_results()['bp_rp']) > 0:
                    # c_targ = j.get_results()['bp_rp'][0]        
                    # g_targ = j.get_results()['phot_g_mean_mag'][0]
                    # p_targ = j.get_results()['parallax'][0]
                    c_targ = r['bp_rp'][0]
                    g_targ = r['phot_g_mean_mag'][0]
                    p_targ = r['parallax'][0]
                    if str(p_targ) == '--':
                        m_targ = np.nan                        
                    else:
                        m_targ = g_targ + 5*(np.log10(p_targ)-2)
                    bprp_arr.append(c_targ)
                    absmag_arr.append(m_targ)
                        
                else:
                    bprp_arr.append(np.nan)
                    absmag_arr.append(np.nan)                    
                

        bprp_arr = np.array(bprp_arr)
        absmag_arr = np.array(absmag_arr)

        np.savetxt(out_dir+'hr-'+f[3:-4]+'.txt', np.array([bprp_arr, absmag_arr]))
        
            

        
