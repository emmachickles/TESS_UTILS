import numpy as np
import pandas as pd 
import os, pdb

atlas_dir = '/data/echickle/BLS_Results/EOF_ATLAS/'
fnames = os.listdir(atlas_dir)
res_atlas = np.empty((0, 10)) # gaiaid, sig, snr, wid, period, period_min, q, phi0, nt, epo
for f in fnames:
    res_atlas = np.append(res_atlas, np.loadtxt(atlas_dir+f, delimiter=','), axis=0)
np.savetxt('/home/echickle/out/bls_results_atlas.txt', res_atlas)

tess_dir = '/data/echickle/BLS_Results/TESS_Cycle_5/'
fnames = os.listdir(tess_dir)
res_tess = np.empty((0, 17)) # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi, sect, cam, ccd
for f in fnames:
    sect, cam, ccd = f.split('-')[1:4]
    sect = int(sect)
    cam = int(cam) 
    ccd = int(ccd.split('.')[0])

    tmp = np.loadtxt(tess_dir+f, delimiter=',')
    tmp = np.column_stack((tmp, np.full((tmp.shape[0], 3), [sect, cam, ccd])))
    res_tess = np.append(res_tess, tmp, axis=0)
np.savetxt('/home/echickle/out/bls_results_tess.txt', res_tess)

wd_cat  = pd.read_csv('/home/echickle/data/WDs.txt', header=None, sep='\s+', dtype='float') # ticid, ra, dec gaiaid, mag
res_match = []
for i in range(len(res_tess)):
    if i % 100 == 0:
        print('Processing TESS result {} of {}'.format(i, len(res_tess)))
    ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi, sector, cam, ccd = res_tess[i]
    ticid = np.int64(ticid)
    if np.count_nonzero( np.int64(wd_cat[0]) == np.int64(ticid) ) > 0:
        idx = np.where( np.int64(wd_cat[0]) == int(ticid) )[0][0]
        gaiaid = np.int64(wd_cat.iloc[idx][3])

        if np.count_nonzero( np.int64(res_atlas[:,0]) == np.int64(gaiaid) ) > 0:
            idx = np.where( np.int64(res_atlas[:,0]) == np.int64(gaiaid) )[0][0]
            gaiaid_a, sig_a, snr_a, wid_a, period_a, period_min_a, q_a, phi0_a, nt_a, epo_a = res_atlas[idx]

            if np.abs(period_min - period_min_a) < 1 or \
                np.abs(period_min - period_min_a/2) < 1 or \
                    np.abs(period_min - period_min_a*2) < 1:
                res_match.append([ticid, gaiaid, period_min, period_min_a, sector, cam, ccd])
                print('TICID: {} GAIAID: {} period: {} period_a: {}'.format(np.int64(ticid), gaiaid, period_min, period_min_a))
        
res_match = np.array(res_match)
np.savetxt('/home/echickle/out/match.txt', res_match)
