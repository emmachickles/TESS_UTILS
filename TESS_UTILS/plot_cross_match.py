import numpy as np
import matplotlib.pyplot as plt

res_match = np.loadtxt('/pool001/echickle/match.txt')
atlas_dir = '/pool001/echickle/ATLAS/'
tess_dir = '/nobackup1c/users/echickle/'
out_dir = '/pool001/echickle/plots/'


def phase_fold(t, P, t0, centr_0=False):
    t = t - t0
    phases = (t%P)/P
    if centr_0:
        inds = np.nonzero(phases > 0.5)
        phases[inds] += -1
    return phases

def binning(t,y,dy,P,t0=0,N=500,cycles=3):

    # remove nans
    inds = np.nonzero(~np.isnan(y))
    t, y, dy = t[inds], y[inds], dy[inds]
    
    binned_LC=[]
    mean_phases=np.linspace(0,1-1/N,N)
    phases=phase_fold(t, P, t0)
    lightcurve=np.array((phases,y,dy)).T
    lightcurve=lightcurve[np.argsort(lightcurve[:,0])]
    for i in mean_phases:
        lightcurve_bin=lightcurve[lightcurve[:,0]>i]
        lightcurve_bin=lightcurve_bin[lightcurve_bin[:,0]<i+1/N]
        weights=1/(lightcurve_bin[:,2]**2)
        weighted_mean_flux=np.sum(lightcurve_bin[:,1]*weights)/np.sum(weights)
        weighted_mean_flux_error=np.sqrt(1/np.sum(weights))
        binned_LC.append((i+0.5/N,weighted_mean_flux,weighted_mean_flux_error))
    binned_LC=np.array(binned_LC)
    binned_LC[:,2]=binned_LC[:,2]/np.nanmedian(binned_LC[:,1])
    binned_LC[:,1]=binned_LC[:,1]/np.nanmedian(binned_LC[:,1])

    if cycles==1:
        binned_LC=binned_LC
    elif cycles==2:
        binned_LC2=np.array((binned_LC[:,0]-1,binned_LC[:,1],binned_LC[:,2])).T
        binned_LC=np.vstack((binned_LC2,binned_LC))
    elif cycles==3:
        binned_LC2=np.array((binned_LC[:,0]-1,binned_LC[:,1],binned_LC[:,2])).T
        binned_LC3=np.array((binned_LC[:,0]+1,binned_LC[:,1],binned_LC[:,2])).T
        binned_LC=np.vstack((binned_LC2,binned_LC))    
        binned_LC=np.vstack((binned_LC,binned_LC3))  	

    return binned_LC


for i in range(len(res_match)):

    ticid = np.int64(res_match[i][0])
    gaiaid = np.int64(res_match[i][1])
    period = res_match[i][2]
    period_a = res_match[i][3]
    sector = np.int64(res_match[i][4])
    cam = np.int64(res_match[i][5])
    ccd = np.int64(res_match[i][6]) 

    t, y, dy=np.loadtxt(atlas_dir+str(gaiaid),usecols=(0,3,4),skiprows=0).T
    binned_atlas = binning(t, y, dy=1, P=period_a, t0=0, N=500, cycles=3)

    t = np.load(tess_dir+'s00'+str(sector)+'-lc/ts-'+str(cam)+'-'+str(ccd)+'.npy')
    ccd_tic = np.load(tess_dir+'s00'+str(sector)+'-lc/id-'+str(cam)+'-'+str(ccd)+'.npy')
    ind = np.nonzero(ccd_tic == ticid)[0][0]
    y = np.load(tess_dir+'s00'+str(sector)+'-lc/lc-'+str(cam)+'-'+str(ccd)+'.npy')
    y = y[ind]

    dy = np.ones(y)*np.std(y)

    binned_tess = binning(t, y, dy=dy, P=period, t0=0, N=500, cycles=3)

    fig, ax = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    ax[0].errorbar(binned_atlas[:,0], binned_atlas[:,1], yerr=binned_atlas[:,2], fmt='o', color='blue', label='ATLAS')
    ax[0].set_ylabel('Relative Flux')
    ax[0].set_title('ATLAS and TESS Phase Folded Light Curves')
    ax[0].legend()
    ax[1].errorbar(binned_tess[:,0], binned_tess[:,1], yerr=binned_tess[:,2], fmt='o', color='red', label='TESS')
    ax[1].set_xlabel('Phase')
    ax[1].set_ylabel('Relative Flux')
    ax[1].legend()
    plt.suptitle('GAIA ID: {} TIC ID: {} Period: {:.2f} Sector: {} Cam: {} Ccd: {}'.format(gaiaid, ticid, period, sector, cam, ccd))
    plt.savefig(out_dir+'phase_folded_light_curve_{}_{}_{}.png'.format(gaiaid, ticid, period))


