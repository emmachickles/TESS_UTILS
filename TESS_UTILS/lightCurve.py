"""
Module: lightCurve
This module provides functionality for extracting and saving light curves with TESS.
"""

import numpy as np

def clean(t, y, dy=None, window_length=0.1, n_std=5):
    
    from wotan import flatten
    y = flatten(t, y, window_length=window_length, method='biweight')

    med = np.nanmedian( y )
    std = np.nanstd( y )
    good_inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    t, y = t[good_inds], y[good_inds]

    if dy is not None:
        dy = dy[good_inds]
        return t, y, dy

    return t, y
