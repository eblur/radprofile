
## radprofile.py - A library of classes and functions for dealing with radprofiles
##
## Initiated 2017.07.14 - lia@astro.wisc.edu
##--------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


PIX2ARCSEC = 0.5  # arcsec / pix

class RadProfile(object):
    def __init__(self, bin_lo, bin_hi, bin_unit, val, val_err, val_unit):
        self.bin_lo   = bin_lo
        self.bin_hi   = bin_hi
        self.val      = val
        self.bin_unit = bin_unit
        self.val_err  = val_err
        self.val_unit = val_unit

    @property
    def bin_mid(self):
        return 0.5 * (self.bin_lo + self.bin_hi)

    def renorm(self, y, yerr=0.0, new_unit=None):
        assert y != 0.0  # Don't ever divide by zero
        x, xerr = self.val, self.val_err
        prop_err_term = (xerr/x)**2 + (yerr/y)**2
        self.val     = x / y
        self.val_err = new_val * np.sqrt(prop_err_term)  # Works out to xerr/y when yerr=0.0
        if new_unit is not None:
            self.val_unit = new_unit

    def plot(self, ax, loglog=True, **kwargs):
        xerr = 0.5 * (self.bin_hi - self.bin_lo)
        if not all(self.val_err == 0.0):
            ax.errorbar(self.bin_mid, self.val, xerr=xerr, yerr=self.val_err, ls='', **kwargs)
        else:
            ax.errorbar(self.bin_mid, self.val, xerr=xerr, ls='', **kwargs)
        ax.set_xlabel('Radius %s' % self.bin_unit)
        ax.set_ylabel(self.val_unit)
        if loglog:
            ax.set_xscale('log', nonposx='clip')
            ax.set_yscale('log', nonposy='clip')

class FITSProfile(object):
    def __init__(self, filename):
        self.data = fits.open(filename)[1].data

    @property
    def columns(self):
        return self.data.columns

    @property
    def bin_lo(self):
        return self.data['R'][:,0] * PIX2ARCSEC # arcsec

    @property
    def bin_hi(self):
        return self.data['R'][:,1] * PIX2ARCSEC # arcsec

    @property
    def bin_mid(self):
        return 0.5 * (self.bin_lo + self.bin_hi)

    def surbri_profile(self):
        surbri      = self.data['SUR_BRI'] / PIX2ARCSEC**2
        surbri_err  = self.data['SUR_BRI_ERR'] / PIX2ARCSEC**2
        surbri_unit = self.columns['SUR_BRI'].unit.rstrip(r"pixel**2") + "arcsec**2"
        result = RadProfile(self.bin_lo, self.bin_hi, 'arcsec', surbri, surbri_err, surbri_unit)
        return result

    def flux_profile(self):
        surbri      = self.data['SUR_FLUX'] / PIX2ARCSEC**2
        surbri_err  = self.data['SUR_FLUX_ERR'] / PIX2ARCSEC**2
        surbri_unit = self.columns['SUR_FLUX'].unit.rstrip(r"pixel**2/s") + r"/s/arcsec**2"
        result = RadProfile(self.bin_lo, self.bin_hi, 'arcsec', surbri, surbri_err, surbri_unit)
        return result

    def profile_from_field(self, fname):
        # Return a radial profile for a particular column name in the fits file
        assert fname in self.columns.names
        value     = self.data[fname]

        err_name = fname + "_ERR"
        if err_name in self.columns.names:
            value_err = self.data[err_name]
        else:
            value_err = np.zeros_like(value)
        
        value_unit = self.columns[fname].unit
        result = RadProfile(self.bin_lo, self.bin_hi, 'arcsec', value, value_err, value_unit)
        return result

#------- Some mathy stuff

def diff(rp1, rp2):
    """
    Subtract two radial profiles
    Returns a radial profile with RP1 - RP2
    """
    assert all(rp1.bin_lo == rp2.bin_lo)
    assert all(rp1.bin_hi == rp2.bin_hi)
    assert rp1.bin_unit == rp2.bin_unit
    assert rp1.val_unit == rp2.val_unit
    value      = rp1.val - rp2.val
    value_err  = np.sqrt(rp1.val_err**2 + rp2.val_err**2)
    result = RadProfile(rp1.bin_lo, rp1.bin_hi, rp1.bin_unit, value, value_err, rp1.val_unit)
    return result

