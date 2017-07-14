
## radprofile.py - A library of classes and functions for dealing with radprofiles
##
## Initiated 2017.07.14 - lia@astro.wisc.edu
##--------------------------------------------------------------------------------

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

    def plot(self, ax, **kwargs):
        xerr = 0.5 * (self.bin_hi - self.bin_lo)
        ax.errorbar(self.bin_mid, self.val, xerr=xerr, yerr=self.val_err, ls='', **kwargs)
        ax.set_xlabel('Radius %s' % self.bin_unit)
        ax.set_ylabel(self.val_unit)
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
