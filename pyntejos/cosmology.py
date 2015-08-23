"""Cosmology module from Neil Crighton's old astro. Deprecated, use
astropy.cosmology instead.
"""

import sys
import numpy as np

try:
    from scipy import integrate
except ImportError:
    print ('WARNING: No Scipy found, so most cosmology functions will '
           'not work.')
from math import sqrt,log10,sin,sinh,pi

# originally by Andrew Becker; becker@astro.washington.edu, some
# changes by nhmc

# Many of these adapted from
# astro-ph/9905116

# CONSTANTS
C       = 299792.458        # speed of light in km/s (exact)
PC      = 3.08567782e16     # 1 parsec in m (from wikipedia, which
                            # cites P. Kenneth Seidelmann,
                            # Ed. (1992). Explanatory Supplement to
                            # the Astronomical Almanac. Sausalito, CA:
                            # University Science Books. p. 716 and
                            # s.v. parsec in Glossary.)

class Cosmology(object):
    """
    Instantiate with optional keyword arguments (defaults shown).

    Om = 0.27:   Omega matter
    Ol = 0.73:   Omega lambda
    w  = -1.0:   pressure/density dark energy ratio
    H0 = 73:     Present day hubble constant in km/s/Mpc

    Derived values
    
    Ok:    Curvature density  (assumes Omega_total=1)
    h:     Dimensionless Hubble parameter (H0 = 70*h km/s/Mpc)
           Used to give H0-independent distances 
    Th:    Hubble time in seconds
    Dh:    Hubble distance in metres

    Methods

    Hz(z):      Hubble constant at redshift z (km/s/Mpc)
    a(z):       scale factor at redshift z
    Tl(z):      Lookback time for redshift z (seconds)
    Dc(z):      Line of sight comoving distance at z (metres)
    Dm(z):      Transverse comoving distance at z (metres)
    Da(z):      Angular diameter distance (metres)
    Da2(z1,z2): Angular diameter distance between objects at z1 and z2
    Dl(z):      Luminosity distance (metres)
    distmod(z): Distance modulus
    X(z):       Absorption distance corresponding to redshift z

    Constants
    
    C  = 2.9979E5:        Speed of light (km/s)
    PC = 3.08567782E16:   Parsec (metres)


    
    Note that the integration routine used, scipy.integrate.quad, will
    start giving nonsense for very large redshifts (>100?).

    Note also the energy density from radiation, Omega_r, is ignored
    (valid for redshifts < ~10).

    Examples
    --------
    
    """
    # Note all the distance units come from Dh, so this returns
    # everything in metres

    def __init__(self, Om=0.27, Ol=0.73, w=-1., H0=73.):
        # The defaults above imply Ok = 0 if we assume Omega total = 1
        # (which we do here).

        # all densities are in units of the critical density
        self.Om = Om                   # matter density
        self.Ol = Ol                   # lambda density
        self.Ok = 1. - Om - Ol         # curvature density
        self.w = w                     # pressure / density for dark energy
        self.H0 = H0                   # present day Hubble constant, km/s/Mpc
        self.h = H0 / 100.             # 100 km/s/Mpc * h = H0
        self.h70 = H0 / 70.            # 70 km/s/Mpc * h70 = H0
        self.Th = 1. / H0 * PC * 1.e3  # Hubble time in s
        self.Dh = C / H0 * PC * 1.e6   # Hubble distance in m

    def __repr__(self):
        return "Cosmology(Om=%(Om)s, Ol=%(Ol)s, Ok=%(Ok)s, H0=%(H0)s, w=%(w)s)" % self.__dict__

    def efunc(self,z):
        """ Function for integration. Eqn 14 from Hogg, adding in 'w'
        to Omega lambda."""
        zp1 = 1. + z
        temp = 3. * (1. + self.w)
        return sqrt(self.Om*zp1**3 + self.Ok*zp1**2 + self.Ol*zp1**temp)

    def inv_efunc(self,z):
        """ Function for integration. Eqn 14 from Hogg, adding in 'w'
        to Omega lambda."""
        zp1 = 1. + z
        temp = 3. * (1. + self.w)
        return 1. / sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                         + self.Ol*zp1**temp)

    def tfunc(self,z):
        """ Function for integration to find lookback time. Eqn 30
        from Hogg, adding in 'w' to Omega lambda."""
        zp1 = 1. + z
        temp = 3.*(1. + self.w)
        return 1. / (zp1*sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                              + self.Ol*zp1**temp))

    def xfunc(self,z):
        """ Function for integration to find the absorption distance.
        """
        zp1 = 1. + z
        temp = 3.*(1. + self.w)
        return zp1**2 / sqrt(self.Om*zp1**3 + self.Ok*zp1**2
                             + self.Ol*zp1**temp)

    def Hz(self,z):
        """ Returns the Hubble constant at redshift z. """
        return self.H0 * self.efunc(z)

    def a(self,z):
        """ Returns the scale factor at z,   1 / (1+z). """
        return 1. / (1.+z)

    def Tl(self,z):
        """ Returns the lookback time in seconds; the difference
        between the age of the Universe now and the age at z."""
        return self.Th * integrate.quad(self.tfunc, 0, z)[0]

    def Dc(self,z):
        """ Returns the comoving distance to an object at redshift z
        in metres. Remains constant with time if the objects are in
        the Hubble flow."""
        return self.Dh * integrate.quad(self.inv_efunc, 0, z)[0]

    def Dm(self,z):
        """ Returns the transverse comoving distance in metres.
        Dm*dtheta gives the transverse comoving distance between two
        objects separated by an angle dtheta radians at redshift
        z. Note Dm = Dc if Ok is zero (as in current lambda CDM
        model)."""
        Ok = self.Ok
        Dc  = self.Dc(z)
        if Ok == 0:
            return Dc
        sqrtOk = sqrt(abs(Ok))
        Dh  = self.Dh
        if Ok > 0:
            return Dh / sqrtOk * sinh(sqrtOk * Dc / Dh)
        else:
            return Dh / sqrtOk * sin(sqrtOk * Dc / Dh)

    def Da(self,z):
        """ Angular diameter distance in metres.  Ratio of an object's
        physical transverse size to its angular size in radians."""
        return self.Dm(z) / (1. + z)

    def Da2(self, z1, z2):
        """ Angular diameter distance in metres between objects at 2
        redshifts.  Useful for gravitational lensing."""
        # does not work for negative curvature
        assert(self.Ok) >= 0

        # z1 < z2
        if (z2 < z1): z1,z2 = z2,z1

        Dm1 = self.Dm(z1)
        Dm2 = self.Dm(z2)
        Ok  = self.Ok
        Dh_2  = self.Dh * self.Dh

        return 1. / (1.+z2) * (Dm2*sqrt(1. + Ok*Dm1**2 / Dh_2) -
                               Dm1*sqrt(1. + Ok*Dm2**2 / Dh_2))

    def Dl(self,z):
        """ Returns the luminosity distance in metres.  (Relationship
        between bolometric flux and bolometric luminosity.)"""
        return (1.+z) * self.Dm(z)

    def distmod(self,z):
        """ Returns the distance modulus (apparent magnitude -
        absolute magnitude for an object at redshift z)."""
        # Remember that Dl is in m
        return 5. * log10(self.Dl(z) / PC / 10.)

    def X(self, z2, z1=0):
        """ Return the absorption distance from redshift
        z2 to z1. Dimensionless (despite being called a distance) """
        return integrate.quad(self.xfunc, z1, z2)[0]

# convenience functions
def kpc_per_arcmin(z, cosmo=None):
    """ Find the number of transverse kpc corresponding to an
    arcminute at the given redshift.
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    arcmin_in_radians = 1 / 60. * pi / 180
    return cosmo.Dm(z) / 1.e3 / PC * arcmin_in_radians

# convenience functions
def arcsec_per_kpc(z, cosmo=None):
    """ Find the number of transverse kpc corresponding to an
    arcminute at the given redshift.
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    arcsec_in_radians = 1 / 3600. * pi / 180
    return 1 / (cosmo.Dm(z) / 1.e3 / PC * arcsec_in_radians)


def distmod(z, cosmo=None):
    """Find the distance modulus at the given redshift.
    """
    if cosmo is None:
        cosmo = Cosmology(Om=0.27, Ol=0.73, w=-1., H0=70.)
    return cosmo.distmod(z)

def to_xyz(ra, dec, r, deg2rad=np.pi/180.):
    """ Convert a comoving distance, ra and dec into 3d comoving
    coordinates.

    Convert a vector pointing from the origin with comoving length
    r in the direction ra,dec (in degrees) to comoving 3-d
    coordinates x,y,z.

    Assumes universe is flat.

    ra = 0 corresponds to the positive x-z half plane. dec = 0
    corresponds to the whole x-y plane. If this is confusing, draw a
    picture :)

    To align the axes such that a vector r with direction (ra0, dec0)
    points along the positive x-axis, use ra-ra0 and dec-dec0 instead
    of ra and dec. In this case increasing dec corresponds to
    increasing z, and increasing ra (from ra=0) corresponds to
    increasing y.
    """
    # radians
    ra1 =  deg2rad * ra
    dec1 =  deg2rad * dec

    cos_dec1 = np.cos(dec1)
    X = r * (cos_dec1 * np.cos(ra1))
    Y = r * (cos_dec1 * np.sin(ra1))
    Z = r * np.sin(dec1)

    return X,Y,Z


if __name__ == '__main__':

    c = Cosmology()

    if len(sys.argv) < 2:
        raise Exception('Usage : cosmology.py z1 [z2]')
    z1 = float(sys.argv[1])

    print 'Cosmology : H0           =', c.H0
    print 'Cosmology : Omega Matter =', c.Om
    print 'Cosmology : Omega Lambda =', c.Ol
    print ''

    print 'Hubble distance                %.2f Mpc' % (c.Dh / PC / 1.e6)
    print 'Hubble time                    %.2f Gyr' % (c.Th/3600./24./
                                                       365.25/1.e9)
    print ''
    print 'For z = %.2f:' % (z1)
    print 'Lookback time                  %.2f Gyr' % (c.Tl(z1)/3600./24./
                                                       365.25/1.e9)
    print 'Scale Factor a                 %.2f' % (c.a(z1))
    print 'Comoving L.O.S. Distance (w)   %.2f Mpc' % (c.Dc(z1)/PC/1.e6)
    print 'Angular diameter distance      %.2f Mpc' % (c.Da(z1)/PC/1.e6)
    print 'Luminosity distance            %.2f Mpc' % (c.Dl(z1)/PC/1.e6)
    print 'Distance modulus               %.2f mag' % (c.distmod(z1))
