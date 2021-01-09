import numpy as np

import BasicFunctions.generalCalculations as gc


def hklroot(h, k, l):
	return (h ** 2 + k ** 2 + l ** 2) ** 0.5


def hklroot2(hkl):
	vals = splitHkl(hkl)
	return hklroot(vals[0], vals[1], vals[2])


def angles2latticeDists(angles, wavelength):
	return wavelength / (2 * gc.sind(angles / 2))


def wavelengths2energies(wavelengths):
	h = 4.135667662e-15  # eVs 6.626070040e-34 Js
	c = 2.99792458e8  # m/s
	return h * c * 10 ** 6 / wavelengths  # keV


def energies2wavelengths(energies):
	h = 4.135667662e-15  # eVs 6.626070040e-34 Js
	c = 2.99792458e8  # m/s
	return h * c * 10 ** 6 / energies  # nm


def energies2latticeDists(energies, diffAngle):
	return energies2wavelengths(energies) / (2 * gc.sind(diffAngle / 2))  # nm


def latticeDists2energies(latticeDists, diffAngle):
	# wavelengths2energies can be also used for lattice distances
	return wavelengths2energies(latticeDists) / (2 * gc.sind(diffAngle / 2))  # keV


def latticeDists2angles(latticeDists, wavelength):
	angles = 2 * gc.asind(wavelength / (2 * latticeDists))
	useful = np.nonzero(np.imag(angles) == 0)
	return angles, useful


def latticeDists2aVals(d, h, k, l):
	return d * hklroot(h, k, l)


def latticeDists2aVals2(d, hkl):
	return (d * hklroot2(hkl))


def aVals2latticeDists2(aVal, hkl):
	return aVal / hklroot2(hkl)


def aVals2latticeDists(aVal, h, k, l):
	return aVal / hklroot(h, k, l)


def strains2dVals(strains, d0):
	return d0 * (strains + 1)


def dVals2strains(dVals, d0):
	return dVals / d0 - 1


def channels2energies(channels, par, par1=None):
	if par1 is None or len(par) == 3:
		energies = par[0] * channels**2 + par[1] * channels + par[2]
	else:
		energies = par * channels + par1
	return energies


def calcModulus(s1, hs2):
	e = 1 / (s1 + hs2)
	nu = hs2 / (s1 + hs2) - 1
	return e, nu


def calcDec(e, nu):
	s1 = -nu / e
	hs2 = (1 + nu) / e
	return s1, hs2


def splitHkl(hkl):  # only useful for positive h, k, l and k and l smaller than 10
	l = hkl % 10
	k = (hkl // 10) % 10
	h = hkl // 100
	return h, k, l


def splitHkl2(hkl):  # only useful for positive h, k, l and k and l smaller than 10
	l = hkl % 10
	k = (hkl // 10) % 10
	h = hkl // 100
	return np.transpose(np.array([h, k, l]))


def mergeHkl(h, k, l):  # works only for positive h, k, l and k and l smaller than 10
	return h * 100 + k * 10 + l


def mergeHkl2(hklVector):  # works only for positive h, k, l and k and l smaller than 10
	return mergeHkl(hklVector[:,0], hklVector[:,1], hklVector[:,2])
