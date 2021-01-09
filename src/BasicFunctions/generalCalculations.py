import numpy as np
import decimal as dec

import BasicFunctions.generalFunctions as gf


def sind(x):
	return np.sin(np.radians(x))


def cosd(x):
	return np.cos(np.radians(x))


def tand(x):
	return np.tan(np.radians(x))


def cotd(x):
	return 1/tand(x)


def asind(x):
	return np.degrees(np.arcsin(x))


def acosd(x):
	return np.degrees(np.arccos(x))


def atand(x):
	return np.degrees(np.arctan(x))


def atan2d(y, x):
	return np.degrees(np.arctan2(y, x))


def rad2deg(x):
	return np.radians(x)


def deg2rad(x):
	return np.degrees(x)


# Creates equidistant values in interval [a,b] with number of n
# a: lowest value, b: highest value, n: number of values
# x: array of equidistant values of number n in interval [a,b]
def createRange(a, b, n):
	return list(range(a, b, (b - a) // (n - 1)))


def distance(a, b):
	return sum(abs(np.array(a) - np.array(b)))


def roundDigits(values, numDecimal=0, exact=False, classicRounding=False):
	if exact and numDecimal < 0:
		if classicRounding:
			return np.array([float(dec.Decimal(i).quantize(dec.Decimal('1'),
				rounding=dec.ROUND_HALF_UP) * dec.Decimal(-numDecimal)) for i in (np.array(values) / -numDecimal)])
		else:
			return np.round(np.array(values) / -numDecimal) * -numDecimal
	else:
		if classicRounding:
			return np.array([float(dec.Decimal(i).quantize(dec.Decimal('1'),
				rounding=dec.ROUND_HALF_UP) / dec.Decimal(10 ** numDecimal)) for i in (10 ** numDecimal * np.array(values))])
		else:
			return np.round(10 ** numDecimal * np.array(values)) / 10 ** numDecimal


def transfInterval(x, a=None, b=None, c=0, d=1):
	# x (scalar or vector) of [a b] will be transformed to [c d]
	if a is None:
		a = np.min(x)
	if b is None:
		b = np.max(x)
	return (x - a) * (d - c) / (b - a) + c


def createWeightedData(x, y, a):
	# check dimensions of values
	if gf.size(x, 0) != gf.size(y, 0) or gf.size(x, 1) != gf.size(y, 1):
		if gf.size(x, 0) == gf.size(a, 0) and gf.size(x, 1) == gf.size(a, 1):
			# try to transpose data
			if gf.size(x, 0) == gf.size(np.transpose(y), 0) and gf.size(x, 1) == gf.size(np.transpose(y), 1):
				y = np.transpose(y)
		elif gf.size(y, 0) == gf.size(a, 0) and gf.size(y, 1) == gf.size(a, 1):
			# try to transpose data
			if gf.size(x, 0) == gf.size(np.transpose(y), 0) and gf.size(x, 1) == gf.size(np.transpose(y), 1):
				y = np.transpose(y)
				a = np.transpose(a)
	# check validity of weigths
	a[a < 0] = 1
	# create new data set according to weights
	minA = np.min(a)
	weights = np.ceil(a / minA)
	n = int(np.sum(weights))
	xNew = gf.ones(n)
	yNew = gf.ones(n)
	pos = 0
	for i in range(len(weights)):
		xNew[pos:pos + int(weights[i])] = x[i]
		yNew[pos:pos + int(weights[i])] = y[i]
		pos = pos + int(weights[i])
	return xNew, yNew


def strrepAll(origStr, oldSubstr, newSubstr):
	strVal = origStr
	for i in range(len(oldSubstr)):
		return strVal.replace(oldSubstr[i], newSubstr)
