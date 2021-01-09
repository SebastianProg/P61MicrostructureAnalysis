import numpy as np

import BasicFunctions.generalFunctions as gf
import BasicFunctions.generalCalculations as gc


# Determine linear regression of given values
# parameters:
# values of both signals
# returns:
# increasing, constant offset, coefficient of determination, determined values
def linReg(x, y, fixParName='', fixPar=0):  # scipy.stats.linregress(x, y=None)
	n = len(x)
	x = np.array(x)
	y = np.array(y)
	# check dimensions of values
	if gf.size(x, 0) != gf.size(y, 0) or gf.size(x, 1) != gf.size(y, 1):
		# try to transpose data
		if gf.size(x,0) == gf.size(np.transpose(y), 0) and gf.size(x, 1) == gf.size(np.transpose(y), 1):
			y = np.transpose(y)
	sumX = sum(x)
	sumY = sum(y)
	meanX = sumX/n
	meanY = sumY/n
	sumXY = sum(x * y)
	sumXX = sum(x ** 2)
	if fixParName == "b":
		b = fixPar
		m = (sumXY - b * sumX) / sumXX
	elif fixParName == "m":
		m = fixPar
		b = meanY - m * meanX
	else:
		m = (n * sumXY - sumX * sumY) / (n * sumXX - sumX**2)
		b = meanY - m * meanX
		#b = (sumXX * sumY - sumX * sumXY) / (n * sumXX - sumX**2)
		#m = (sumXY - b * sumX) / sumXX
	yD = m * x + b
	CoD = 1 - sum((y - yD) ** 2) / sum((y - meanY) ** 2)
	err = sum((y - yD) ** 2) / (n * (n - 1)) # variance
	mErr = err / sum((x - meanX)**2)
	bErr = err / n * sum(x**2) / sum((x - meanX)**2)
	return m, b, CoD, yD, err, mErr, bErr


# Determine linear regression of given values with accuracy
# parameters:
# values of both signals and accuracy of values
# returns:
# increasing, constant offset, coefficient of determination, determined values
def linRegWeighted(x, y, a, fixParName='', fixPar=0):
	xNew, yNew = gc.createWeightedData(x, y, a)
	# linear regression with new data set
	m, b, CoD, yD, err, mErr, bErr = linReg(xNew, yNew, fixParName, fixPar)
	yD = m * x + b
	return m, b, CoD, yD, err, mErr, bErr


def distributionAnalysis(x, y):
	expVal = sum(y * x) / sum(y)
	stdVal = (sum(y * (x - expVal) ** 2) / sum(y)) ** 0.5
	return expVal, stdVal
