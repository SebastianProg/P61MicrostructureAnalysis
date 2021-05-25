import math
import numpy as np
import scipy.optimize as op
import scipy.special as spe
import os
import shutil as util

import basics.functions as bf
import diffraction.calculations as dc
import diffraction.conversions as conv
import plotting.specific as ps
import basics.datamodels as dm
import filehandling.general as fg
import filehandling.specific as fs


def initPeakfitSettings(*args):
	# default values
	settings = {"peakFunction": "Gauss", "withGraphic": 1, "peakZones": 1, "withCorrection": 0, "withFourier": 0,
		"bestWindow": 0, "onlyValues": 0, "toEnergies": 0, "par": 0, "parDeadTime": 0, "saveExcel": 0,
		"saveText": 1, "confInterval": 0}  # 'Gauss', 'Pearson7', 'Lorentz', 'PseudoVoigt'
	if len(args) % 2 == 0:
		# set specified values
		for i in range(0, len(args), 2):
			if bf.hasKey(settings, args[i]):  # adding new entries is currently not allowed
				settings[args[i]] = args[i + 1]
	elif len(args) == 1:
		# file with settings is specified
		lines = fg.readLines(args[0])
		for line in lines:
			lineData = line.split('\t')
			if lineData[1].startswith('['):
				settings[lineData[0]] = bf.strArray2npArray(lineData[1])
			else:
				settings[lineData[0]] = eval(lineData[1])
	return settings


def estimate_peak_params(xx, yy, test_ind=5, consider_background=True):
	# test_ind defines distance to maximum value to estimate sigma parameter
	max_ind = np.argmax(yy)
	max_pos = xx[max_ind]
	max_val = yy[max_ind]
	l_sigma = 0
	r_sigma = 0
	if max_ind > test_ind:
		l_sigma = np.abs(xx[max_ind - test_ind] - xx[max_ind]) / (
				2 * (np.log(max_val) - np.log(yy[max_ind - test_ind]))) ** 0.5
		est_sigma = l_sigma
		l_used = (xx > max_pos - 3.1 * l_sigma) & (xx < max_pos - 2.9 * l_sigma)
		if np.any(l_used):
			x_used = xx[l_used]
			lx = x_used[int(len(x_used) / 2)]
			l_val = np.mean(yy[l_used])
		else:
			lx = xx[int(test_ind / 2)]
			l_val = np.mean(yy[0:test_ind])
	if max_ind < len(xx):
		r_sigma = np.abs(xx[max_ind + test_ind] - xx[max_ind]) / (
					2 * (np.log(max_val) - np.log(yy[max_ind + test_ind]))) ** 0.5
		est_sigma = r_sigma
		r_used = (xx > max_pos + 2.9 * r_sigma) & (xx < max_pos + 3.1 * r_sigma)
		if np.any(r_used):
			x_used = xx[r_used]
			rx = x_used[int(len(x_used) / 2)]
			r_val = np.mean(yy[r_used])
		else:
			rx = xx[len(xx) - 1 - int(test_ind / 2)]
			r_val = np.mean(yy[-test_ind - 1:-1])
	if l_sigma > 0 and r_sigma > 0:
		est_sigma = np.mean([l_sigma, r_sigma])
	if consider_background:
		m_val = (r_val - l_val) / (rx - lx)
		b_val = np.mean([l_val - m_val * lx, r_val - m_val * rx])
		max_val, max_pos, est_sigma = estimate_peak_params(xx, yy - (m_val * xx + b_val), test_ind, False)
		return max_val, max_pos, est_sigma, m_val, b_val
	else:
		return max_val, max_pos, est_sigma


def evalPeakData(settings, peakData, n=1):
	# check values of settings -> error handling, default values
	handle = -1  # handle of figure
	# IO_Korrektur Z. 25
	# Fourierfilter Z. 29
	# evaluate each peak
	analysisData = np.zeros((bf.size(settings["peakZones"], 0), 11))
	for i in range(bf.size(settings["peakZones"], 0)):
		x = np.array(peakData[:, 0])  # channels, energies or angles
		y = np.array(peakData[:, 1])  # intensities
		#determine start and end value in actual data
		if bf.length(settings["peakZones"]) == 1 and bf.size(settings["peakZones"], 0) == bf.size(settings["peakZones"], 1):
			startI = 0
			stopI = len(peakData[:, 1]) - 1  # all data
		else:
			startI = np.argwhere(x >= settings["peakZones"][i, 0])[0]  # start value
			startI = startI[0]
			stopI = np.argwhere(x <= settings["peakZones"][i, 1])[-1]  # end value
			stopI = stopI[0]
		xpart = x[startI:stopI+1]
		ypart = y[startI:stopI+1]
		#fit current peak
		filename = 'Daten' + str(n) + '.xls'
		if settings["peakFunction"] == 'Gauss':
			func = gaussfit
		elif settings["peakFunction"] == 'Pearson7':
			func = pearson7fit
		elif settings["peakFunction"] == 'Lorentz':
			func = lorentzfit
		elif settings["peakFunction"] == 'PseudoVoigt':
			func = pseudovoigtfit
		if settings["bestWindow"] == 1:
			[par, startI, stopI] = dc.bestwindowfit(func, x, y, startI, stopI)
			xpart = x[startI:stopI+1]  # new relevant x values
			ypart = y[startI:stopI+1]  # new relevant y values
			[par, par0, yMax, yFit, peakQuality, err, meanNoise] = func(xpart, ypart, par)
		else:
			[par, par0, yMax, yFit, peakQuality, err, meanNoise] = func(xpart, ypart)
		analysisData[i, :] = np.concatenate(([xpart[0], xpart[-1]], par[0:3], [yMax, peakQuality, err], par[3:6]))
		if settings["saveText"] == 1 or settings["saveExcel"] == 1:
			#save results in Excel
			headerStr = "Peak-Start\tPeak-End\tPeakarea-A\tPeakwidth-s\tPeakpos-xc\tPeakintMax\tPeakquality\tDeltaVal\tBackSlope-m\tBackVal-b\tpar6"
			fg.excelsave(filename, headerStr, [analysisData[i, :]])
		#show results
		if settings["withGraphic"] == 1:
			handle = ps.picture(bf.size(settings["peakZones"], 0), yFit, xpart, ypart, n, i)
	return handle, analysisData


def evalPeakFile(settings, fileName, n=1):
	#import data
	[curData, curHeader] = fs.getSpectrumFromFile(fileName)
	if settings["toEnergies"] != 0:
		if settings["parDeadTime"] != 0 and settings["par"] != 0:
			#determine one parameter depending on dead time value
			livetime = (bf.replace(curHeader[2], 'TLIVE=', '')).float
			realtime = (bf.replace(curHeader[3], 'TREAL=', '')).float
			deadTime = (1 - livetime / realtime) * 100 #in %
			settings["par"][2] = dc.calcEnergyCalibDeadTime(settings["parDeadTime"], deadTime)
		#convert channels to energies
		curData[:, 0] = conv.channels2energies(curData[:, 0], settings["par"])
	elif settings["onlyValues"] != 0:
		curData[:, 0] = np.array(range(1, len(curData[:, 0]) + 1))
	measureData = curData
	#check values of settings -> error handling, default values
	[h, analysisData] = evalPeakData(settings, curData, n)
	return h, analysisData, measureData, curHeader


def evalPeakFiles(settings, fileNames, saveFile=""):
	# check values of settings -> error handling
	#measureData = zeros((length(fileNames),8192,2))
	measureDataCell = []	#cell(len(fileNames),1)
	analysisData = np.zeros((len(fileNames), bf.size(settings["peakZones"], 0), 11))
	headers = []	#cell(len(fileNames),1)
	for i in range(len(fileNames)):
		# get current file
		curFileName = fileNames[i]
		[h, curAnalysisData, curMeasureData, curHead] = evalPeakFile(settings, curFileName, i)
		measureDataCell.append(curMeasureData)
		analysisData[i, :, :] = curAnalysisData
		headers.append(curHead)
		'''# wait one second and then close current figure with peak fits
		if settings["withGraphic"] == "1" and len(fileNames) > 1:
			plt.pause(1)
			plt.close(h)
		'''
	#transform measure data
	measureData = []
	#measureData = np.zeros((len(measureDataCell),bf.size(measureDataCell[0])))
#	for i = 1:length(measureDataCell)
#		measureData(i,:,:) = measureDataCell{i}
	#save results to files
	if bf.isempty(saveFile):
		saveFile = fg.requestSaveFile((("All files", "*.*"),), 'Save as...')
		if bf.isempty(saveFile):
			return
	resFiles = []
	for i in range(len(fileNames)):
		sourceFile = 'Daten' + str(i) + '.xls'
		if fg.existsFile(sourceFile):
			pathSrc, nameSrc = fg.fileparts(fileNames[i])
			pathDest, nameDest = fg.fileparts(saveFile)
			if settings["saveExcel"] == 1:
				resFiles.append(nameDest + nameSrc)
			if settings["saveText"] == 1:
				resFiles.append(nameDest + nameSrc.replace(".xls", ".txt"))
			if len(pathDest) > 0 and pathDest != "/":
				resFiles[i] = pathDest + '/' + resFiles[i]
			util.copyfile(sourceFile, resFiles[i])
			os.remove(sourceFile)
	return analysisData, measureData, headers, resFiles


def gaussfit(xdata, ydata, x0=None):
	[m, b, CoD, yD, err, mErr, bErr] = dm.linReg(np.concatenate((xdata[0:5], xdata[-5:])), np.concatenate((ydata[0:5],
		ydata[-5:])))
#     [yMax i] = max(ydata); % search maximum value
	[yMax, i] = bf.max(ydata - (m * xdata + b))
	if x0 is None:
		yLeft = m * xdata[0] + b
		yRight = m * xdata[-1] + b
		ySum = np.sum(ydata * np.abs(np.gradient(xdata))) - (0.5 * np.abs((xdata[-1] - xdata[0]) * (yRight - yLeft)) +
			np.abs(bf.min(yLeft, yRight)[0] * (xdata[-1] - xdata[0])))  # estimated area of peak
		xMax = xdata[i]  # find belonging x value
		sigma = ySum / (yMax * (2 * np.pi)**0.5)  # estimate width of peak
		x0 = [ySum, sigma, xMax, m, b, 1.0]  # start values Gauss
	x_par = op.leastsq(gaussopt, x0, (xdata, ydata))
	x = x_par[0]
	# x = op.curve_fit(gaussfunction2, xdata, ydata, x0)
	yFit = gaussfunction(x, xdata)
	# check quality of peak
	meanNoise = (yLeft + yRight) / 2
	peakQuality = (yMax + meanNoise)/(meanNoise + 3 * meanNoise**0.5)
	err = 0
	# adapt dimensions of x if x0 has different dimensions
	if bf.size(x0) != bf.size(x):
		x = np.transpose(x)
	return x, x0, yMax, yFit, peakQuality, err, meanNoise


def pearson7fit(xdata, ydata, x0=None):
	[m, b, CoD, yD, err, mErr, bErr] = dm.linReg(np.concatenate((xdata[0:5], xdata[-5:])),
		np.concatenate((ydata[0:5], ydata[-5:])))
	[yMax, i] = bf.max(ydata - (m * xdata + b))
	if x0 is None:
		yLeft = m * xdata[1] + b
		yRight = m * xdata[-1] + b
		ySum = np.sum(ydata * np.abs(np.gradient(xdata))) - (0.5 * np.abs((xdata[-1] - xdata[0]) * (yRight - yLeft)) +
			np.abs(bf.min(yLeft, yRight)[0] * (xdata[-1] - xdata[0])))  # estimated area of peak
		xMax = xdata[i]  # find belonging x value
		sigma = ySum / (yMax * (2*np.pi)**0.5) # estimate width of peak
		x0 = [ySum, sigma, xMax, m, b, 2.0]  # start values Pearson7
	else:
		x0 = [x0[0], x0[1], x0[2], x0[3], x0[4], 1]  # start values Pearson7
	x_par = op.leastsq(pearson7opt, x0, (xdata, ydata))
	x = x_par[0]
	yFit = pearson7function(x, xdata)
	# adapt x to fit peak description
	x = [x[0], x[1]/(2 * x[5] - 3)**0.5, x[2], x[3], x[4], x[5]]
	# check quality of peak
	meanNoise = (yLeft + yRight) / 2
	peakQuality = (yMax + meanNoise)/(meanNoise + 3 * meanNoise**0.5)
	err = 0
	# adapt dimensions of x if x0 has different dimensions
	if bf.size(x0) != bf.size(x):
		x = np.transpose(x)
	return x, x0, yMax, yFit, peakQuality, err, meanNoise


def lorentzfit(xdata, ydata, x0=None):
	[m, b, CoD, yD, err, mErr, bErr] = dm.linReg(np.concatenate((xdata[0:5], xdata[-5:])),
		np.concatenate((ydata[0:5], ydata[-5:])))
	[yMax, i] = bf.max(ydata - (m * xdata + b))
	if x0 is None:
		yLeft = m * xdata[0] + b
		yRight = m * xdata[-1] + b
		ySum = sum(ydata * np.abs(np.gradient(xdata))) - (0.5 * np.abs((xdata[-1] - xdata[0]))) * (yRight - yLeft) + \
			np.abs(bf.min(yLeft, yRight)[0] * (xdata[-1] - xdata[0]))  # estimated area of peak
		xMax = xdata[i]  # find belonging x value
		sigma = ySum / (yMax * np.pi)  # estimate width of peak
		x0 = [ySum, sigma, xMax, m, b, 1.0]  # start values Gauss
	x_par = op.leastsq(lorentzopt, x0, (xdata, ydata))
	x = x_par[0]
	yFit = lorentzfunction(x, xdata)
	# check quality of peak
	meanNoise = (yLeft + yRight) / 2
	peakQuality = (yMax + meanNoise)/(meanNoise + 3 * meanNoise**0.5)
	err = 0
	# adapt dimensions of x if x0 has different dimensions
	if bf.size(x0) != bf.size(x):
		x = np.transpose(x)
	return x, x0, yMax, yFit, peakQuality, err, meanNoise


def pseudovoigtfit(xdata, ydata, x0=None):
	[m, b, CoD, yD, err, mErr, bErr] = dm.linReg(np.concatenate((xdata[0:5], xdata[-5:])),
		np.concatenate((ydata[0:5], ydata[-5:])))
	[yMax, i] = bf.max(ydata - (m * xdata + b))
	if x0 is None:
		yLeft = m * xdata[0] + b
		yRight = m * xdata[-1] + b
		ySum = np.sum(ydata * np.abs(np.gradient(xdata))) - (0.5 * np.abs((xdata[-1] - xdata[0]) * (yRight - yLeft)) +
			np.abs(bf.min(yLeft, yRight)[0] * (xdata[-1] - xdata[0])))  # estimated area of peak
		xMax = xdata[i]  # find belonging x value
		sigma = ySum / (yMax * np.pi)  # estimate width of peak
		x0 = [ySum, sigma, xMax, m, b, 0.5]  # start values Pseudo-Voigt
	x_par = op.leastsq(pseudovoigtopt, x0, (xdata, ydata))
	x = x_par[0]
	yFit = pseudovoigtfunction(x, xdata)
	# check quality of peak
	meanNoise = (yLeft + yRight) / 2
	peakQuality = (yMax + meanNoise)/(meanNoise + 3 * meanNoise**0.5)
	err = 0
	# adapt dimensions of x if x0 has different dimensions
	if bf.size(x0) != bf.size(x):
		x = np.transpose(x)
	return x, x0, yMax, yFit, peakQuality, err, meanNoise


def gaussfunction2(x, par1, par2, par3, par4, par5, par6):
	return gaussfunction([par1, par2, par3, par4, par5, par6], x) 


def gaussfunction(par, x):
	if bf.size(par, 1) < 6:
		par = np.transpose(par)
	f = lambda para, x: (para[0] / (para[1] * ((2 * math.pi)**0.5))) * np.exp(-0.5 * (((x - para[2])**2) /
		(para[1]**2))) + (para[3] * x + para[4])  # gauss function
	return f(par, x)


def lorentzfunction(par, x):
	if bf.size(par, 1) < 6:
		par = np.transpose(par)
	# lorentz function
	f = lambda para, x: (para[0] * para[1] / (np.pi * (para[1]**2 + (x - para[2])**2)) + para[3] * x + para[4])
	return f(par, x)


def pseudovoigtfunction(par, x):
	if bf.size(par, 1) < 6:
		par = np.transpose(par)
	# lorentz function
	fL = lambda para, x: (para[0] * para[1] / (math.pi * (para[1]**2 + (x - para[2])**2)))
	# gauss function
	fG = lambda para, x: (para[0] / (para[1] * ((2 * math.pi)**0.5))) * np.exp(-0.5 * (((x - para[2])**2) / (para[1]**2)))
	f = lambda para, x: para[5] * fL(para,x) + (1 - para[5]) * fG(para,x) + (para[3] * x + para[4])
	return f(par, x)


def pearson7function(par, x):
	if bf.size(par, 1) < 6:
		par = np.transpose(par)
	f = lambda para, x: para[0] * (1 + ((x - para[2]) / para[1])**2)**(-para[5]) / (para[1] *
		spe.beta(para[5] - 1/2, 1/2)) + (para[3] * x + para[4])
	return f(par, x)


def gaussopt(par, x, y):
	return gaussfunction(par, x) - y


def lorentzopt(par, x, y):
	return lorentzfunction(par, x) - y


def pseudovoigtopt(par, x, y):
	return pseudovoigtfunction(par, x) - y


def pearson7opt(par, x, y):
	return pearson7function(par, x) - y
