# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:53:13 2019

@author: Sebastian
"""

from io import StringIO  # for handling unicode strings
import numpy as np
import os

import basics.functions as bf
import filehandling.general as fg
import diffraction.conversions as conv
import diffraction.calculations as dc

FILE_SAMPLE_MATERIALS = 'system_files/Probenmaterial.txt'
FILE_XRAY_SOURCES_ENERGIES = 'system_files/StrahlungsquellenEnergien_keV.txt'
FILE_XRAY_SOURCES_WAVELENGTHS = 'system_files/StrahlungsquellenWellenlaengen_nm.txt'


def writeSpectrumFile(file, data, header=None, preheader='', newlineStr=None):
	if newlineStr is None:
		newlineStr = os.linesep
	fid = fg.openFile(file, 'w')
	if header is not None:
		for i in range(len(header)):
			fg.writeLine(fid, ('%s' % (preheader + header[i])))
	fid.close()
	fid = open(file, 'ab')
	np.savetxt(fid, data, delimiter='\t', newline=newlineStr, fmt='%g')
	fid.close()


def getSpectrumFromFile(specFile="", dataDelim='\t', headerLines=5):
	data = 0
	heading = 0
	if len(specFile) == 0:
		specFile = fg.requestFiles((('Text file', '*.txt'),), 'Spektrumdatei auswaehlen', 'off')
	if len(specFile) > 0:
		#Read data from file
		heading = fg.readLines(specFile, headerLines)
		data = fg.dlmread(specFile, dataDelim, headerLines)
	return data, heading


def accumSpectra(resFile='', files=None, method='sum', dataDelim='\t', headerLines=5):
	# method: 'sum', 'mean', 'max', 'min'
	results = None
	header = []
	if files is None:
		files = fg.requestFiles((('Text file', '*.txt'),), 'Spektrumdateien auswaehlen', 'on')
	if len(files) > 1:
		for i in range(len(files)):
			[data, heading] = getSpectrumFromFile(files[i], dataDelim, headerLines)
			# initialize values if not existing
			if results is None:
				if method == 'min':
					results = np.ones(data.shape) * 1e8
					tlive = 1e6
					treal = 1e6
				else:
					results = np.zeros(data.shape)
					tlive = 0
					treal = 0
				results[:, 0] = data[:, 0]
			# get current live and real time
			if len(heading) > 2 and 'TLIVE' in heading[2]:
				tliveCur = float(bf.replace(heading[2], 'TLIVE=', ''))
			else:
				tliveCur = tlive
			if len(heading) > 3 and 'TREAL' in heading[3]:
				trealCur = float(bf.replace(heading[3], 'TREAL=', ''))
			else:
				trealCur = treal
			# update values
			if method == 'sum' or method == 'mean':
				results[:, 1] = results[:, 1] + data[:, 1]
				tlive = tlive + tliveCur
				treal = treal + trealCur
			elif method == 'max':
				results[:, 1] = np.max((results[:, 1], data[:, 1]), axis=0)
				tlive = np.max((tlive, tliveCur))
				treal = np.max((treal, trealCur))
			elif method == 'min':
				results[:, 1] = np.min((results[:, 1], data[:, 1]), axis=0)
				tlive = np.min((tlive, tliveCur))
				treal = np.min((treal, trealCur))
		if method == 'mean':
			results[:, 1] = results[:, 1] / len(files)
			tlive = tlive / len(files)
			treal = treal / len(files)
		if headerLines > 0:
			header = ['DATE=' + bf.datetimestr('%d-%m-%Y'), 'TIME=' + bf.datetimestr('%H:%M:%S'), 'TLIVE=' + str(tlive),
				'TREAL=' + str(treal), heading[4]]
		if len(resFile) > 0:
			if headerLines > 0:
				writeSpectrumFile(resFile, results, header[:(min(len(header), headerLines))])
			else:
				writeSpectrumFile(resFile, results)
	return results, header, files




def getDecVals(decFile, hklVals=None):
	decVals = 0
	if decFile is None or len(decFile) == 0:
		fileNames = fg.requestFiles((('Text file', '*.txt'),), 'Select DEC file', "off")
		decFile = fileNames[0]
	if decFile is not None:
		# Convert data to matrix
		data = fg.dlmread(decFile, '\t', 1)
		if hklVals is not None:
			#select relevant values
			return dc.selectDecVals(data, hklVals)
		else:
			return data


def loadXrayData():
	#Kalpha1, Kalpha2, Kbeta1, KalphaMean, AbsorptionK, Voltage
	elementWavelength = fg.readLines(FILE_XRAY_SOURCES_WAVELENGTHS)
	header = elementWavelength[0].strip().split("\t")
	xRayData = {}
	for line in elementWavelength[1:]:
		newElementWavelength = line.split("\t")
		xRayData[newElementWavelength[0]] = {header[1]: int(newElementWavelength[1]),
			header[6]: float(newElementWavelength[6]), header[7]: float(newElementWavelength[7])}
		xRayData[newElementWavelength[0]][header[2]] = {"Wavelength": float(newElementWavelength[2])}
		xRayData[newElementWavelength[0]][header[3]] = {"Wavelength": float(newElementWavelength[3])}
		xRayData[newElementWavelength[0]][header[4]] = {"Wavelength": float(newElementWavelength[4])}
		xRayData[newElementWavelength[0]][header[5]] = {"Wavelength": float(newElementWavelength[5])}
	elementEnergies = fg.readLines(FILE_XRAY_SOURCES_ENERGIES)
	headerEg = elementEnergies[0].strip().split("\t")
	for lineEg in elementEnergies[1:]:
		newElementEnergies = lineEg.split("\t")
		#xRayData[newElementEnergies[0]][headerEg[1]] = int(newElementEnergies[1])
		xRayData[newElementEnergies[0]][headerEg[8]] = float(newElementEnergies[8])
		xRayData[newElementEnergies[0]][headerEg[9]] = float(newElementEnergies[9])
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[2], 'Energy', float(newElementEnergies[2]))
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[3], 'Energy', float(newElementEnergies[3]))
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[4], 'Energy', float(newElementEnergies[4]))
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[5], 'Energy', float(newElementEnergies[5]))
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[6], 'Energy', float(newElementEnergies[6]))
		bf.setdict2(xRayData[newElementEnergies[0]], headerEg[7], 'Energy', float(newElementEnergies[7]))
		xRayData[newElementEnergies[0]]['KalphaMean']['Energy'] = (xRayData[newElementEnergies[0]]['Kalpha1']['Energy']
			+ xRayData[newElementEnergies[0]]['Kalpha2']['Energy']) / 2
	return xRayData


def loadMaterialData():
	#Name Lattice parameter Lattice file
	data = fg.readLines(FILE_SAMPLE_MATERIALS)
	materialData = {}
	for line in data[0:]:
		lineData = line.split("\t")
		materialData[lineData[0]] = {"a": float(lineData[1])}
		materialData[lineData[0]]["dec"] = fg.dlmread(lineData[2].strip(), "\t", 1)
	return materialData


def importNumericalDataWithHeader(importFile='', delim='\t', withMetaLine=False, skipCols=0,
		fileType=(('Data file', '*.txt'), ('Data file', '*.dat')), diagTitle='Select data file', multiSel='off'):
	data = dict()
	header = 0
	if len(importFile) == 0:
		importFile = fg.requestFiles(fileType, diagTitle, multiSel)
		importFile = importFile[0]
	if len(importFile) > 0:
		# read header data from file
		heading, metaInfo = importFileHeader(importFile, delim, withMetaLine, fileType, diagTitle, multiSel)
		# the first column is text
		if withMetaLine:
			skipRows = 2
		else:
			skipRows = 1
		rawData = fg.dlmread(importFile, delim, skipRows, usedCols=range(skipCols, len(heading)))
		# convert data to dictionary
		data = bf.matrix2dict(rawData, heading[skipCols:])
	return data, metaInfo


def importFileHeader(importFile='', delim='\t', withMetaLine=False,
		fileType=(('Data file', '*.txt'), ('Data file', '*.dat')), diagTitle='Select data file', multiSel='off'):
	if len(importFile) == 0:
		importFile = fg.requestFiles(fileType, diagTitle, multiSel)
		importFile = importFile[0]
	if len(importFile) > 0:
		# read header from file
		if withMetaLine:
			heading = fg.readLines(importFile, 2)
			# first line consists of comments or metadata
			metaInfo = heading[0]
			# second line contains header information
			heading = bf.split(heading[1], delim)
		else:
			heading = fg.readLines(importFile, 1)
			metaInfo = ''
			heading = bf.split(heading[0], '\t')
	return heading, metaInfo


def loadFileP61A(fileP61A=""):
	rawData = 0
	if len(fileP61A) == 0:
		fileP61A = fg.requestFiles((('Data file', '*.dat'),), 'Select P61A data file', 'off')
		fileP61A = fileP61A[0]
	if len(fileP61A) > 0:
		# read header data from file
		heading, metaInfo = importFileHeader(fileP61A, '\t', True, (('Data file', '*.dat'),),
			'Select P61A data file', 'off')
		# the first column is text
		rawData = fg.dlmread(fileP61A, '\t', 2, usedCols=range(1, len(heading)))
	return rawData, heading, metaInfo


def loadFileP61A2(fileP61A=''):
	return importNumericalDataWithHeader(fileP61A, '\t', True, 1,
		(('Data file', '*.dat'),), 'Select P61A data file', 'off')


def sin2PsiHeader():
	fileHead = ['hkl', 'tau_um', 'dStar100_nm', 'dStar100Dev', 's11-s33', 'dev_s11', 's22-s33', 'dev_s22',
		's13', 'dev_s13', 's23', 'dev_s23', 's12', 'dev_s12', 's33', 'dev_s33', 'IBs_keV']
	return bf.stringList2string(fileHead, delim='\t') + bf.newlineChar()


def sin2PsiResults(data):
	text = StringIO()
	text.write('%g\t%.3f\t%.8f\t%.8f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f' %
		(data['hklVal'], data['tauMean'], data['dStar100'], data['dStar100Err'], data['stresses'][0],
		data['accuracy'][0], data['stresses'][1], data['accuracy'][1], data['stresses'][2],
		data['accuracy'][2], data['stresses'][3], data['accuracy'][3], data['stresses'][4],
		data['accuracy'][4], data['stresses'][5], data['accuracy'][5]))
	valuesIb = data['integralWidth']
	for i in range(len(valuesIb)):
		text.write('\t%.6f' % valuesIb[i])
	fg.writeLine(text, '')
	return text.getvalue()


def multiWavelengthResults(data):
	# extract relevant data
	hklList = data['hklList']
	tauMean = data['tauMean']
	dStarVals = data['dStar100']
	dStarErrVals = data['dStar100Err']
	stresses = data['stresses']
	accuracy = data['accuracy']
	integralWidth = data['integralWidth']
	# create tau file content
	text = StringIO()
	for i in range(len(hklList)):
		text.write(sin2PsiResults({'hklVal': hklList[i], 'tauMean': tauMean[i],
			'dStar100': dStarVals[i], 'dStar100Err': dStarErrVals[i], 'stresses': stresses[i],
			'accuracy': accuracy[i], 'integralWidth': integralWidth[i]}))
	return text.getvalue()


def universalplotHeader():
	fileHead = ['tau_um', 's11-s33', 'dev_s11', 's22-s33', 'dev_s22', 's13', 'dev_s13', 's23', 'dev_s23', 'hkl', 'psi']
	return bf.stringList2string(fileHead, delim='\t') + bf.newlineChar()


def universalPlotResults(data):
	# extract relevant data
	tauVals = data['tauVals']
	hklVals = data['hklVals']
	psiVals = data['psiVals']
	stresses = data['stresses']
	accuracy = data['accuracy']
	# create upl file content
	text = StringIO()
	for i in range(len(tauVals)):
		text.write('%.3f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%g\t%.3f' % (tauVals[i],
		stresses[i][0], accuracy[i][0], stresses[i][1], accuracy[i][1], stresses[i][2], accuracy[i][2],
		stresses[i][3], accuracy[i][3], hklVals[i], psiVals[i]))
		fg.writeLine(text, '')
	return text.getvalue()


def universalplotS33Header():
	fileHead = ['tau_um', 'dStar100_nm', 'dStar100Dev', 's33', 'dev_s33', 'hkl']
	return bf.stringList2string(fileHead, delim='\t') + bf.newlineChar()


def universalPlotS33Results(data):
	# create upl file content
	text = StringIO()
	text.write('%.3f\t%.8f\t%.8f\t%.0f\t%.0f\t%g' % (data['tauMean'], data['dStar100'], data['dStar100Err'],
		data['s33'], data['dev_s33'], data['hklVal']))
	fg.writeLine(text, '')
	return text.getvalue()


def multiUniversalPlotS33Results(data):
	# extract relevant data
	hklList = data['hklList']
	tauMean = data['tauMean']
	dStarVals = data['dStar100']
	dStarErrVals = data['dStar100Err']
	stresses = data['s33']
	accuracy = data['dev_s33']
	# create upl s33 file content
	text = StringIO()
	for i in range(len(hklList)):
		text.write(universalPlotS33Results({'hklVal': hklList[i], 'tauMean': tauMean[i],
			'dStar100': dStarVals[i], 'dStar100Err': dStarErrVals[i], 's33': stresses[i],
			'dev_s33': accuracy[i]}))
	return text.getvalue()
