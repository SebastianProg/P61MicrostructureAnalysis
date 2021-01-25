# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:53:13 2019

@author: Sebastian
"""

from io import StringIO  # for handling unicode strings

import basics.functions as bf
import filehandling.general as fg


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
