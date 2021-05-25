import numpy as np

import basics.functions as bf
import diffraction.calculations as dc
import filehandling.general as fg
import filehandling.specific as fs
import plotting.specific as ps

# define settings and import the data files
# plotResMwl = False
plotResMwl = True
showErr = True
# showErr = False
fileNames = fg.requestFiles((("Data files", "*.dat"),), "Select P61A data file", "on")
# import all data
allData = dict()
allMetaInfo = dict()
for fileName in fileNames:
	data, metaInfo = fs.loadFileP61A2(fileName)
	allData[fileName] = data
	allMetaInfo[fileName] = metaInfo
# combine all data for one analysis
combinedData = allData[fileNames[0]]
keyList = bf.getKeyList(combinedData)
for fileName in fileNames[1:]:
	for key in keyList:
		combinedData[key] = np.concatenate((combinedData[key], allData[fileName][key]))
# prepare data for multi wavelength and universal plot analysis
a0Val = 0.289
tthVal = 7
# tthVal = 15
inputData = bf.combineDictionaries({'a0Val': a0Val, 'tth': tthVal}, combinedData)
maxPsi = 45
# perform multi wavelength analysis
resDataMwl, plotDataMwl = dc.multiWavelengthAnalysis(inputData, maxPsi)
# plot the results
if plotResMwl:
	ps.plotMultiWavelength(plotDataMwl, showErr)
	#ps.plotIntegralWidth(resDataMwl, showErr)
	ps.plotStrainFreeLatticeSpacing(resDataMwl, showErr)
	ps.plotStresses(resDataMwl, showErr)
# write results to file
fg.export(fs.sin2PsiHeader() + fs.multiWavelengthResults(resDataMwl))

# # import result files of multi wavelength analysis and plot the results
# plotResMwl = True
# fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select MWL result file", "off")
# allData = fg.dlmread(fileNames[0], '\t', 1)
# resDataMwl = {'hklList': allData[:, 0], 's1Dec': None, 'hs2Dec': None,
# 	'tauMean': allData[:, 1], 'dStar100': allData[:, 2], 'stresses': allData[:, 4:15:2],
# 	'accuracy': allData[:, 5:16:2], 'integralWidth': allData[:, 16:22]}
# # plot the results
# if plotResMwl:
# 	ps.plotMultiWavelength(plotDataMwl, showErr)
# 	#ps.plotIntegralWidth(resDataMwl, showErr)
# 	ps.plotStrainFreeLatticeSpacing(resDataMwl, showErr)
# 	ps.plotStresses(resDataMwl, showErr)
