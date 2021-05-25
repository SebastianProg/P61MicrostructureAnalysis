import numpy as np

import basics.functions as bf
import diffraction.calculations as dc
import filehandling.general as fg
import filehandling.specific as fs
import plotting.specific as ps

# define settings and import the data files
# plotResUvp = False
plotResUvp = True
showErr = True
# showErr = False
checkUvpVals = True
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
# perform universal plot
if checkUvpVals:
	minDistPsiStar = 0.15
	minValPsiNormal = 0.08
	minValPsiShear = 0.8
else:
	minDistPsiStar = None
	minValPsiNormal = None
	minValPsiShear = None
resDataUvp, resDataS33 = dc.multiUniversalPlotAnalysis(inputData, maxPsi, minDistPsiStar,
	minValPsiNormal, minValPsiShear)
# plot the results
if plotResUvp:
	ps.plotUniversalPlot(resDataUvp, showErr)
	#ps.plotMultiUniversalPlot(resDataUvp, showErr)
	ps.plotStrainFreeLatticeSpacing(resDataS33, showErr)
	ps.plotStresses(resDataS33, showErr)
# write results to files
fg.export(fs.universalplotHeader() + fs.universalPlotResults(resDataUvp))
fg.export(fs.universalplotS33Header() + fs.multiUniversalPlotS33Results(resDataS33))

# # import result files of universal plot analysis and plot the results
# plotResUvp = True
# fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select UVP result file", "off")
# allData = fg.dlmread(fileNames[0], '\t', 1)
# resDataUvp = {'tauVals': allData[:, 0], 'stresses': allData[:, 1:9:2], 'accuracy': allData[:, 2:10:2],
# 	'hklList': allData[:, 9], 'psiVals': allData[:, 10]}
# fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select UVP-s33 result file", "off")
# allData = fg.dlmread(fileNames[0], '\t', 1)
# resDataS33 = {'tauMean': allData[:, 0], 'aStarVals': allData[:, 1], 'dStar100': allData[:, 2],
# 	'stresses': allData[:, 3], 'accuracy': allData[:, 4], 'hklList': allData[:, 5]}
# # plot the results
# if plotResUvp:
# 	ps.plotUniversalPlot(resDataUvp, showErr)
# 	#ps.plotMultiUniversalPlot(resDataUvp, showErr)
# 	ps.plotStrainFreeLatticeSpacing(resDataS33, showErr)
# 	ps.plotStresses(resDataS33, showErr)
