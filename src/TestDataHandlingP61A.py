import numpy as np

import basics.functions as bf
import diffraction.calculations as dc
import diffraction.conversions as conv
import filehandling.general as fg
import filehandling.specific as fs
import plotting.specific as ps

# test hkl generator and energy calculation of material
bf.init()
phase = 'bcc'
aVal = 0.28665
tth = 7
energyRange = [5, 300]
hklData = dc.hklGenerator(phase, aVal, tth, energyRange)
print(hklData[0:16, :])

# peak positions of LaB6
a0ValLab6 = 0.41569
phaseTextsLab6 = 'LaB6'
hklValsLab6 = np.transpose(np.array([100, 110, 111, 200, 210, 211, 220, 300, 310, 311, 222, 320, 321, 400, 410, 411,
	331, 420, 421, 332, 422, 430, 510]))
dVals = conv.aVals2latticeDists2(a0ValLab6, hklValsLab6)
energiesDet2 = conv.latticeDists2energies(dVals, 11.5)
energiesDet1 = conv.latticeDists2energies(dVals, 15)
for i in range(len(hklValsLab6)):
	print('%d\t%.3f\t%.3f' % (hklValsLab6[i], energiesDet2[i], energiesDet1[i]))
dc.calc3Gamma2(hklValsLab6)

# define settings and import the data files
# plotResMwl = False
# plotResUvp = False
plotResMwl = True
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
# perform multi wavelength analysis
maxPsi = 45
resDataMwl, plotDataMwl = dc.multiWavelengthAnalysis(inputData, maxPsi)
# plot the results
if plotResMwl:
	ps.plotMultiWavelength(plotDataMwl, showErr)
	#ps.plotIntegralWidth(resDataMwl, showErr)
	ps.plotStrainFreeLatticeSpacing(resDataMwl, showErr)
	ps.plotStresses(resDataMwl, showErr)
# write results to file
fg.export(fs.sin2PsiHeader() + fs.multiWavelengthResults(resDataMwl))
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

# import result files of multi wavelength analysis and plot the results
plotResMwl = True
fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select MWL result file", "off")
allData = fg.dlmread(fileNames[0], '\t', 1)
resDataMwl = {'hklList': allData[:, 0], 's1Dec': None, 'hs2Dec': None,
	'tauMean': allData[:, 1], 'dStar100': allData[:, 2], 'stresses': allData[:, 4:15:2],
	'accuracy': allData[:, 5:16:2], 'integralWidth': allData[:, 16:22]}
# plot the results
if plotResMwl:
	ps.plotMultiWavelength(plotDataMwl, showErr)
	#ps.plotIntegralWidth(resDataMwl, showErr)
	ps.plotStrainFreeLatticeSpacing(resDataMwl, showErr)
	ps.plotStresses(resDataMwl, showErr)

# import result files of universal plot analysis and plot the results
plotResUvp = True
fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select UVP result file", "off")
allData = fg.dlmread(fileNames[0], '\t', 1)
resDataUvp = {'tauVals': allData[:, 0], 'stresses': allData[:, 1:9:2], 'accuracy': allData[:, 2:10:2],
	'hklList': allData[:, 9], 'psiVals': allData[:, 10]}
fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select UVP-s33 result file", "off")
allData = fg.dlmread(fileNames[0], '\t', 1)
resDataS33 = {'tauMean': allData[:, 0], 'aStarVals': allData[:, 1], 'dStar100': allData[:, 2],
	'stresses': allData[:, 3], 'accuracy': allData[:, 4], 'hklList': allData[:, 5]}
# plot the results
if plotResUvp:
	ps.plotUniversalPlot(resDataUvp, showErr)
	#ps.plotMultiUniversalPlot(resDataUvp, showErr)
	ps.plotStrainFreeLatticeSpacing(resDataS33, showErr)
	ps.plotStresses(resDataS33, showErr)
