import numpy as np

import BasicFunctions.generalFunctions as gf
import Diffraction.diffractionCalculations as dc
import Diffraction.conversions as conv
import FileHandling.generalFileHandling as gfh
import FileHandling.specificFileHandling as sfh
import Plotting.specificPlotting as spl

# test hkl generator and energy calculation of material
phase = 'bcc'
aVal = 0.28665
tth = 10
energyRange = [5, 300]
hklData = dc.hklGenerator(phase, aVal, tth, energyRange)
print(hklData[0:6, :])

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

# second version can be also applied to first test files (modified header) using new functions working with dictionaries
# plotResMwl = False
# plotResUvp = False
plotResMwl = True
plotResUvp = True
showErr = True
checkUvpVals = True
fileNames = gfh.requestFiles((("Data files", "*.dat"),), "Select P61A data file", "on")
# import all data
allData = dict()
allMetaInfo = dict()
for fileName in fileNames:
	data, metaInfo = sfh.loadFileP61A2(fileName)
	allData[fileName] = data
	allMetaInfo[fileName] = metaInfo
# combine all data for one analysis
combinedData = allData[fileNames[0]]
keyList = gf.getKeyList(combinedData)
for fileName in fileNames[1:]:
	for key in keyList:
		combinedData[key] = np.concatenate((combinedData[key], allData[fileName][key]))
# prepare data for multi wavelength and universal plot analysis
a0Val = 0.289
tthVal = 7
inputData = gf.combineDictionaries({'a0Val': a0Val, 'tth': tthVal}, combinedData)
# perform multi wavelength analysis
maxPsi = 45
resDataMwl, plotDataMwl = dc.multiWavelengthAnalysis(inputData, maxPsi)
# plot the results
if plotResMwl:
	spl.plotMultiWavelength(plotDataMwl, showErr)
	#spl.plotIntegralWidth(resDataMwl, showErr)
	spl.plotStrainFreeLatticeSpacing(resDataMwl, showErr)
	spl.plotStresses(resDataMwl, showErr)
# write results to file
gfh.export(sfh.sin2PsiHeader() + sfh.multiWavelengthResults(resDataMwl))
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
	spl.plotUniversalPlot(resDataUvp, showErr)
	#spl.plotMultiUniversalPlot(resDataUvp, showErr)
	spl.plotStrainFreeLatticeSpacing(resDataS33, showErr)
	spl.plotStresses(resDataS33, showErr)
# write results to files
gfh.export(sfh.universalplotHeader() + sfh.universalPlotResults(resDataUvp))
gfh.export(sfh.universalplotS33Header() + sfh.multiUniversalPlotS33Results(resDataS33))
