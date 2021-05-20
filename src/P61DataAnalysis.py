import numpy as np
import matplotlib.pyplot as plt

import basics.functions as bf
import basics.fitting as fi
import filehandling.general as fg
import filehandling.specific as fs
import plotting.general as pg
import plotting.specific as ps
import diffraction.conversions as conv
import diffraction.calculations as dc


# load radiation data of materials
xrayData = fs.loadXrayData()

# load diffraction data of materials
# materialData = fs.loadMaterialData()
# print(materialData['Fe']['a'])

# define energy calibration of ED detector: E = a * channel^2 + b * channel + c
# par = np.array([9.45569976605392e-10, 0.00803690845845005, -0.0315429348522124])
# specify used diffraction angle
# tth = 15  # for detector 0
tth = 9.72  # for detector 0
# par = np.array([0, 1, 0])  # if no conversion is needed
# par = np.array([])  # if no conversion is needed
# par = np.array([4.07335065180808e-06, 0.0384242754987209, 8.56168090832901])  # for detector 0 (fluorescence)
# par = np.array([2.10860821094973e-09, 0.0505652975494051, -0.31834244817141])  # for detector 0 (diffraction tth=15)
par = np.array([0, 0.0505710892366929, -0.321876037032894])  # for detector 0 (diffraction tth=15, linear)
# tth = 6.785  # for detector 1
tth = 10.23  # for detector 1
# specify diffraction data of used material
# Fe
hklVals = np.array([110, 200, 211, 220, 310, 222, 321, 400, 411, 420, 332, 422, 431, 521])  # bcc phase
hklVals = np.array([110, 200, 211, 220, 310, 222, 321, 400, 411, 420, 332, 422, 431,
	521, 440, 530, 600, 532, 620, 541, 622, 631, 444, 550, 604])  # bcc phase long
a0Val = 0.28665  # lattice parameter in nm
# FeMnAlNiX
a0Val = 0.28977
# Austenite
hklVals = np.array([111, 200, 220, 311, 222, 400, 331, 420, 422, 333, 440, 531, 600, 620, 533, 622, 444, 551, 640, 642,
	731])  # fcc phase
a0Val = 0.35909
# 316L
a0Val = 0.3595
# IN718
a0Val = 0.3605

# # LaB6
# hklVals = np.array([100, 110, 111, 200, 210, 211, 220, 300, 221, 310, 311, 222,
# 	320, 321, 400, 410, 322, 411, 330, 331, 420, 421, 332, 422, 500, 430, 510,
# 	431, 511, 333, 520, 432, 521, 440, 522, 441, 530, 433, 531, 442, 600, 610])
# hklVals = np.array([100, 110, 111, 200, 210, 211, 220, 300, 310, 311, 222,
#    320, 321, 400, 410, 411, 331, 420, 421, 332, 422, 430, 510])
# a0Val = 0.41569

# Load files with measured data for peakfit
fileNames = fg.requestFiles((("Text files", '*.txt'),), 'Messdateien auswaehlen', 'on')
if len(fileNames) > 0:
	pathName, file = fg.fileparts(fileNames[0])

# accumulate measured spectra
# dataDelim = '\t'
# headerLines = 5
# results, header, files = fs.accumSpectra(pathName + '/outSum_values.txt', fileNames, method='sum')
# results, header, files = fs.accumSpectra(pathName + '/outMean_values.txt', fileNames, method='mean')
# results, header, files = fs.accumSpectra(pathName + '/outMin_values.txt', fileNames, method='min')
# results, header, files = fs.accumSpectra(pathName + '/outMax_values.txt', fileNames, method='max')

# check one example spectrum
[curData, header] = fs.getSpectrumFromFile(fileNames[0])
dVals = conv.aVals2latticeDists2(a0Val, hklVals)
energies = conv.latticeDists2energies(dVals, tth)
adaptTextPosY = False
adaptTextPosY = True
showLines = False
showLines = True
vtext = False
vtext = True
ps.plotPeakPositions1(curData[:, 0], curData[:, 1], np.concatenate(([xrayData['W']['Kalpha1']['Energy'],
	xrayData['W']['Kalpha2']['Energy'], xrayData['W']['Kbeta1']['Energy'],
	xrayData['Pb']['Kalpha1']['Energy'], xrayData['Pb']['Kalpha2']['Energy'],
	xrayData['Pb']['Kbeta1']['Energy']], energies)), np.concatenate((['W-Ka1', 'W-Ka2', 'W-Kb1',
	'Pb-Ka1', 'Pb-Ka2', 'Pb-Kb1'], hklVals)), adaptPosY=adaptTextPosY, plotLines=showLines, verticaltext=vtext)
# ps.plotPeakPositions1(curData[:, 0], curData[:, 1], np.concatenate(([xrayData['W']['Kalpha1']['Energy'],
# 	xrayData['W']['Kalpha2']['Energy'], xrayData['W']['Kbeta1']['Energy'],
# 	xrayData['Pb']['Kalpha1']['Energy'], xrayData['Pb']['Kalpha2']['Energy'],
# 	xrayData['Pb']['Kbeta1']['Energy']], energies)), np.concatenate((['W-Kalpha1', 'W-Kalpha2', 'W-Kbeta1',
# 	'Pb-Kalpha1', 'Pb-Kalpha2', 'Pb-Kbeta1'], hklVals)), adaptPosY=adaptTextPosY, plotLines=showLines)
# # plot with own channel to energy conversion
ps.plotPeakPositions1(conv.channels2energies(np.array(range(len(curData[:, 0]))), par), curData[:, 1],
	np.concatenate(([xrayData['W']['Kalpha1']['Energy'], xrayData['W']['Kalpha2']['Energy'],
	xrayData['W']['Kbeta1']['Energy'], xrayData['Pb']['Kalpha1']['Energy'], xrayData['Pb']['Kalpha2']['Energy'],
	xrayData['Pb']['Kbeta1']['Energy']], energies)), np.concatenate((['W-Ka1', 'W-Ka2', 'W-Kb1',
	'Pb-Ka1', 'Pb-Ka2', 'Pb-Kb1'], hklVals)), adaptPosY=adaptTextPosY, plotLines=showLines, verticaltext=vtext)
# ps.plotPeakPositions1(dc.channels2energies(np.array(range(len(curData[:, 0]))), par), curData[:, 1],
# 	np.concatenate(([xrayData['W']['Kalpha1']['Energy'], xrayData['W']['Kalpha2']['Energy'],
# 	xrayData['W']['Kbeta1']['Energy'], xrayData['Pb']['Kalpha1']['Energy'], xrayData['Pb']['Kalpha2']['Energy'],
# 	xrayData['Pb']['Kbeta1']['Energy']], energies)), np.concatenate((['W-Kalpha1', 'W-Kalpha2', 'W-Kbeta1',
# 	'Pb-Kalpha1', 'Pb-Kalpha2', 'Pb-Kbeta1'], hklVals)))

# plot only intensities over channels
ps.plotPeakPositions1(np.array(range(len(curData[:, 0]))), curData[:, 1])

# specify peak ranges in spectrum
# Fe peak positions
# transmission, tth = 6.785: 110, Wka2, Wka1, Wkb1, Pbka2/200, Pbka1, Pbkb1, 211, 220, 310, 222, 321
peakZones = np.array([[50.1, 53.3], [57.3, 58.5], [58.6, 60], [65.8, 68.2], [71.8, 74], [74, 76],
	[83.5, 85.9], [88.5, 90.5], [102.5, 104.4], [114, 117], [125.6, 127.4], [135.5, 138]])
# transmission, tth = 15: ???, ???, 200, 211, 220, 310, 321, 411, 332, 422, 521, 440, 530, 600, 532
peakZones = np.array([[450, 510], [520, 575], [640, 685], [780, 845], [900, 965], [1010, 1075], [1205, 1270],
	[1360, 1425], [1520, 1565], [1585, 1640], [1775, 1830], [1835, 1885], [1890, 1945], [1950, 1995], [2000, 2055]])
# transmission, tth = 6.785, IN718: 111, 200/Wka2/Wka1, Wkb1, Pbka2, Pbka1, 220, Pbkb1, 311, 222, 400
peakZones = np.array([[49.5, 51.5], [56, 60.5], [65.8, 68.2], [71.7, 73.9], [73.9, 76.2],
	[81.2, 83.4], [83.5, 85.9], [94.4, 98.4], [99.2, 102.5], [114.3, 118.5]])
# reflection with carbides, tth = 6.785: 110, Wka2, Wka1, Wkb1, Pbka2/200, Pbka1, Pbkb1, 211, 220, 310, 222, 321,
# 400, 411, 420 (bad shape), 332, 422, 431, 521
peakZones = np.array([[50.1, 53.3], [57.3, 58.5], [58.6, 60], [65.8, 68.2], [71, 74], [74, 76],
	[83.5, 85.9], [88, 91.5], [101.5, 105], [113.7, 117.5], [125.4, 128], [134.5, 139], [144.4, 148],
	[153, 157], [160.5, 166], [169.5, 173.2], [177.5, 181], [184.6, 187.68], [197.2, 204.2]])
# reflection without carbides, tth = 6.785: 110, Wka2, Wka1, Wkb1, Pbka2/200, Pbka1, Pbkb1, 211, 220, 310, 222, 321,
# 400, 411, 420 (bad shape), 332, 422, 431, 521
peakZones = np.array([[50.1, 53.3], [57.3, 58.5], [58.6, 60], [65.8, 68.2], [71, 74], [74, 76],
	[83.5, 86], [88, 91.5], [101.5, 105.5], [113.7, 118], [125.4, 128.4], [134.5, 139.3], [144.4, 148],
	[152.7, 158], [160.5, 166.8], [169.5, 173.8], [177.1, 181.3], [183.2, 190], [197.2, 204.2]])
# reflection, tth = 9.72: 110, 200, Wka2, Wka1, 211, Wkb, 220+Pbka2, Pbka1, 310, Pbkb, 321, 411, 420, 332, 422, 431, 521
peakZones = np.array([[686, 752], [980, 1050], [1128, 1163], [1163, 1205], [1210, 1280], [1300, 1356], [1410, 1467],
	[1467, 1518], [1560, 1635], [1645, 1712], [1865, 1935], [2124, 2184], [2236, 2295], [2328, 2415], [2450, 2518],
	[2544, 2625], [2738, 2813]])
# reflection, tth = 9.72: 110, 200, Wka2, Wka1, 211, Wkb, Pbka1, 310, Pbkb, 321, 411, 420, 332, 422, 431, 521
peakZones = np.array([[33.26, 38.67], [47.82, 53.69], [56.52, 58.49], [58.49, 60.21], [60.46, 64.11], [65.77, 68.35],
	[73.71, 76.19], [79.02, 81.65], [83.22, 86.15], [93.84, 96.88], [106.6, 109.5], [112.6, 115.5], [117.9, 121.3], [123.2, 126.5],
	[128.1, 132], [138, 141.3]])
# reflection, tth = 10.23: 110, 200, Wka1+211, Wkb, 220, Pbka2, Pbka1, 310, 321, 411, 420, 332, 422, 431, 521
peakZones = np.array([[33.48, 35.13], [47.67, 49.38], [58.67, 60.23], [66.13, 67.97], [67.88, 69.58], [71.92, 73.53],
	[74.03, 75.83], [75.83, 77.58], [89.83, 91.63], [102, 103.7], [107.5, 109.4], [112.8, 114.6], [117.5, 120.1],
	[122.4, 125], [131.6, 134.2]])
# reflection, tth = 6.225: 110, Wka2, Wka1, Wkb, Pbka2, Pbka1, 200, Pbkb, 211, 220, 310, 222, 321, 411
peakZones = np.array([[54.58, 57.52], [57.52, 58.48], [58.48, 60.13], [65.97, 68.08], [71.67, 73.83], [73.83, 76.42],
	[78.22, 81.3], [83.38, 86.17], [96, 99.08], [111.2, 114.1], [124.6, 127.3], [136.9, 139.2], [147.7, 150.5],
	[167.7, 170.3]])
# transmission, IN718, tth = 10.23: 111, 200, 220, Wka2, Wka1, 311, Pbka2, Pbka1, 400, 422, 333, 440, 620, 642, 553
peakZones = np.array([[32.6, 34.08], [37.73, 39.23], [53.13, 55.58], [57.13, 58.63], [58.58, 60.1], [62.38, 65.1],
	[71.63, 73.83], [73.83, 75.88], [75.83, 78.22], [93, 95.42], [98.47, 101.5], [107.8, 110.2], [120.5, 123.2],
	[142.2, 145.6], [146, 149.6]])
# FeMnAlNiX peak positions
# reflection, tth = 6.785: 110, Wka2, Wka1, Wkb1, Pbka2/200, Pbka1, Pbkb1, 211, 220, 310, 222, 321,
# 400, 411, 420 (bad shape), 332, 422, 431, 521
peakZones = np.array([[50.3, 52.2], [57.3, 58.5], [58.6, 60], [65.8, 68.2], [71, 74], [74, 76],
	[83.7, 85.9], [88, 91.5], [101.5, 105], [113.7, 117.5], [125.4, 128], [134.5, 139], [144.4, 148],
	[153, 157], [160.5, 166], [169.5, 173.2], [177.5, 181], [184.6, 187.68], [197.2, 204.2]])
# reflection, tth = 6.785, Schott sample: austenite 311, ferrite 321
peakZones = np.array([[94, 100], [132.5, 141.5]])
# peakZones = np.array([[1421, 1471], [1468, 1521], [1646, 1724]])  # Channels of Pb fluorescence of Detector 1

# specify peak fit settings
peakfitSettings = fi.initPeakfitSettings()
# peakfitSettings["peakFunction"] = 'PseudoVoigt'  # 'Gauss', 'Pearson7', 'Lorentz', 'PseudoVoigt'
# peakfitSettings["withGraphic"] = 1
peakfitSettings["withGraphic"] = 0
# peakfitSettings["onlyValues"] = 1 # for first detector for calibration
peakfitSettings["peakZones"] = peakZones
# import settings from file
filesSettings = fg.requestFiles((("Text files", "*.txt"), ("All files", "*.*")), "Select peak fit settings file", "off")
peakfitSettings = fi.initPeakfitSettings(filesSettings[0])
# export settings to file
fg.saveSettings(peakfitSettings)

# perform peak fit
# [analysisData, measureData, headers, peakResFiles] = fi.evalPeakFiles(peakfitSettings, (fileNames[2],), pathName + '/results_')
[analysisData, measureData, headers, peakResFiles] = fi.evalPeakFiles(peakfitSettings, fileNames, pathName + '/results_')

# plot result of peakfit in complete spectrum
selFileNames = fg.requestFiles((("Text files", '*.txt'),), 'Messdateien auswaehlen', 'on')
# selFileNames = (fileNames[22],)
# selFileNames = fileNames[0:5:]
# selFileNames = fileNames
for i in range(len(selFileNames)):
	[data, heading] = fs.getSpectrumFromFile(selFileNames[i])
	if peakfitSettings["onlyValues"]:
		data[:, 0] = np.array(range(1, len(data[:, 0]) + 1))
	pathName, file = fg.fileparts(selFileNames[i])
	peakData = fg.dlmread(pathName + '/results_' + file, '\t', 1)
	pg.figure()
	plt.plot(data[:, 0], data[:, 1], 'k.-')
	for j in range(peakData.shape[0]):
		xVals = data[(data[:, 0] >= peakData[j, 0]) & (data[:, 0] <= peakData[j, 1]), 0]
		plt.plot(xVals, fi.gaussfunction(peakData[j, [2, 3, 4, 8, 9]], xVals), 'r.:')
	plt.grid("on")
	plt.tight_layout()  # layout without overlapping
	plt.show()

# remove positions of positions files where spectra are missing
headLine = bf.stringList2string(['om', 'eta', 'tth0', 'tth1', 'psi0', 'psi1', 'phi', 'x', 'y', 'z'], '\t')
checkAxes = {'_eu.x_': 7, '_eu.chi_': 5}
fileNamesAxes = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select positions files", "on")
peakResFiles = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select peakfit files", "on")
for axisFile in fileNamesAxes:
	curPosData = fg.dlmread(axisFile, '\t', 1)
	lineExists = np.zeros(curPosData.shape[0])
	for i in range(curPosData.shape[0]):
		selList = list(checkAxes.keys())
		curAxesVals = curPosData[i, list(checkAxes.values())]
		[selList.insert(j * 2 + 1, str(curAxesVals[j])) for j in range(len(curAxesVals))]
		selStr = bf.stringList2string(selList)
		lineExists[i] = sum(1 for f in peakResFiles if selStr in f)
	fg.dlmwrite(bf.replace(axisFile, 'motpos', 'motpos_sel'), curPosData[lineExists > 0, :],
		head=headLine, format='%f')


# specify final file settings
finalFileSettings = fs.initFinalFileSettings()
# finalFileSettings["type"] = 'P61_0'
# finalFileSettings["type"] = 'P61_1'
# finalFileSettings["unused"] = np.array([1, 2])  # detector 0, tth=15, bcc, starts with one
# finalFileSettings["unused"] = np.array([3, 4, 6, 7, 8, 10])  # detector 0, tth=9.72, bcc
# finalFileSettings["unused"] = np.array([3, 4, 6, 7])  # detector 1, tth=10.23, bcc
# finalFileSettings["unused"] = np.array([2, 3, 4, 5, 6, 7])  # detector 1, tth=6.785, bcc
# finalFileSettings["unused"] = np.array([3, 4, 5, 7])  # detector 1, tth=6.785, fcc
# finalFileSettings["unused"] = np.array([1])
# finalFileSettings["unused"] = np.array([])
# finalFileSettings["energy"] = 1 # not used now and actually not needed!!!
# finalFileSettings["pars"] = np.array([])
# finalFileSettings["pars"] = par
# import settings from file
filesSettings = fg.requestFiles((("Text files", "*.txt"), ("All files", "*.*")), "Select final file settings file", "off")
finalFileSettings = fs.initFinalFileSettings(filesSettings[0])
# export settings to file
fg.saveSettings(finalFileSettings)

# create full measurement result file
fileNamesAxes = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select positions file", "off")
peakResFiles = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select peakfit files", "on")
data, resFile = fs.createPeakMeasurementFile(finalFileSettings, fileNamesAxes, peakResFiles)

# create full measurement result files for transmission for line scans (automatic)
# posFilePrefix = "motpos_sel_"
posFilePrefix = "motpos_"
# scanAxis = "_eu.x_"
# scanAxis = "_eu.y_"
scanAxis = "_eu.z_"
fileNamesAxes = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select positions files", "on")
peakResFiles = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select peakfit files", "on")
for axisFile in fileNamesAxes:
	pathName, file = fg.fileparts(axisFile)
	selStr = scanAxis + bf.replace(bf.replace(file, posFilePrefix), ".txt") + "_"
	curPeakResFiles = tuple(f for f in peakResFiles if selStr in f)
	data, resFile = fs.createPeakMeasurementFile(finalFileSettings, (axisFile,), curPeakResFiles)

# create full measurement result files for transmission for multidimensional scans (automatic)
sepStr = '_'
# posFilePrefix = "motpos_sel_"
posFilePrefix = "motpos_"
scanAxis1 = 'x'
scanAxis2 = 'y'
scanAxis3 = 'z'
# scanAxis1 = None
# scanAxis2 = None
# scanAxis3 = None
fileNamesAxes = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select positions files", "on")
peakResFiles = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select peakfit files", "on")
selStr1 = None
selStr2 = None
selStr3 = None
for axisFile in fileNamesAxes:
	pathName, file = fg.fileparts(axisFile)
	fileParts = file.replace('.txt', '').split(sepStr)
	if scanAxis1 is not None:
		selStr1 = [fileParts[i] + sepStr + fileParts[i+1] for i in range(len(fileParts)) if scanAxis1 in fileParts[i]]
	if scanAxis2 is not None:
		selStr2 = [fileParts[i] + sepStr + fileParts[i+1] for i in range(len(fileParts)) if scanAxis2 in fileParts[i]]
	if scanAxis3 is not None:
		selStr3 = [fileParts[i] + sepStr + fileParts[i+1] for i in range(len(fileParts)) if scanAxis3 in fileParts[i]]
	curPeakResFiles = tuple(f for f in peakResFiles if (selStr1 is None or len(selStr1) == 0 or selStr1[0] in f) and
		(selStr2 is None or len(selStr2) == 0 or selStr2[0] in f) and
		(selStr3 is None or len(selStr3) == 0 or selStr3[0] in f))
	if len(curPeakResFiles) > 0:
		data, resFile = fs.createPeakMeasurementFile(finalFileSettings, (axisFile,), curPeakResFiles)

# correct data according reference data
fileNamesRef = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select reference final file", "off")
fileNamesSample = fg.requestFiles((("Text files", "*.txt"), ("all files", "*.*")), "Select sample final files", "on")
dataRef = fg.dlmread(fileNamesRef[0], '\t', 1)
fileHead = bf.stringList2string(['LNr', 'dVal[nm]', 'dVar[nm]', 'Iint', 'Integralb', 'tth', 'phiP', 'psiP', 'etaP', 'Ringstr', 'RT', 'DT',
	'xdiff', 'ydiff', 'zdiff', 'motor1', 'motor2', 'motor3', 'temp1', 'temp2', 'wavelength', 'deltatime'], '\t')
for curFile in fileNamesSample:
	corrFile = curFile.replace('.txt', '_corr.txt')
	dataCur = fg.dlmread(curFile, '\t', 1)
	for i in range(dataCur.shape[0]):
		curPeakNum = dataCur[i, 0]
		curPhi = dataCur[i, 6]
		curPsi = dataCur[i, 7]
		dataCur[i, 1] -= (dataRef[(dataRef[:, 0] == curPeakNum) & (dataRef[:, 6] == curPhi) & (dataRef[:, 7] == curPsi), 1] -
			dataRef[(dataRef[:, 0] == curPeakNum) & (dataRef[:, 6] == 0) & (dataRef[:, 7] == 0), 1])
	fg.dlmwrite(corrFile, dataCur, '\t', head=fileHead)

# specify stress analysis settings
stressAnalysisSettings = dc.initStressAnalysisSettings()
stressAnalysisSettings["showPlots"] = True
stressAnalysisSettings["showDeviation"] = False
stressAnalysisSettings["method"] = 'ED'
stressAnalysisSettings["decList"] = np.array([[110, -1.27e-06, 5.8e-06], [211, -1.27e-06, 5.8e-06], [220, -1.27e-06, 5.8e-06],
	[310, -1.67e-06, 7.02e-06], [222, -1.05e-06, 5.17e-06], [321, -1.27e-06, 5.8e-06]])  # transmission, tth = 6.785
stressAnalysisSettings["decList"] = np.array([[200, -1.90e-06, 7.70e-06], [211, -1.27e-06, 5.8e-06],
	[220, -1.27e-06, 5.8e-06], [310, -1.67e-06, 7.02e-06], [321, -1.27e-06, 5.8e-06], [411, -1.64e-06, 6.93e-06],
	[332, -1.10e-06, 5.30e-06], [422, -1.27e-06, 5.8e-06], [521, -1.537e-06, 6.613e-06], [440, -1.27e-06, 5.8e-06],
	[530, -1.41e-06, 6.22e-06], [600, -1.90e-06, 7.70e-06], [532, -1.27e-06, 5.8e-06]])  # transmission, tth = 15 (detector 0)
stressAnalysisSettings["decList"] = np.array([[110, -1.27e-06, 5.8e-06], [211, -1.27e-06, 5.8e-06],
	[220, -1.27e-06, 5.8e-06], [310, -1.67e-06, 7.02e-06], [222, -1.05e-06, 5.17e-06], [321, -1.27e-06, 5.8e-06],
	[400, -1.90e-06, 7.70e-06], [411, -1.64e-06, 6.93e-06], [420, -1.50e-06, 6.49e-06], [332, -1.10e-06, 5.30e-06],
	[422, -1.27e-06, 5.8e-06], [431, -1.27e-06, 5.8e-06], [521, -1.537e-06, 6.613e-06]])  # reflection, tth = 6.785
stressAnalysisSettings["decList"] = np.array([[311, -1.8e-06, 7.52e-06]])  # reflection, Austenite, tth = 6.785
stressAnalysisSettings["decList"] = np.array([[110, -1.27e-06, 5.8e-06], [200, -1.90e-06, 7.70e-06],
	[211, -1.27e-06, 5.8e-06], [310, -1.67e-06, 7.02e-06], [321, -1.27e-06, 5.8e-06], [411, -1.64e-06, 6.93e-06],
	[420, -1.50e-06, 6.49e-06], [332, -1.10e-06, 5.30e-06], [422, -1.27e-06, 5.8e-06], [431, -1.27e-06, 5.8e-06],
	[521, -1.537e-06, 6.613e-06]])  # reflection, tth = 9.72 (detector 0)
stressAnalysisSettings["decList"] = np.array([[110, -1.27e-06, 5.8e-06], [200, -1.90e-06, 7.70e-06],
	[220, -1.27e-06, 5.8e-06], [310, -1.67e-06, 7.02e-06], [321, -1.27e-06, 5.8e-06], [411, -1.64e-06, 6.93e-06],
	[420, -1.50e-06, 6.49e-06], [332, -1.10e-06, 5.30e-06], [422, -1.27e-06, 5.8e-06], [431, -1.27e-06, 5.8e-06],
	[521, -1.537e-06, 6.613e-06]])  # reflection, tth = 10.23 (detector 1)
# calculate DEC values for 521
# s1Calc = np.interp(0.43, [0, 0.27, 0.75, 1], [-1.90e-06, -1.67e-06, -1.27e-06, -1.05e-06])
# hs2Calc = np.interp(0.43, [0, 0.27, 0.75, 1], [7.70e-06, 7.02e-06, 5.8e-06, 5.17e-06])
stressAnalysisSettings["material"] = 'Fe'
stressAnalysisSettings["a0Val"] = a0Val
stressAnalysisSettings["maxPsi"] = 45
# stressAnalysisSettings["maxPsi"] = 60
# stressAnalysisSettings["maxPsi"] = 90
# import settings from file
files = fg.requestFiles((("Text files", "*.txt"), ("All files", "*.*")), "Select stress analysis settings file", "off")
stressAnalysisSettings = dc.initStressAnalysisSettings(files[0])
# export settings to file
fg.saveSettings(stressAnalysisSettings)

'''
# perform sin2psi and universalplot analysis
# fileNames = fg.requestFiles((("text files","*.txt"),("all files","*.*")), "Select final file", "off")
# resFile = fileNames[0]
stressResFile, results, dStarVals0, dStarVals90, regVals, regVals0, regVals90, maxPsi, sMain, d0 = \
	dc.createSin2PsiFile(stressAnalysisSettings, resFile)
# if stressAnalysisSettings["a0Val"] is not None and stressAnalysisSettings["a0Val"] > 0:
# 	resFileS33, resData = dc.createS33File(stressAnalysisSettings["a0Val"], stressAnalysisSettings["decList"], resFile)
# 	sorted, s33Data = dc.createUniversalplotFile(stressAnalysisSettings["method"], stressAnalysisSettings["decList"],
# 		stressAnalysisSettings["material"], stressAnalysisSettings["a0Val"], resFile)
'''

# perform sin2psi and universal plot analysis based on dictionaries
# plotResMwl = False
# plotResUvp = False
plotResMwl = True
plotResUvp = True
showErr = False
# showErr = True
checkUvpVals = True
# import all data
# fileNames = fg.requestFiles((("text files", "*.txt"), ("all files", "*.*")), "Select final file", "off")
# resFile = fileNames[0]
allData = fg.dlmread(resFile, '\t', 1)
# prepare data for multi wavelength and universal plot analysis
peakNum = int(max(allData[:, 0]))
minCount = 100000
maxCount = 0
for i in range(peakNum):
	minCount = bf.min(minCount, sum(allData[:, 0] == i + 1))[0]
	maxCount = bf.max(maxCount, sum(allData[:, 0] == i + 1))[0]
if minCount == maxCount:
	combinedData = {'x': allData[allData[:, 0] == 1, 12], 'y': allData[allData[:, 0] == 1, 13],
		'z': allData[allData[:, 0] == 1, 14], 'eta': allData[allData[:, 0] == 1, 8], 'psi': allData[allData[:, 0] == 1, 7],
		'phi': allData[allData[:, 0] == 1, 6], 'pv_amplitude': 0}
	for i in range(peakNum):
		combinedData['pv' + str(i) + '_amplitude'] = allData[allData[:, 0] == i + 1, 3]
		combinedData['pv' + str(i) + '_amplitude_err'] = np.zeros(maxCount)
		combinedData['pv' + str(i) + '_fraction'] = np.ones(maxCount)
		combinedData['pv' + str(i) + '_fraction_err'] = np.zeros(maxCount)
		combinedData['pv' + str(i) + '_sigma'] = allData[allData[:, 0] == i + 1, 4]
		combinedData['pv' + str(i) + '_sigma_err'] = np.zeros(maxCount)
		dVals = allData[allData[:, 0] == i + 1, 1]
		energies = conv.latticeDists2energies(dVals, allData[0, 5])
		energyDiffs = conv.latticeDists2energies(dVals + allData[allData[:, 0] == i + 1, 2], allData[0, 5])
		combinedData['pv' + str(i) + '_center'] = energies
		combinedData['pv' + str(i) + '_center_err'] = np.abs(energyDiffs - energies)
		mueVals = dc.calcMue(energies, stressAnalysisSettings["material"])
		tauVals = dc.calcTauPsiEta(mueVals, allData[0, 5], allData[allData[:, 0] == 1, 7], allData[allData[:, 0] == 1, 8])
		combinedData['pv' + str(i) + '_depth'] = tauVals
		combinedData['pv' + str(i) + '_dspac'] = dVals * 10
		combinedData['pv' + str(i) + '_lp'] = conv.latticeDists2aVals2(allData[allData[:, 0] == i + 1, 1],
			stressAnalysisSettings["decList"][i, 0])
		h, k, l = conv.splitHkl(stressAnalysisSettings["decList"][i, 0])
		combinedData['pv' + str(i) + '_h'] = np.ones(maxCount) * int(h)
		combinedData['pv' + str(i) + '_k'] = np.ones(maxCount) * int(k)
		combinedData['pv' + str(i) + '_l'] = np.ones(maxCount) * int(l)
		combinedData['pv' + str(i) + '_s1'] = np.ones(maxCount) * stressAnalysisSettings["decList"][i, 1]
		combinedData['pv' + str(i) + '_hs2'] = np.ones(maxCount) * stressAnalysisSettings["decList"][i, 2]
		combinedData['pv' + str(i) + '_3gamma'] = np.ones(maxCount) * dc.calc3Gamma2(stressAnalysisSettings["decList"][i, 0])
	inputData = bf.combineDictionaries({'a0Val': a0Val, 'tth': allData[0, 5]}, combinedData)
else:
	print('There are missing values, handling is yet not implemented')

# perform multi wavelength analysis
maxPsi = 45
# maxPsi = 60
# maxPsi = 90
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

# combine stress results of different points (mappings, time and line scans)
# baseName = 'sin2alpha_stresses_'
baseName = 'stresses_'
# define relevant axis
relevantAxis = ['x', 'y', 'z']
# define number of evaluated stresses
# stressCount = 5
stressCount = 6
# stressCount = 7
# stressCount = 8
# stressCount = 11
# stressCount = 13
# stressCount = 14
fileNames = fg.requestFiles((("Text files", '*.txt'),), 'Select single stress files', 'on')
# define variables for relevant values
posVals = np.zeros((len(fileNames), len(relevantAxis)))  # position values
aStarData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
stressData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
stressVarData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
stressShearData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
stressShearVarData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
widthData = np.zeros((len(fileNames), stressCount))  # depths steps x peak number
# iterate over all files
for i in range(len(fileNames)):
	pathName, file = fg.fileparts(fileNames[i])
	fileParts = file.replace('.txt', '').split(sepStr)
	for j in range(len(relevantAxis)):
		posVals[i, j] = float([fileParts[i + 1] for i in range(len(fileParts)) if relevantAxis[j] == fileParts[i]][0])
	# import the data of current file
	curStressData = fg.dlmread(fileNames[i], '\t', 1)
	aStarData[i, :] = curStressData[:, 2]
	stressData[i, :] = curStressData[:, 4]
	stressVarData[i, :] = curStressData[:, 5]
	stressShearData[i, :] = curStressData[:, 8]
	stressShearVarData[i, :] = curStressData[:, 9]
	widthData[i, :] = curStressData[:, 16]
for i in range(curStressData.shape[0]):
	fg.dlmwrite(pathName + '/' + baseName + str(int(curStressData[i, 0])) + '.txt', np.concatenate((posVals,
		np.transpose([stressData[:, i], stressVarData[:, i], widthData[:, i]])), 1), format='%f')
