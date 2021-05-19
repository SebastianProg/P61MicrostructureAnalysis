import numpy as np
import matplotlib.pyplot as plt

import basics.functions as bf
import basics.calculations as bc
import plotting.general as pg


def stereoProjPlot(phi, psi, int=None, lineSpec='', titleText=None, maxPsi=90, fig=None, partPlot=False, thetaStep=None,
		rStep=None, rMax=None, cbarOri='vertical'):
	# only allow defined psi values
	phi = phi[(psi >= -maxPsi) & (psi <= maxPsi)]
	if int is not None:
		int = int[(psi >= 0) & (psi <= maxPsi)]
	psi = psi[(psi >= 0) & (psi <= maxPsi)]
	#h = polarplot(np.deg2rad(phi),tand(0.5 * psi),lineSpec)
	if fig is None:
		fig = plt.figure()
	ax = fig.add_subplot(111, projection='polar')
	if int is None:
		#c = ax.scatter(np.deg2rad(phi), bc.tand(0.5 * psi))
		if len(lineSpec) > 0:
			c = plt.polar(np.deg2rad(phi), bc.tand(0.5 * psi), lineSpec)
		else:
			c = plt.polar(np.deg2rad(phi), bc.tand(0.5 * psi))
	else:
		c = ax.scatter(np.deg2rad(phi), bc.tand(0.5 * psi), c = int, cmap='jet')
		fig.colorbar(c, orientation=cbarOri)
	if thetaStep is not None:
		plt.thetagrids(np.array(range(0, 360, thetaStep)))
	if rMax is not None and rStep is not None:
		plt.rgrids(np.array(range(0, rMax, rStep)))
	if not partPlot:
		c[0].axes.set_ylim(0, 1)
	if titleText is not None:
		ax.set_title(titleText)
	plt.tight_layout()  # layout without overlapping
	plt.show()
	return fig, c, ax


def changeStereoProjPlot(plotObj, phi, psi, int=None, maxPsi=90):
	# only allow defined psi values
	phi = phi[(psi >= -maxPsi) & (psi <= maxPsi)]
	if int is not None:
		int = int[(psi >= 0) & (psi <= maxPsi)]
	psi = psi[(psi >= 0) & (psi <= maxPsi)]
	# set new values of plot
	if int is not None:
		pg.setPlotDataPolar(plotObj, bc.tand(0.5 * psi), np.deg2rad(phi), int)
	else:
		pg.setPlotData(plotObj, bc.tand(0.5 * psi), np.deg2rad(phi), int)


# This function gets a spectrum plot (over energy, lattice distances,
# diffraction angles) and adds lines of specified peak positions (can be
# fluorescence or diffraction peaks)
def plotPeakPositions1(x, y, peakPositions=None, peakNames=None, lineCol='r', adaptPosY=True, peakRange=[0.99, 1.01]):
	h = plt.figure()
	plt.plot(x, y, 'k.-')
	maxVal = max(y)
	if peakPositions is not None:
		# select only positions inside data range
		peakNames = peakNames[peakPositions <= max(x)]
		peakPositions = peakPositions[peakPositions <= max(x)]
		plt.vlines([peakPositions, peakPositions], 0, maxVal, lineCol)
		if peakNames is not None:
			for i in range(len(peakPositions)):
				if peakNames[i]:
					peakNames[i] = str(peakNames[i])
				if adaptPosY:
					relVals = (x >= peakRange[0] * peakPositions[i]) & (x <= peakRange[1] * peakPositions[i])
					maxVal = 1.01 * max(y[relVals])
					# plt.vlines([peakPositions[i], peakPositions[i]], 0, maxVal, lineCol)
				plt.text(peakPositions[i], maxVal, peakNames[i])
	plt.grid(True)
	plt.show()


def plotSin2Psi(data, showErr=True):
	symbols = ['s', '^', 'p', 'd', 'v', 'o', '+', 'x', '*', 'h', '<', '>', '.']
	colors = ['r', 'g', 'b', 'c', 'm', 'y', 'r', 'g', 'b', 'c', 'm']
	# extract relevant data
	dVals = data['dVals']
	dErrVals = data['dErr']
	phiVals = data['phiVals']
	psiVals = data['psiVals']
	hklVal = data['hklVal']
	mVals = data['mVals']
	bVals = data['bVals']
	errVals = data['errVals']
	valsAll = data['meanVals']
	sinpsi2Star = data['sinpsi2Star']
	# define derived values
	phiUni = np.sort(np.unique(phiVals))
	sinpsi2 = bc.sind(psiVals) ** 2
	psiUni = np.unique(psiVals)
	psiUniWithoutZero = psiUni[psiUni != 0]
	psiSign = np.sign(psiUniWithoutZero[-1])  # sign of last psi value
	psiUni = psiUni[(np.sign(psiUni) == psiSign) | (psiUni == 0)]  # only negative or positive values
	maxUsedPsi = np.max(np.abs(psiUni))
	sinpsi2Uni = bc.sind(psiUni) ** 2
	sinpsi2Distr = np.arange(0, 1.001, 0.01)
	t = plt.figure()
	for i in range(len(phiUni)):  # for each phi value plot data points
		used = phiVals == phiUni[i]
		if showErr:
			# plt.errorbar(sinpsi2[used], dVals[used], dErrVals[used], fmt=colors[i] + symbols[i], ecolor='k')
			plt.errorbar(sinpsi2[used], dVals[used], dErrVals[used], fmt=colors[i] + symbols[i],
				label='Phi=' + str(phiUni[i]) + '°')
		else:
			plt.plot(sinpsi2[used], dVals[used], colors[i] + symbols[i], label='Phi=' + str(phiUni[i]) + '°')
	# plot regression results
	if mVals[5] != 0 or bVals[5] != 0:  # mean
		xMean = data['xMean']
		if showErr:
			# plt.errorbar(xMean, valsAll[:, 0], valsAll[:, 0] * valsAll[:, 1],
			# fmt=colors[4] + symbols[4], ecolor='k')
			plt.errorbar(xMean, valsAll[:, 0], valsAll[:, 0] * valsAll[:, 1], fmt=colors[4] + symbols[4],
				label='mean values')
		else:
			plt.plot(xMean, valsAll[:, 0], colors[4] + symbols[4], label='mean values')
		plt.plot(xMean, data['yMean'], 'k.')
	if mVals[4] != 0 or bVals[4] != 0:  # s12
		plt.plot(data['x12'], data['y12'], 'k.')
	if mVals[0] != 0 or bVals[0] != 0:  # s11
		plt.plot(data['x1'], data['y1'], 'k.')
	if mVals[1] != 0 or bVals[1] != 0:  # s22
		plt.plot(data['x2'], data['y2'], 'k.')
	if mVals[5] != 0 or bVals[5] != 0:  # mean
		plt.plot([0, 1], [bVals[5], mVals[5] + bVals[5]], 'k-')
	if mVals[4] != 0 or bVals[4] != 0:  # s12
		plt.plot([0, 1], [bVals[4], mVals[4] + bVals[4]], 'k-')
	if mVals[0] != 0 or bVals[0] != 0:  # s11
		plt.plot([0, 1], [bVals[0], mVals[0] + bVals[0]], 'k-')
	if mVals[1] != 0 or bVals[1] != 0:  # s22
		plt.plot([0, 1], [bVals[1], mVals[1] + bVals[1]], 'k-')
	if mVals[5] != 0 or bVals[5] != 0:  # mean
		plt.plot(sinpsi2Distr, mVals[0] * sinpsi2Distr + mVals[2] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[5], 'k:')
		plt.plot(sinpsi2Distr, mVals[0] * sinpsi2Distr - mVals[2] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[5], 'k:')
		plt.plot(sinpsi2Distr, mVals[1] * sinpsi2Distr + mVals[3] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[5], 'k:')
		plt.plot(sinpsi2Distr, mVals[1] * sinpsi2Distr - mVals[3] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[5], 'k:')
	else:
		if mVals[0] != 0 or bVals[0] != 0:  # s11
			plt.plot(sinpsi2Distr,
				mVals[0] * sinpsi2Distr + mVals[2] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[0], 'k:')
			plt.plot(sinpsi2Distr,
				mVals[0] * sinpsi2Distr - mVals[2] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[0], 'k:')
		if mVals[1] != 0 or bVals[1] != 0:  # s22
			plt.plot(sinpsi2Distr,
				mVals[1] * sinpsi2Distr + mVals[3] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[1], 'k:')
			plt.plot(sinpsi2Distr,
				mVals[1] * sinpsi2Distr - mVals[3] * 2 * (sinpsi2Distr * (1 - sinpsi2Distr)) ** 0.5 + bVals[1], 'k:')
	# plot line at sin2psiStar
	# plt.vlines(bf.ones(2, 1) * sinpsi2Star, np.min(dVals), np.max(dVals), 'k', 'dashed')
	plt.grid()
	plt.xlabel('sin^2 psi')
	plt.ylabel('d in nm')
	plt.legend()
	plt.title('sin^2 psi curve for (' + str(hklVal) + ') peak')
	plt.tight_layout()  # layout without overlapping
	plt.show()  # saveas(gcf,[pathName,'Auswertung\Sin2Psi_',num2str(p),'.fig'])


def plotMultiWavelength(data, showErr=True):
	hklNames = bf.getKeyList(data)
	for hkl in hklNames:
		plotSin2Psi(data[hkl], showErr)


def plotUniversalPlot(data, showErr=True):
	symbols = ['s', '^', 'p', 'd', 'v', 'o', '+', 'x', '*', 'h', '<', '>', '.']
	colors = ['r', 'g', 'b', 'c', 'm', 'y', 'r', 'g', 'b', 'c', 'm']
	# extract relevant data
	tauVals = data['tauVals']
	psiVals = data['psiVals']
	stresses = bf.getDictValOrDef(data, 'stresses')
	accuracy = bf.getDictValOrDef(data, 'accuracy')
	stressNames = ['s11-s33', 's22-s33', 's13', 's23']
	if stresses is not None and accuracy is not None:
		for i in range(bf.size(stresses, 1)):
			curStresses = np.round(stresses[:, i])
			if sum(curStresses) != 0 or max(curStresses) != 0 or min(curStresses) != 0:
				if showErr:
					# pg.plotErrData(curStresses, np.round(accuracy[:, i]), tauVals, 'rp', 'on', 'Information depths in um',
					# 	'Residual stresses in MPa', 'Residual stresses ' + stressNames[i], 'k')
					pg.plotErrData(curStresses, np.round(accuracy[:, i]), tauVals, 'rp', 'on', 'Information depths in um',
						'Residual stresses in MPa', 'Residual stresses ' + stressNames[i])
				else:
					pg.plotData(curStresses, tauVals, 'rp', 'on', 'Information depths in um',
						'Residual stresses in MPa', 'Residual stresses ' + stressNames[i])
	else:
		for stressName in stressNames:
			stressVals = bf.getDictValOrDef(data, stressName)
			accuracyVals = bf.getDictValOrDef(data, 'dev_' + stressName)
			if stressVals is not None and accuracyVals is not None:
				if sum(stressVals) != 0 or max(stressVals) != 0 or min(stressVals) != 0:
					if showErr:
						# pg.plotErrData(np.round(stressVals), np.round(accuracyVals), tauVals, 'rp', 'on', 'Information depths in um',
						# 	'Residual stresses in MPa', 'Residual stresses ' + stressName, 'k')
						pg.plotErrData(np.round(stressVals), np.round(accuracyVals), tauVals, 'rp', 'on',
							'Information depths in um', 'Residual stresses in MPa',
							'Residual stresses ' + stressName)
					else:
						pg.plotData(np.round(stressVals), tauVals, 'rp', 'on', 'Information depths in um',
							'Residual stresses in MPa', 'Residual stresses ' + stressName)


def plotStresses(data, showErr=True):
	hklList = data['hklList']
	tauMean = data['tauMean']
	stresses = bf.getDictValOrDef(data, 'stresses')
	accuracy = bf.getDictValOrDef(data, 'accuracy')
	stressNames = ['s11-s33', 's22-s33', 's13', 's23', 's12', 's33']
	if stresses is not None and accuracy is not None:
		for i in range(bf.size(stresses, 1)):
			curStresses = np.round(stresses[:, i])
			if sum(curStresses) != 0 or max(curStresses) != 0 or min(curStresses) != 0:
				if showErr:
					# pg.plotErrData(curStresses, np.round(accuracy[:, i]), tauMean, 'ro-', 'on', 'Information depths in um',
					# 	'Residual stresses in MPa', 'Residual stresses ' + stressNames[i], 'k')
					pg.plotErrData(curStresses, np.round(accuracy[:, i]), tauMean, 'ro-', 'on',
						'Information depths in um', 'Residual stresses in MPa',
						'Residual stresses ' + stressNames[i])
				else:
					pg.plotData(curStresses, tauMean, 'ro-', 'on', 'Information depths in um',
						'Residual stresses in MPa', 'Residual stresses ' + stressNames[i])
	else:
		for stressName in stressNames:
			stressVals = bf.getDictValOrDef(data, stressName)
			accuracyVals = bf.getDictValOrDef(data, 'dev_' + stressName)
			if stressVals is not None and accuracyVals is not None:
				if sum(stressVals) != 0 or max(stressVals) != 0 or min(stressVals) != 0:
					if showErr:
						# pg.plotErrData(np.round(stressVals), np.round(accuracyVals), tauMean, 'ro-', 'on', 'Information depths in um',
						# 	'Residual stresses in MPa', 'Residual stresses ' + stressName, 'k')
						pg.plotErrData(np.round(stressVals), np.round(accuracyVals), tauMean,
							'ro-', 'on', 'Information depths in um', 'Residual stresses in MPa',
							'Residual stresses ' + stressName)
					else:
						pg.plotData(np.round(stressVals), tauMean, 'rp-', 'on', 'Information depths in um',
							'Residual stresses in MPa', 'Residual stresses ' + stressName)


def plotStrainFreeLatticeSpacing(data, showErr=True):
	hklList = data['hklList']  # perhaps used to plot as text at each data point
	tauMean = data['tauMean']
	aStarVals = data['dStar100']
	if showErr:
		aStarErrVals = data['dStar100Err']
		# pg.plotErrData(aStarVals, aStarErrVals, tauMean, 'ro-', 'on', 'Information depths in um',
		# 	'a* in nm', 'k')
		pg.plotErrData(aStarVals, aStarErrVals, tauMean, 'ro-', 'on', 'Information depths in um',
			'a* in nm')
	else:
		pg.plotData(aStarVals, tauMean, 'ro-', 'on', 'Information depths in um', 'a* in nm')


def plotInverseLaplaceResults(tauVals, stressVals, tauStresses, zStresses, zVals=None, stressValsErr=None, title=None):
	h = pg.figure(7.2, 4.6)
	if stressValsErr is None:
		plt.plot(tauVals, stressVals, 'r*')
	else:
		pg.plotErrData(stressVals, stressValsErr, tauVals, 'r*', 'on')
		#pg.plotErrData(stressVals, stressValsErr, tauVals, 'r*', 'on', ecol='k')
	plt.plot(tauVals, tauStresses, 'r-')
	if zVals is None:
		plt.plot(tauVals, zStresses, 'bo-')
	else:
		plt.plot(zVals, zStresses, 'bo-')
	plt.grid()
	plt.xlabel('Information depth and sample depth in um')
	plt.ylabel('Residual stresses in MPa')
	plt.legend(('measured', 'fitted', 'calculated'))
	if title is not None:
		plt.title(title)
	plt.tight_layout()  # layout without overlapping
	plt.show()
