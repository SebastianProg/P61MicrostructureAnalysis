import math
import numpy as np
import scipy.optimize as op

import basics.functions as bf
import basics.calculations as bc
import basics.datamodels as dm
import filehandling.general as fg
import diffraction.conversions as conv


# calculate tau specifying psi and eta (equation of scattering vector-method)
def calcTauPsiEta(mu, tth, psi, eta):
	th = tth / 2
	return (bc.sind(th) ** 2 - bc.sind(psi) ** 2 + bc.cosd(th) ** 2 * bc.sind(psi) ** 2 * bc.sind(eta) ** 2) / \
		   (2 * mu * bc.sind(th) * bc.cosd(psi))


def calc3Gamma(h, k, l):
	return 3 * (h ** 2 * k ** 2 + k ** 2 * l ** 2 + l ** 2 * h ** 2) / (h ** 2 + k ** 2 + l ** 2) ** 2


def calc3Gamma2(hkl):
	vals = conv.splitHkl(hkl)
	return calc3Gamma(vals[0], vals[1], vals[2])


def calcEtaTau(mu, tau, tth, psi):
	th = tth / 2
	return bc.asind(2 * mu * tau * bc.sind(th) * bc.cosd(psi) - bc.sind(th) ** 2 + bc.sind(psi) ** 2) ** 0.5 / \
		   (bc.cosd(th) * bc.sind(psi))


def calcEtaMin(tth, psi):
	th = tth / 2
	#	if psi <= th:
	#		etaMin = 0
	#	else:
	#		etaMin = asind((sind(psi)**2 - sind(th)**2)**0.5 / (bc.cosd(th) * sind(psi)))
	if psi == 0:
		return 0
	else:
		return bc.asind((max(0, bc.sind(psi) ** 2 - bc.sind(th) ** 2)) ** 0.5 / (bc.cosd(th) * bc.sind(psi)))


# calculate min and max tau for given psi value
def calcValidTauRange(mu, tth, psi):
	th = tth / 2
	etaMin = calcEtaMin(tth, psi)
	tauMin = calcTauPsiEta(mu, tth, psi, etaMin)
	# Psi mode
	tauMax = (bc.sind(th) ** 2 - bc.sind(psi) ** 2 + bc.cosd(th) ** 2 * bc.sind(psi) ** 2) / \
			 (2 * mu * bc.sind(th) * bc.cosd(psi))
	return etaMin, tauMax


# calculate min and max tau for given psi value
def calcTauRange(mu, tth, psi, etaMin):
	th = tth / 2
	tauMin = calcTauPsiEta(mu, tth, psi, etaMin)
	if abs(psi) > th:
		tauMin = None  # NaN heißt nichts?
	else:
		# Omega mode
		tauMin = (bc.sind(th) ** 2 - bc.sind(psi) ** 2 + bc.cosd(th) ** 2 * bc.sind(psi) ** 2) / \
				 (2 * mu * bc.sind(th) * bc.cosd(psi))
	# Psi mode
	tauMax = (bc.sind(th) ** 2 - bc.sind(psi) ** 2 + bc.cosd(th) ** 2 * bc.sind(psi) ** 2) / \
			 (2 * mu * bc.sind(th) * bc.cosd(psi))
	return tauMin, tauMax


# calculate tau specifying alpha and beta (general equation)
def calcTauAlphaBeta(mu, alpha, beta):
	return bc.sind(alpha) * bc.sind(beta) / (mu * (bc.sind(alpha) + bc.sind(beta)))


def calcSigma33(d, d0, s1, hs2):
	return (d / d0 - 1) / (hs2 + 3 * s1)


def calcMue(energy, material):
	# energy in keV
	# mue in 1/um
	mue = np.zeros(bf.size(energy))
	if bf.size(material, 0) > 1:
		# matrix of absorption data
		mue = bf.ones(bf.size(energy)) * np.nan
		for i in range(bf.size(material, 0)):
			mue[energy >= material[i, 0] & energy < material[i, 1]] = material[i, 2] * energy[energy >= material[i, 0] &
				energy < material[i, 1]] ** material[i, 3] / 10000
	else:
		# name of material
		if material == 'Fe':
			#mue = 0.0006 .* (energy ./ 1000).^-2.736 .* 7.89 / 10000; % in um
			mue[energy < 7.11138755] = 77083 * energy[energy < 7.11138755] ** -2.66 / 10000 # in um
			mue[energy >= 7.11138755] = 739152 * energy[energy >= 7.11138755] ** -2.755 / 10000 # in um
		elif material == 'Ni':
			mue[energy < 8.3328] = 95880 * energy[energy < 8.3328] ** -2.556 / 10000 # in um
			mue[energy >= 8.3328] = 929809 * energy[energy >= 8.3328] ** -2.708 / 10000 # in um
		else:
			mue = 10 ** 6
	if bf.size(mue,0) < bf.size(mue,1):
		mue = np.transpose(mue)
	return mue


def calcF33(s1, hs2):
	return 3 * s1 + hs2


def calcF33_2(dek):
	return calcF33(dek[1], dek[2])


def s33FreeOrientation(s1, hs2, dPsi0, d0, dStar, dPlus):
	return bc.asind((dPsi0 - d0) * calcF33(s1, hs2) / hs2 * (dStar - dPlus)) ** 0.5


def s11s22FreeOrientation(s1, hs2):
	return bc.asind(-2 * s1 / hs2) ** 0.5


def latticeSpacings(spacings, angles, h, k, l):
	# check if single values are arrays or not
	if bf.length(spacings) == 1:
		if len(bf.size(spacings)) > 0:
			spacing = spacings[0]
		else:
			spacing = spacings
	if bf.length(angles) == 1:
		if len(bf.size(angles)) > 0:
			angle = angles[0]
		else:
			angle = angles
	# determine lattice spacings
	if bf.length(spacings) == 1 and bf.length(angles) == 1 and angle == 90:
		# cubic
		dVals = conv.aVals2latticeDists(spacing, h, k, l)
	elif bf.length(spacings) == 1 and bf.length(angles) == 1 and angle != 90:
		# rhomboedric
		dVals = (spacing ** 2 * (1 - 3 * bc.cosd(angle) ** 2 + 2 * bc.cosd(angle) ** 3) / ((h ** 2 + k ** 2 + l ** 2) *
																						   bc.sind(angle) ** 2 + 2 * (h * k + k * l + h * l) * (bc.cosd(angle) ** 2 - bc.cosd(angle)))) ** 0.5
	elif bf.length(spacings) == 2 and bf.length(angles) == 1 and angle == 90:
		# tetragonal
		#dVals = spacings[0] / (h ** 2 + k ** 2 + l ** 2 * (spacings[0] ** 2 / spacings[1] ** 2)) ** 0.5
		dVals = 1 / ((h ** 2 + k ** 2) / spacings[0] ** 2 + l ** 2 / spacings[1] ** 2) ** 0.5
	elif bf.length(spacings) == 2 and bf.length(angles) == 2 and angles[0] == 90 and angles[1] == 120:
		# hexagonal
		dVals = spacings[0] / (4 / 3 * (h ** 2 + h * k + k ** 2) + l ** 2 * spacings[0] ** 2 / spacings[1] ** 2) ** 0.5
	elif bf.length(spacings) == 3 and bf.length(angles) == 1 and angle == 90:
		# orthorhombic
		dVals = 1 / ((h / spacings[0]) ** 2 + (k / spacings[1]) ** 2 + (l / spacings[2]) ** 2) ** 0.5
	elif bf.length(spacings) == 3 and bf.length(angles) == 2 and angles[0] == 90:  # and angles[1] != 90
		# monoklin
		dVals = 1 / (h ** 2 / (spacings[0] ** 2 * bc.sind(angles[1]) ** 2) + k ** 2 / spacings[1] ** 2 + l ** 2 /
					 (spacings[2] ** 2 * bc.sind(angles[1]) ** 2) - 2 * h * l * bc.cosd(angles[1]) / (spacings[0] * spacings[2] *
																									  bc.sind(angles[1]) ** 2)) ** 0.5
	elif bf.length(spacings) == 3 and bf.length(angles) == 3:
		# triklin
		v = spacings[0] * spacings[1] * spacings[2] * (1 - bc.cosd(angles[0]) ** 2 - bc.cosd(angles[1]) ** 2 -
													   bc.cosd(angles[2]) ** 2 + 2 * bc.cosd(angles[0]) * bc.cosd(angles[1]) * bc.cosd(angles[2])) ** 0.5
		s11 = spacings[1] ** 2 * spacings[2] ** 2 * bc.sind(angles[0]) ** 2
		s22 = spacings[0] ** 2 * spacings[2] ** 2 * bc.sind(angles[1]) ** 2
		s33 = spacings[0] ** 2 * spacings[1] ** 2 * bc.sind(angles[2]) ** 2
		s12 = spacings[0] * spacings[1] * spacings[2] ** 2 * (bc.cosd(angles[0]) * bc.cosd(angles[1]) - bc.cosd(angles[2]))
		s23 = spacings[0] ** 2 * spacings[1] * spacings[2] * (bc.cosd(angles[1]) * bc.cosd(angles[2]) - bc.cosd(angles[0]))
		s13 = spacings[0] * spacings[1] ** 2 * spacings[2] * (bc.cosd(angles[2]) * bc.cosd(angles[0]) - bc.cosd(angles[1]))
		dVals = (v ** 2 / (s11 * h ** 2 + s22 * k ** 2 + s33 * l ** 2 + 2 * s12 * h * k + 2 * s23 * k * l + 2 *
			s13 * h * l)) ** 0.5
	return dVals


# spacings = [a b c];
# angles = [aAngle, bAngle, cAngle];
# spacings = [0.28665 0.28665 0.28665];
# angles = [90 90 90];
# spacings = 0.28665;
# angles = 90;
# hklVals = 110;
# dVals = latticeSpacings2(spacings, angles, hklVals);
def latticeSpacings2(spacings, angles, hklVals):
	# split hkl values
	d = conv.splitHkl(hklVals)
	return latticeSpacings(spacings, angles, d[0], d[1], d[2])


def latticeSpacings3(spacings, angles, hklVector):
	# split hkl values
	return latticeSpacings(spacings, angles, hklVector[:, 0], hklVector[:, 1], hklVector[:, 2])


def calcDeadTime(realTime, liveTime):
	return (1 - liveTime / realTime) * 100


def sin2PsiAnalysis(data, maxPsi=None):
	# data: dVals, errVals, tauVals, phiVals, psiVals, hklVal, s1Val, hs2Val, ibVals???
	# extract needed data
	dVals = data['dVals']
	dErrVals = data['dErr']
	tauVals = data['tauVals']
	phiVals = data['phiVals']
	psiVals = data['psiVals']
	ibVals = data['ibVals']
	hklVal = data['hklVal']
	s1Val = data['s1Val']
	hs2Val = data['hs2Val']
	a0Val = bf.getDictValOrDef(data, 'a0Val')
	# define derived values
	phiUni = np.sort(np.unique(phiVals))
	sinpsi2 = bc.sind(psiVals) ** 2
	sinpsi2Distr = np.arange(0, 1.001, 0.01)
	psiUni = np.unique(psiVals)
	psiUniWithoutZero = psiUni[psiUni != 0]
	psiSign = np.sign(psiUniWithoutZero[-1])  # sign of last psi value
	psiUni = psiUni[(np.sign(psiUni) == psiSign) | (psiUni == 0)]  # only negative or positive values
	if maxPsi is not None:
		psiUni = psiUni[(psiUni <= maxPsi) & (psiUni >= -maxPsi)]  # take only values smaller than or equal to |maxPsi°|
	maxUsedPsi = np.max(np.abs(psiUni))
	sinpsi2Uni = bc.sind(psiUni) ** 2
	sin2psiUni = bc.sind(np.abs(2 * psiUni))
	# define global regression values
	mVals = np.zeros(6)
	bVals = np.zeros(6)
	errVals = np.zeros(6)
	errValsStar = 0
	# perform linear regressions for one peak
	sinpsi2StarSingle = s1Val / hs2Val
	sinpsi2Star = -2 * s1Val / hs2Val
	sinpsi2Plus = 1 + s1Val / hs2Val
	dValsValid = (dVals > 0) & (np.isnan(dVals) == False)
	tauValsValid = (tauVals >= 0) & (tauVals < 1e6)
	ibValsValid = (ibVals > 0) & (ibVals < 1000)
	if maxPsi is None:
		tauMean = np.sum(tauVals[tauValsValid]) / len(tauVals[tauValsValid])
		# tauMean = (max(tauVals[tauValsValid]) + min(tauVals[tauValsValid])) / 2
	else:
		# take only values smaller than or equal to | maxPsi° |
		curTau = tauVals[(psiVals <= maxPsi) & (psiVals >= -maxPsi) & tauValsValid]
		tauMean = np.sum(curTau) / len(curTau)
		# tauMean = (max(curTau) + min(curTau)) / 2
	# determine mean values of IB
	valuesIb = np.zeros(len(phiUni))
	for j in range(len(phiUni)):
		if maxPsi is None:
			valuesIb[j] = np.mean(ibVals[phiVals == phiUni[j] & ibValsValid])
		else:
			valuesIb[j] = np.mean(ibVals[(phiVals == phiUni[j]) & (psiVals <= maxPsi) & (psiVals >= -maxPsi) & ibValsValid])
	# perform linear regression of all phi values
	valsAll = np.zeros((len(psiUni), 2))
	vals45 = np.zeros((len(psiUni), 2))
	vals0_180 = np.zeros((len(psiUni), 3))
	vals90_270 = np.zeros((len(psiUni), 3))
	for i in range(len(psiUni)):
		if bf.max(bf.containsItems(psiVals, psiUni[i])[0])[0] and bf.max(bf.containsItems(psiVals, -psiUni[i])[0])[
			0]:
			# positive and negative psi values
			val0 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 0) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 0) & dValsValid]])
			val90 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 90) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 90) & dValsValid]])
			val180 = np.array([dVals[(psiVals == -psiUni[i]) & (phiVals == 0) & dValsValid],
				dErrVals[(psiVals == -psiUni[i]) & (phiVals == 0) & dValsValid]])
			val270 = np.array([dVals[(psiVals == -psiUni[i]) & (phiVals == 90) & dValsValid],
				dErrVals[(psiVals == -psiUni[i]) & (phiVals == 90) & dValsValid]])
		else:
			val0 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 0) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 0) & dValsValid]])
			val90 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 90) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 90) & dValsValid]])
			val180 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 180) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 180) & dValsValid]])
			val270 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 270) & dValsValid],
				dErrVals[(psiVals == psiUni[i]) & (phiVals == 270) & dValsValid]])
		val45 = np.array([dVals[(psiVals == psiUni[i]) & (phiVals == 45) & dValsValid],
			dErrVals[(psiVals == psiUni[i]) & (phiVals == 45) & dValsValid]])
		if val45.size > 0:
			vals45[i, :] = val45
		if val0.size > 0 and val180.size > 0:
			vals0_180[i, 0] = 0.5 * (val0[0] + val180[0])
			vals0_180[i, 1] = 0.5 * (val0[0] - val180[0])
			vals0_180[i, 2] = 0.5 * (val0[1] + val180[1])
		elif val0.size > 0:
			vals0_180[i, 0] = val0[0]
			vals0_180[i, 1] = 0
			vals0_180[i, 2] = val0[1]
		elif val180.size > 0:
			vals0_180[i, 0] = val180[0]
			vals0_180[i, 1] = 0
			vals0_180[i, 2] = val180[1]
		if val90.size > 0 and val270.size > 0:
			vals90_270[i, 0] = 0.5 * (val90[0] + val270[0])
			vals90_270[i, 1] = 0.5 * (val90[0] - val270[0])
			vals90_270[i, 2] = 0.5 * (val90[1] + val270[1])
		elif val90.size > 0:
			vals90_270[i, 0] = val90[0]
			vals90_270[i, 1] = 0
			vals90_270[i, 2] = val90[1]
		elif val270.size > 0:
			vals90_270[i, 0] = val270[0]
			vals90_270[i, 1] = 0
			vals90_270[i, 2] = val270[1]
		if vals0_180[i, 0] != 0 and vals90_270[i, 0] != 0:
			valsAll[i, :] = 0.5 * (vals0_180[i, [0, 2]] + vals90_270[i, [0, 2]])
	# check for validity
	# perform linear regression for all phi values
	usedIndex = np.array(range(len(psiUni)))
	used = usedIndex[valsAll[:, 0] != 0]
	if len(used) > 0:
		if np.all(valsAll[used, 1] != 0):
			[m, b, CoD, yD, err, err33, bErr] = dm.linRegWeighted(sinpsi2Uni[used], valsAll[used, 0], 1 / valsAll[used, 1])
		else:
			[m, b, CoD, yD, err, err33, bErr] = dm.linReg(sinpsi2Uni[used], valsAll[used, 0])
		mVals[5] = m
		bVals[5] = b
		errVals[5] = err33
		errValsStar = bErr
	# perform linear regression of phi values with phi = 45° (if exist)
	used12 = usedIndex[vals45[:, 0] != 0]
	if len(used12) > 0:
		if len(used) > 0:
			if np.all(vals45[used12, 1] != 0):
				[m12, b12, CoD12, yD12, err, err12, bErr] = dm.linRegWeighted(sinpsi2Uni[used12], vals45[used12, 0],
					1 / vals45[used12, 1], 'b', b)
			else:
				[m12, b12, CoD12, yD12, err, err12, bErr] = dm.linReg(sinpsi2Uni[used12], vals45[used12, 0], 'b', b)
		else:
			if np.all(vals45[used12, 1] != 0):
				[m12, b12, CoD12, yD12, err, err12, bErr] = dm.linRegWeighted(sinpsi2Uni[used12], vals45[used12, 0],
					1 / vals45[used12, 1])
			else:
				[m12, b12, CoD12, yD12, err, err12, bErr] = dm.linReg(sinpsi2Uni[used12], vals45[used12, 0])
			mVals[5] = m12
			bVals[5] = b12
		mVals[4] = m12
		bVals[4] = b12
		errVals[4] = err12
	# perform linear regression of phi values with phi = 0°/180°
	used1 = usedIndex[vals0_180[:, 0] != 0]
	if len(used1) > 0:
		if len(used) > 0:
			if np.all(vals0_180[used1, 2] != 0):
				[m1, b1, CoD1, yD1, err, err1, bErr] = dm.linRegWeighted(sinpsi2Uni[used1], vals0_180[used1, 0],
					1 / vals0_180[used1, 2], 'b', b)
			else:
				[m1, b1, CoD1, yD1, err, err1, bErr] = dm.linReg(sinpsi2Uni[used1], vals0_180[used1, 0], 'b', b)
		else:
			if np.all(vals0_180[used1, 2] != 0):
				[m1, b1, CoD1, yD1, err, err1, bErr] = dm.linRegWeighted(sinpsi2Uni[used1], vals0_180[used1, 0],
					1 / vals0_180[used1, 2])
			else:
				[m1, b1, CoD1, yD1, err, err1, bErr] = dm.linReg(sinpsi2Uni[used1], vals0_180[used1, 0])
			mVals[5] = m1
			bVals[5] = b1
		# [m1f, err1f] = lsqcurvefit(linM,m1,sinpsi2Uni(used1),vals0_180(used1,1));
		mVals[0] = m1
		bVals[0] = b1
		errVals[0] = err1
		if np.all(vals0_180[used1, 2] != 0):
			[m13, b13, CoD13, yD13, err, err13, bErr] = dm.linRegWeighted(sin2psiUni[used1], vals0_180[used1, 1],
				1 / vals0_180[used1, 2], 'b', 0)
		else:
			[m13, b13, CoD13, yD13, err, err13, bErr] = dm.linReg(sin2psiUni[used1], vals0_180[used1, 1], 'b', 0)
		# [m1f, err1f] = lsqcurvefit(linBase,m1,sin2psiUni(used1),vals0_180(used1,2));
		mVals[2] = m13
		bVals[2] = b13
		errVals[2] = err13
	# perform linear regression of phi values with phi = 90°/270°
	used2 = usedIndex[vals90_270[:, 0] != 0]
	if len(used2) > 0:
		if len(used) > 0:
			if np.all(vals90_270[used2, 2] != 0):
				[m2, b2, CoD2, yD2, err, err2, bErr] = dm.linRegWeighted(sinpsi2Uni[used2], vals90_270[used2, 0],
					1 / vals90_270[used2, 2], 'b', b)
			else:
				[m2, b2, CoD2, yD2, err, err2, bErr] = dm.linReg(sinpsi2Uni[used2], vals90_270[used2, 0], 'b', b)
		else:
			if np.all(vals90_270[used2, 2] != 0):
				[m2, b2, CoD2, yD2, err, err2, bErr] = dm.linRegWeighted(sinpsi2Uni[used2], vals90_270[used2, 0],
					1 / vals90_270[used2, 2])
			else:
				[m2, b2, CoD2, yD2, err, err2, bErr] = dm.linReg(sinpsi2Uni[used2], vals90_270[used2, 0])
			mVals[5] = m2
			bVals[5] = b2
		# [m2f, err2f] = lsqcurvefit(linM,m2,sinpsi2Uni(used2),vals90_270(used2,1));
		mVals[1] = m2
		bVals[1] = b2
		errVals[1] = err2
		if np.all(vals90_270[used2, 2] != 0):
			[m23, b23, CoD23, yD23, err, err23, bErr] = dm.linRegWeighted(sin2psiUni[used2], vals90_270[used2, 1],
				1 / vals90_270[used2, 2], 'b', 0)
		else:
			[m23, b23, CoD23, yD23, err, err23, bErr] = dm.linReg(sin2psiUni[used2], vals90_270[used2, 1], 'b', 0)
		# [m2f, err2f] = lsqcurvefit(linBase,m2,sin2psiUni(used2),vals90_270(used2,2));
		mVals[3] = m23
		bVals[3] = b23
		errVals[3] = err23
	# determine stresses
	dStarVal0 = mVals[0] * sinpsi2StarSingle + bVals[0]
	dStarVal90 = mVals[1] * sinpsi2StarSingle + bVals[1]
	dStarVal = mVals[5] * sinpsi2Star + bVals[5]
	dPlusVal = mVals[5] * sinpsi2Plus + bVals[5]
	dComVal = mVals[5] * 2 / 3 + bVals[5]
	s11_s33 = mVals[0] / (hs2Val * dStarVal)
	ds11 = 2 * errVals[0] ** 0.5 / (hs2Val * dStarVal)
	s22_s33 = mVals[1] / (hs2Val * dStarVal)
	ds22 = 2 * errVals[1] ** 0.5 / (hs2Val * dStarVal)
	s13 = mVals[2] / (hs2Val * dStarVal)
	ds13 = 2 * errVals[2] ** 0.5 / (hs2Val * dStarVal)
	s23 = mVals[3] / (hs2Val * dStarVal)
	ds23 = 2 * errVals[3] ** 0.5 / (hs2Val * dStarVal)
	if np.any((phiVals == 45) | (phiVals == -45) | (phiVals == 225)):
		s12 = mVals[4] / (hs2Val * dStarVal) - 0.5 * (s11_s33 + s22_s33)
	else:
		s12 = mVals[4] / (hs2Val * dStarVal)
	ds12 = 2 * errVals[4] ** 0.5 / (hs2Val * dStarVal)
	aStarVal = conv.latticeDists2aVals2(dStarVal, hklVal)
	if a0Val is None or a0Val == 0:
		s33 = 0
		ds33 = 0
	else:
		s33 = calcSigma33(aStarVal, a0Val, s1Val, hs2Val)
		ds33 = 2 * errVals[5] ** 0.5 / (hs2Val * dStarVal)
	# d0_1 = regVals(:,2) ./ (dekList(:,2) .* (s11_s33 + s22_s33) + s33 .* f33 + 1);
	#     [sMain, d0, quality] = stressesWithS33(dekList, dStarVals, [s11_s33 s22_s33 s12 s13 s23], ...
	#         [data(:,1) dVals error data(:,[7 8])], plotData);
	# combine results
	stresses = np.array([s11_s33, s22_s33, s13, s23, s12, s33])
	accuracy = np.array([ds11, ds22, ds13, ds23, ds12, ds33])
	# resData: tauMean, dStar, stresses, accuracy, mVals???, bVals???
	resData = {'tauMean': tauMean, 'dStar100': aStarVal, 'dStar100Err': 2 * errValsStar ** 0.5,
		'stresses': stresses, 'accuracy': accuracy, 'meanIB': valuesIb}
	plotData = {'dVals': data['dVals'], 'dErr': data['dErr'], 'phiVals': data['phiVals'],
		'psiVals': data['psiVals'], 'hklVal': hklVal, 'mVals': mVals, 'bVals': bVals, 'errVals': errVals,
		'sinpsi2Star': sinpsi2Star, 'meanVals': valsAll[used,:]}
	if len(used) > 0:
		plotData['xMean'] = sinpsi2Uni[used]
		plotData['yMean'] = yD
	if len(used12) > 0:
		plotData['x12'] = sinpsi2Uni[used12]
		plotData['y12'] = yD12
	if len(used1) > 0:
		plotData['x1'] = sinpsi2Uni[used1]
		plotData['y1'] = yD1
	if len(used2) > 0:
		plotData['x2'] = sinpsi2Uni[used2]
		plotData['y2'] = yD2
	return resData, plotData


#
def multiWavelengthAnalysis(data, maxPsi=None):
	keyList = bf.getKeyList(data)
	peakCount = int(bf.replace(keyList[-1].split('_')[0], 'pv')) + 1  # must be adapted in further versions!!!
	phiVals = data['phi']
	psiVals = data['psi']
	tauMean = bf.zeros(peakCount)
	hklList = bf.zeros(peakCount, dtypeVal=np.int)
	s1Dec = bf.zeros(peakCount)
	hs2Dec = bf.zeros(peakCount)
	dStarVals = bf.zeros(peakCount)
	dStarErrVals = bf.zeros(peakCount)
	stresses = np.zeros((peakCount, 6))  # s11-33 s22-33 s13 s23 s12 s33
	accuracy = np.zeros((peakCount, 6))
	integralWidth = np.zeros((peakCount, 6))  # phi0 phi90 phi180 phi270 phi45 phi225
	plotData = dict()
	for p in range(peakCount):  # perform this for all peaks
		ibVals = data['pv' + str(p) + '_sigma']
		centerVals = data['pv' + str(p) + '_center']
		centerErrVals = data['pv' + str(p) + '_center_err']
		dMinVals = conv.energies2latticeDists(centerVals - centerErrVals, data['tth'])
		dMaxVals = conv.energies2latticeDists(centerVals + centerErrVals, data['tth'])
		dErrVals = np.abs(dMaxVals - dMinVals) / 2
		tauVals = data['pv' + str(p) + '_depth']
		dVals = data['pv' + str(p) + '_dspac'] / 10
		#hklVals = data['pv' + str(p) + '_hklList']  # first version with adaptions
		#s1Vals = data['pv' + str(p) + '_s1List']  # first version with adaptions
		#hs2Vals = data['pv' + str(p) + '_s2List']  # first version with adaptions
		#hklList[p] = hklVals[0]
		hVals = data['pv' + str(p) + '_h']  # second version
		kVals = data['pv' + str(p) + '_k']  # second version
		lVals = data['pv' + str(p) + '_l']  # second version
		s1Vals = data['pv' + str(p) + '_s1']  # second version
		#hs2Vals = data['pv' + str(p) + '_s2']  # second version
		hs2Vals = data['pv' + str(p) + '_hs2']  # third version
		hklList[p] = conv.mergeHkl(hVals[0], kVals[0], lVals[0])
		s1Dec[p] = s1Vals[0]
		#hs2Dec[p] = hs2Vals[0] * 0.5  # test valid for first and second version!!!!!
		hs2Dec[p] = hs2Vals[0]
		curData = {'dVals': dVals, 'dErr': dErrVals, 'tauVals': tauVals, 'phiVals': phiVals, 'psiVals': psiVals,
			'hklVal': hklList[p], 's1Val': s1Dec[p], 'hs2Val': hs2Dec[p], 'ibVals': ibVals}
		bf.extendDictionary(curData, data, ('a0Val',))
		# perform sin2psi analysis for current peak data
		curResData, curPlotData = sin2PsiAnalysis(curData, maxPsi)
		# remember results
		tauMean[p] = curResData['tauMean']
		dStarVals[p] = curResData['dStar100']
		dStarErrVals[p] = curResData['dStar100Err']
		stresses[p] = curResData['stresses']
		accuracy[p] = curResData['accuracy']
		curMeanIB = curResData['meanIB']
		integralWidth[p, 0:len(curMeanIB)] = curMeanIB
		plotData[str(hklList[p])] = curPlotData
	# resData: hklList, s1Dec, hs2Dec, tauMean, aStar, aStarErr, stresses, accuracy, mVals???, bVals???
	resData = {'hklList': hklList, 's1Dec': s1Dec, 'hs2Dec': hs2Dec, 'tauMean': tauMean, 'dStar100': dStarVals,
		'dStar100Err': dStarErrVals, 'stresses': stresses, 'accuracy': accuracy, 'integralWidth': integralWidth}
	return resData, plotData


def universalPlotAnalysis(data, maxPsi=None, minDistPsiStar=0.15, minValPsiNormal=0.08,
		minValPsiShear=0.8):
	# extract needed data
	a0Val = bf.getDictValOrDef(data, 'a0Val')
	psiUni = data['psiUni']
	sin2psiUni = data['sin2psiUni']
	sinpsi2Uni = data['sinpsi2Uni']
	psiVals = data['psiVals']
	phiVals = data['phiVals']
	dVals = data['dVals']
	dErrVals = data['dErrVals']
	tauVals = data['tauVals']
	hklVal = data['hklVal']
	s1Val = data['s1Val']
	hs2Val = data['hs2Val']
	phi4 = data['phi4']
	# define needed variables
	sinpsi2Star = -2 * s1Val / hs2Val
	psiStar = bc.asind(sinpsi2Star ** 0.5)
	valsAll = np.zeros((len(psiUni), 3))
	fplus = np.zeros((len(psiUni), 3))
	fminus = np.zeros((len(psiUni), 3))
	f13 = np.zeros((len(psiUni), 3))
	f23 = np.zeros((len(psiUni), 3))
	stresses = np.zeros((len(psiUni), 4))
	errVals = np.zeros((len(psiUni), 4))
	tauRes = np.zeros(len(psiUni))
	resData = dict()
	validCounter = 0
	for i in range(len(psiUni)):
		if phi4:
			cond0 = (psiVals == psiUni[i]) & (phiVals == 0)
			cond90 = (psiVals == psiUni[i]) & (phiVals == 90)
			cond180 = (psiVals == psiUni[i]) & (phiVals == 180)
			cond270 = (psiVals == psiUni[i]) & (phiVals == 270)
		else:
			cond0 = (psiVals == psiUni[i]) & (phiVals == 0)
			cond90 = (psiVals == psiUni[i]) & (phiVals == 90)
			cond180 = (psiVals == -psiUni[i]) & (phiVals == 0)
			cond270 = (psiVals == -psiUni[i]) & (phiVals == 90)
		val0 = np.concatenate((dVals[cond0], dVals[cond0] - dErrVals[cond0],
			dVals[cond0] + dErrVals[cond0]))
		val90 = np.concatenate((dVals[cond90], dVals[cond90] - dErrVals[cond90],
			dVals[cond90] + dErrVals[cond90]))
		val180 = np.concatenate((dVals[cond180], dVals[cond180] - dErrVals[cond180],
			dVals[cond180] + dErrVals[cond180]))
		val270 = np.concatenate((dVals[cond270], dVals[cond270] - dErrVals[cond270],
			dVals[cond270] + dErrVals[cond270]))
		if len(val0) > 0 and len(val90) > 0 and len(val180) > 0 and len(val270) > 0:
			fplus[i, :] = val0 + val90 + val180 + val270
			fminus[i, :] = (val0 + val180) - (val90 + val270)
			f13[i, :] = val0 - val180
			f23[i, :] = val90 - val270
			valsAll[i, :] = 0.25 * fplus[i, :]
			tauRes[i] = np.mean(tauVals[psiVals == psiUni[i]])
			validCounter += 1
	# perform linear regression for all phi values to get dstar
	used = valsAll[:, 0] != 0
	if len(valsAll[used, 0]) > 1:
		if maxPsi is not None:
			usedReg = (sinpsi2Uni <= bc.sind(maxPsi) ** 2) & used
		else:
			usedReg = used
		[m, b, CoD, yD, err, mErr, bErr] = dm.linRegWeighted(sinpsi2Uni[usedReg], valsAll[usedReg, 0],
			1 / valsAll[usedReg, 1])
		dStarVal = m * sinpsi2Star + b
		# calculate fplus
		denom = hs2Val * sinpsi2Uni[used] + 2 * s1Val
		fplus[used, 0] = (0.25 * fplus[used, 0] / dStarVal - 1) / denom
		fplus[used, 1] = (0.25 * fplus[used, 1] / dStarVal - 1) / denom
		fplus[used, 2] = (0.25 * fplus[used, 2] / dStarVal - 1) / denom
		# indicate or correct invalid values
		if minDistPsiStar is not None:
			invalidVals = used & (np.abs(sinpsi2Uni - sinpsi2Star) <= minDistPsiStar)  # around psiStar
			fplus[invalidVals, 0] = np.NAN
			fplus[invalidVals, 1] = np.NAN
			fplus[invalidVals, 2] = np.NAN
		# fplus[used, 0] = (0.25 * fplus[used, 0] / dStarVal - 1) / (np.sign(denom) * np.array(
		# 	[np.max((np.abs(i), 5e-7)) for i in denom]))  # constrained denominator
		# fplus[used, 1] = (0.25 * fplus[used, 1] / dStarVal - 1) / (np.sign(denom) * np.array(
		# 	[np.max((np.abs(i), 5e-7)) for i in denom]))  # constrained denominator
		# fplus[used, 2] = (0.25 * fplus[used, 2] / dStarVal - 1) / (np.sign(denom) * np.array(
		# 	[np.max((np.abs(i), 5e-7)) for i in denom]))  # constrained denominator
		# calculate fminus
		fminus[used, 0] = 0.25 * (fminus[used, 0] / dStarVal) / (hs2Val * sinpsi2Uni[used])
		fminus[used, 1] = 0.25 * (fminus[used, 1] / dStarVal) / (hs2Val * sinpsi2Uni[used])
		fminus[used, 2] = 0.25 * (fminus[used, 2] / dStarVal) / (hs2Val * sinpsi2Uni[used])
		# indicate invalid values
		if minValPsiNormal is not None:
			invalidVals = used & (sinpsi2Uni <= minValPsiNormal)  # small psi values
			fminus[invalidVals, 0] = np.NAN
			fminus[invalidVals, 1] = np.NAN
			fminus[invalidVals, 2] = np.NAN
		# calculate f13 and f23
		f13[used, 0] = 0.5 * (f13[used, 0] / dStarVal) / (hs2Val * sin2psiUni[used])
		f13[used, 1] = 0.5 * (f13[used, 1] / dStarVal) / (hs2Val * sin2psiUni[used])
		f13[used, 2] = 0.5 * (f13[used, 2] / dStarVal) / (hs2Val * sin2psiUni[used])
		f23[used, 0] = 0.5 * (f23[used, 0] / dStarVal) / (hs2Val * sin2psiUni[used])
		f23[used, 1] = 0.5 * (f23[used, 1] / dStarVal) / (hs2Val * sin2psiUni[used])
		f23[used, 2] = 0.5 * (f23[used, 2] / dStarVal) / (hs2Val * sin2psiUni[used])
		# indicate invalid values
		if minValPsiShear is not None:
			invalidVals = used & (sin2psiUni <= minValPsiShear)  # small and high psi values
			f13[invalidVals, 0] = np.NAN
			f13[invalidVals, 1] = np.NAN
			f13[invalidVals, 2] = np.NAN
			f23[invalidVals, 0] = np.NAN
			f23[invalidVals, 1] = np.NAN
			f23[invalidVals, 2] = np.NAN
		# determine stresses
		stresses[used, 0] = fplus[used, 0] + fminus[used, 0]
		stresses[used, 1] = fplus[used, 0] - fminus[used, 0]
		stresses[used, 2] = f13[used, 0]
		stresses[used, 3] = f23[used, 0]
		# determine error values
		errVals[used, 0] = np.abs((fplus[used, 2] + fminus[used, 2]) - (fplus[used, 1] + fminus[used, 1])) / 2
		errVals[used, 1] = np.abs((fplus[used, 2] - fminus[used, 2]) - (fplus[used, 1] - fminus[used, 1])) / 2
		errVals[used, 2] = np.abs(f13[used, 2] - f13[used, 1]) / 2
		errVals[used, 3] = np.abs(f23[used, 2] - f23[used, 1]) / 2
		# also determine s33
		minVal, minPos = bf.min(abs(psiVals - psiStar))
		tauS33 = np.mean(tauVals[minPos])  # tau equivalent to condition at psiStar
		#leftVal, leftPos = bf.max(psiVals[psiVals < psiStar])
		#rightVal, rightPos = bf.min(psiVals[psiVals > psiStar])
		#tauS33 = np.interp(psiStar, [leftVal, rightVal], [tauVals[leftPos], tauVals[rightPos]])
		aStarVal = conv.latticeDists2aVals2(dStarVal, hklVal)
		if a0Val is None or a0Val == 0:
			s33 = 0
			ds33 = 0
		else:
			s33 = calcSigma33(aStarVal, a0Val, s1Val, hs2Val)
			ds33 = 2 * mErr ** 0.5 / (hs2Val * dStarVal)
		resData = {'tauRes': tauRes, 'stresses': stresses, 'errVals': errVals, 'validCounter': validCounter,
			'tauS33': tauS33, 'dStar100': aStarVal, 'dStar100Err': 2 * bErr ** 0.5, 's33': s33, 'dev_s33': ds33}
	return resData


def multiUniversalPlotAnalysis(data, maxPsi=None, minDistPsiStar=0.15,
		minValPsiNormal=0.08, minValPsiShear=0.8):
	keyList = bf.getKeyList(data)
	peakCount = int(bf.replace(keyList[-1].split('_')[0], 'pv')) + 1  # must be adapted in further versions!!!
	tthVal = data['tth']
	phiVals = data['phi']
	psiVals = data['psi']
	psiUni = np.unique(psiVals)
	psiUni = psiUni[psiUni != 0]  # no zero value
	psiSign = np.sign(psiUni[-1])  # sign of last psi value
	psiUni = psiUni[np.sign(psiUni) == psiSign]  # only negative or positive values
	sinpsi2Uni = bc.sind(psiUni) ** 2
	sin2psiUni = bc.sind(np.abs(2 * psiUni))
	tauRes = np.zeros((peakCount, len(psiUni)))
	hklRes = np.zeros((peakCount, len(psiUni)))
	psiRes = np.zeros((peakCount, len(psiUni)))
	stresses = np.zeros((peakCount, len(psiUni), 4))
	errVals = np.zeros((peakCount, len(psiUni), 4))
	tauS33 = np.zeros(peakCount)
	aStarVals = np.zeros(peakCount)
	aStarErrVals = np.zeros(peakCount)
	s33 = np.zeros(peakCount)
	dev_s33 = np.zeros(peakCount)
	hklList = np.zeros(peakCount)
	phi4 = len(np.unique(phiVals)) == 4
	validCounter = 0
	for p in range(peakCount):  # for all peaks create one plot
		centerVals = data['pv' + str(p) + '_center']
		centerErrVals = data['pv' + str(p) + '_center_err']
		dMinVals = conv.energies2latticeDists(centerVals - centerErrVals, tthVal)
		dMaxVals = conv.energies2latticeDists(centerVals + centerErrVals, tthVal)
		dErrVals = np.abs(dMaxVals - dMinVals) / 2
		tauVals = data['pv' + str(p) + '_depth']
		dVals = data['pv' + str(p) + '_dspac'] / 10  # in nm
		#hklVals = data['pv' + str(p) + '_hklList']  # first version with adaptions
		#s1Vals = data['pv' + str(p) + '_s1List']  # first version with adaptions
		#hs2Vals = data['pv' + str(p) + '_s2List']  # first version with adaptions
		#hklVal = hklVals[0]
		hVals = data['pv' + str(p) + '_h']  # second version
		kVals = data['pv' + str(p) + '_k']  # second version
		lVals = data['pv' + str(p) + '_l']  # second version
		s1Vals = data['pv' + str(p) + '_s1']  # second version
		#hs2Vals = data['pv' + str(p) + '_s2']  # second version
		hs2Vals = data['pv' + str(p) + '_hs2']  # third version
		hklVal = conv.mergeHkl(hVals[0], kVals[0], lVals[0])
		hklList[p] = hklVal
		hklRes[p] = hklVal * np.ones(len(psiUni))
		psiRes[p] = psiUni
		s1Val = s1Vals[0]
		#hs2Val = hs2Vals[0] * 0.5  # test valid for first and second version!!!!!
		hs2Val = hs2Vals[0]
		curData = {'tauVals': tauVals, 'dVals': dVals, 'dErrVals': dErrVals, 'psiVals': psiVals,
			'phiVals': phiVals, 'psiUni': psiUni, 'sin2psiUni': sin2psiUni, 'sinpsi2Uni': sinpsi2Uni, 'phi4': phi4,
			'hklVal': hklVal, 's1Val': s1Val, 'hs2Val': hs2Val}
		bf.extendDictionary(curData, data, ('a0Val',))
		# perform universal plot analysis for current peak data
		curResData = universalPlotAnalysis(curData, maxPsi, minDistPsiStar,
			minValPsiNormal, minValPsiShear)
		# remember results
		tauRes[p] = curResData['tauRes']
		stresses[p] = curResData['stresses']
		errVals[p] = curResData['errVals']
		aStarVals[p] = curResData['dStar100']
		aStarErrVals[p] = curResData['dStar100Err']
		tauS33[p] = curResData['tauS33']
		s33[p] = curResData['s33']
		dev_s33[p] = curResData['dev_s33']
		validCounter += curResData['validCounter']
	# reshape data
	tauRes = np.reshape(tauRes, np.prod(bf.size(tauRes)))
	hklRes = np.reshape(hklRes, np.prod(bf.size(hklRes)))
	psiRes = np.reshape(psiRes, np.prod(bf.size(psiRes)))
	stresses = np.reshape(stresses, (bf.size(stresses, 0) * bf.size(stresses,1), bf.size(stresses, 2)))
	errVals = np.reshape(errVals, (bf.size(errVals, 0) * bf.size(errVals, 1), bf.size(errVals, 2)))
	# remove values with tau = 0
	hklRes = hklRes[tauRes > 0]
	psiRes = psiRes[tauRes > 0]
	stresses = stresses[tauRes > 0]
	errVals = errVals[tauRes > 0]
	tauRes = tauRes[tauRes > 0]
	# sort data concerning increasing information depth
	hklRes = hklRes[np.argsort(tauRes)]
	psiRes = psiRes[np.argsort(tauRes)]
	stresses = stresses[np.argsort(tauRes)]
	errVals = errVals[np.argsort(tauRes)]
	tauRes = tauRes[np.argsort(tauRes)]
	resData = {'tauVals': tauRes, 'stresses': stresses, 'accuracy': errVals, 'hklVals': hklRes,
		'psiVals': psiRes, 'validCount': validCounter}
	resDataS33 = {'tauMean': tauS33, 'dStar100': aStarVals, 'dStar100Err': aStarErrVals, 's33': s33,
		'dev_s33': dev_s33, 'hklList': hklList}
	return resData, resDataS33


def initStressAnalysisSettings(*args):
	# default values
	settings = {"method": "ED", "showPlots": True, "showDeviation": False, "anode": "", "material": "Fe", "a0Val": 0.28665,
		"decList": np.array([]), "maxPsi": 45}  # 'ED', 'AD'
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


# par contains the factors of polynom a0, a1, ... an, b
def stressDampedPolynomTau(par, x):
	y = 0
	for i in range(len(par) - 1):
		y = y + par[i] / (par[-1] + 1 / x)**(i + 1)
	y = y / x
	return y


# par contains the factors of polynom a0, a1, ... an, b
def stressDampedPolynomReal(par, x):
	y = 0
	for i in range(len(par) - 1):
		y = y + par[i] / (math.factorial(i) * x**i)
	y = y * np.exp(-par[-1] * x)
	return y


def stressCurve(damping, n, tauVals, stresses, zVals=None):
	if zVals is None:
		zVals = tauVals
	if bf.size(stresses, 0) > 1 and bf.size(stresses, 1) > 1:
		if bf.size(stresses, 0) == 2 and bf.size(stresses, 1) == bf.length(tauVals):
			stressVals = stresses[0, :]
			weigthVals = stresses[1, :]
		elif bf.size(stresses, 0) == bf.length(tauVals) and bf.size(stresses, 1) == 2:
			stressVals = stresses[:, 0]
			weigthVals = stresses[:, 1]
		# stress accuracy values used to weigthen the stress values
		[tauValsFit, stressVals] = bc.createWeightedData(tauVals, stressVals, weigthVals)
	else:
		tauValsFit = tauVals
		stressVals = stresses
	if damping:
		stressDampedPolynomTauOpt = lambda par, x, y: stressDampedPolynomTau(par, x) - y
		res = op.leastsq(stressDampedPolynomTauOpt, bf.ones(n), (tauValsFit, stressVals))
		aVals = res[0]
		err = res[1]
		tauStresses = stressDampedPolynomTau(aVals, tauVals)
		zStresses = stressDampedPolynomReal(aVals, zVals)
	else:
		aVals = np.polyfit(tauValsFit, stressVals, n)
		#err = res.normr
		tauStresses = np.polyval(aVals, tauVals)
		bVals = aVals
		for i in range(len(aVals),0,-1):
			bVals[i - 1] = bVals[i - 1] / math.factorial(i - 1)
		zStresses = np.polyval(bVals, zVals)
	return tauStresses, zStresses, aVals#, err


# provided by Guilherme and adapted
def hklGenerator(phase, latticePar, tth, energyRange=None):
	latticePar = np.multiply(latticePar, 10)  # adapt lattice parameter to match Angstroem unit (used in calculations)
	hklTable = np.array([])
	if phase == 'hcp' or phase == 'HCP' or phase == 'Hcp':
		if bf.length(latticePar) > 2:
			a = latticePar[0]
			c = latticePar[2]
		else:
			a = latticePar[0]
			c = latticePar[1]
		for h in bf.createRange(0, 1, 10):
			for k in bf.createRange(0, 1, h):
				for l in bf.createRange(0, 1, 10):
					# hkl: l = even or h+2k ~= 3n
					if (h + k + l) > 0 and ((l % 2) == 0 or ((h + 2 * k) % 3) * ((h + 2 * k) % 3) != 0):
						# dVal = latticeSpacings([a, c], [90, 120], h, k, l)
						dVal = (1 / ((4 / 3 * (h ** 2 + h * k + k ** 2) / a ** 2) + (l ** 2 / c ** 2))) ** 0.5
						energy = 12.39842 / (2 * dVal * bc.sind(tth / 2))
						if energyRange is not None:
							if energy > energyRange[0] and energy < energyRange[1]:
								if len(hklTable) > 0:
									hklTable = np.vstack([hklTable, [h, k, l, dVal, energy]])
								else:
									hklTable = np.hstack([hklTable, [h, k, l, dVal, energy]])
	elif phase == 'fcc' or phase == 'FCC' or phase == 'Fcc':
		if bf.length(latticePar) > 1:
			a = latticePar[0]
		else:
			a = latticePar
		for h in bf.createRange(0, 1, 10):
			for k in bf.createRange(0, 1, h):
				for l in bf.createRange(0, 1, k):
					# hkl: h,k,l = odd or h,k,l = even
					if (h + k + l) > 0 and ((h % 2) * (k % 2) * (l % 2) != 0 or (h % 2) + (k % 2) + (l % 2) == 0):
						dVal = conv.aVals2latticeDists(a, h, k, l)
						energy = 12.39842 / (2 * dVal * bc.sind(tth / 2))
						if energyRange is not None:
							if energy > energyRange[0] and energy < energyRange[1]:
								if len(hklTable) > 0:
									hklTable = np.vstack([hklTable, [h, k, l, dVal, energy]])
								else:
									hklTable = np.hstack([hklTable, [h, k, l, dVal, energy]])
	elif phase == 'bcc' or phase == 'BCC' or phase == 'Bcc':
		if bf.length(latticePar) > 1:
			a = latticePar[0]
		else:
			a = latticePar
		for h in bf.createRange(0, 1, 10):
			for k in bf.createRange(0, 1, h):
				for l in bf.createRange(0, 1, k):
					#hkl: h+k+l = even
					if (h + k + l) > 0 and (h + k + l) % 2 == 0:
						dVal = conv.aVals2latticeDists(a, h, k, l)
						energy = 12.39842 / (2 * dVal * bc.sind(tth / 2))
						if energyRange is not None:
							if energy > energyRange[0] and energy < energyRange[1]:
								if len(hklTable) > 0:
									hklTable = np.vstack([hklTable, [h, k, l, dVal, energy]])
								else:
									hklTable = np.hstack([hklTable, [h, k, l, dVal, energy]])
	# sort hkl table
	hklTable = hklTable[hklTable[:, 4].argsort()]
	return hklTable
