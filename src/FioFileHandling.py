import numpy as np
import itertools

import basics.functions as bf
import filehandling.general as fg
import filehandling.specific as fs


# Fio file processing for ROI analysis
fileNames = fg.requestFiles((("Fio files", '*.fio'),), 'Motorpositionsdateien auswaehlen', 'on')
for i in range(len(fileNames)):
    fioData = fs.read_fio(fileNames[i])
    fg.dlmwrite(bf.replace(fileNames[i], '.fio', '_roi.txt'), fioData.values,
        head=bf.stringList2string(list(fioData.columns), '\t'), format='%f')

# header for position files: om eta tth0 tth1 psi0 psi1 phi x y z
header = ['om', 'eta', 'tth0', 'tth1', 'psi0', 'psi1', 'phi', 'x', 'y', 'z']
headLine = bf.stringList2string(header, '\t')

# Fio file processing for residual stress analysis
# phioffset = 45
phioffset = 0
fileNames = fg.requestFiles((("Fio files", '*.fio'),), 'Motorpositionsdateien auswaehlen', 'on')
for i in range(len(fileNames)):
    scanData, fixData = fs.read_fio(fileNames[i])
    motVals = np.zeros((scanData.shape[0], 10))
    if 'eu.omg' in scanData.columns:
        motVals[:, 0] = np.array(scanData['eu.omg'])
    if 'eu.eta' in scanData.columns:
        motVals[:, 1] = np.array(scanData['eu.eta'])
    # motVals[:, 2] = 15
    # motVals[:, 3] = 6.785
    motVals[:, 2] = 9.72
    motVals[:, 3] = 10.23
    if 'eu.chi' in scanData.columns:
        motVals[:, 4] = -90 + np.array(scanData['eu.chi'])
        motVals[:, 5] = np.array(scanData['eu.chi'])
    if 'eu.phi' in scanData.columns:
        motVals[:, 6] = np.array(scanData['eu.phi']) - phioffset
    if 'eu.x' in scanData.columns:
        motVals[:, 7] = np.array(scanData['eu.x'])
    if 'eu.y' in scanData.columns:
        motVals[:, 8] = np.array(scanData['eu.y'])
    if 'eu.z' in scanData.columns:
        motVals[:, 9] = np.array(scanData['eu.z'])
    fg.dlmwrite(bf.replace(fileNames[i], '.fio', '.txt'), motVals, head=headLine, format='%f')

# split motor positions files into separate ones according to one axis (e. g. different measurement points)
# splitCol = 2  # tth0
# splitCol = 3  # tth1
splitCol = 7  # x
# splitCol = 8  # y
# splitCol = 9  # z
fileNames = fg.requestFiles((("Text files", '*.txt'),), 'Select positions files', 'on')
for curFile in fileNames:
    motposData = fg.dlmread(curFile, '\t', 1)
    uniVals = np.unique(motposData[:, splitCol])
    for val in uniVals:
        fg.dlmwrite(bf.replace(curFile, '.txt', '_' + header[splitCol] + '_' + str(val) + '.txt'),
            motposData[motposData[:, splitCol] == val, :], head=headLine, format='%f')

# split motor positions files into separate ones according to three axis (e. g. different measurement points)
splitCol1 = 7  # x
splitCol2 = 8  # y
splitCol3 = 9  # z
fileNames = fg.requestFiles((("Text files", '*.txt'),), 'Select positions files', 'on')
for curFile in fileNames:
    motposData = fg.dlmread(curFile, '\t', 1)
    uniVals1 = np.unique(motposData[:, splitCol1])
    uniVals2 = np.unique(motposData[:, splitCol2])
    uniVals3 = np.unique(motposData[:, splitCol3])
    for val1 in uniVals1:
        str1 = '_' + header[splitCol1] + '_' + str(val1)
        for val2 in uniVals2:
            str2 = '_' + header[splitCol2] + '_' + str(val2)
            for val3 in uniVals3:
                str3 = '_' + header[splitCol3] + '_' + str(val3)
                fg.dlmwrite(bf.replace(curFile, '.txt', str1 + str2 + str3 + '.txt'),
                    motposData[(motposData[:, splitCol1] == val1) & (motposData[:, splitCol2] == val2) &
                               (motposData[:, splitCol3] == val3), :], head=headLine, format='%f')
