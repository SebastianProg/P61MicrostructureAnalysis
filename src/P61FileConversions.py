import numpy as np

from P61ANexusReader import P61ANexusReader
import basics.functions as bf
import filehandling.general as fg
import filehandling.specific as fs


# Nexus file conversion
P61ANexusReader().main()

# Standard conversion of fio files
# header for position files: om eta tth0 tth1 psi0 psi1 phi x y z petracurrent
header = ['om', 'eta', 'tth0', 'tth1', 'psi0', 'psi1', 'phi', 'x', 'y', 'z', 'petraCur']
headLine = bf.stringList2string(header, '\t')
# Fio file processing for residual stress analysis or texture analysis
# phioffset = 45
phioffset = 0
fileNames = fg.requestFiles((("Fio files", '*.fio'),), 'Select motor positions files', 'on')
for i in range(len(fileNames)):
    scanData, fixData = fs.read_fio(fileNames[i])
    motVals = np.zeros((scanData.shape[0], len(header)))
    # set constant values
    if 'eu.omg' in fixData.keys():
        motVals[:, 0] = fixData['eu.omg']
    if 'eu.eta' in fixData.keys():
        motVals[:, 1] = fixData['eu.eta']
    if 'eu.alp' in fixData.keys():
        motVals[:, 0] = fixData['eu.alp']
    if 'eu.bet' in fixData.keys():
        motVals[:, 1] = fixData['eu.bet']
    if 'eu.chi' in fixData.keys():
        motVals[:, 4] = -90 + fixData['eu.chi']
        motVals[:, 5] = fixData['eu.chi']
    if 'eu.phi' in fixData.keys():
        motVals[:, 6] = fixData['eu.phi'] - phioffset
    if 'eu.x' in fixData.keys():
        motVals[:, 7] = fixData['eu.x']
    if 'eu.y' in fixData.keys():
        motVals[:, 8] = fixData['eu.y']
    if 'eu.z' in fixData.keys():
        motVals[:, 9] = fixData['eu.z']
    if 'petracurrent' in fixData.keys():
        motVals[:, 10] = fixData['petracurrent']
    # set user determined values -> tth values
    # motVals[:, 2] = 15  # detector 0
    # motVals[:, 3] = 6.785  # detector 1
    # motVals[:, 2] = 9.72  # detector 0
    # motVals[:, 3] = 10.23  # detector 1
    # motVals[:, 3] = 6.225  # detector 1
    motVals[:, 2] = 8.14  # detector 0
    motVals[:, 3] = 5.1  # detector 1
    # set variable values
    if 'eu.omg' in scanData.columns:
        motVals[:, 0] = np.array(scanData['eu.omg'])
    if 'eu.eta' in scanData.columns:
        motVals[:, 1] = np.array(scanData['eu.eta'])
    if 'eu.alp' in scanData.columns:
        motVals[:, 0] = np.array(scanData['eu.alp'])
    if 'eu.bet' in scanData.columns:
        motVals[:, 1] = np.array(scanData['eu.bet'])
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
    if 'petracurrent' in scanData.columns:
        motVals[:, 10] = np.array(scanData['petracurrent'])
    # export new text file
    fg.dlmwrite(bf.replace(fileNames[i], '.fio', '.txt'), motVals, head=headLine, format='%f')
