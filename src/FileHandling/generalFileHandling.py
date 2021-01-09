# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:53:13 2019

@author: Sebastian
"""

import tkinter
from tkinter import filedialog
import numpy as np
import os
import nexusformat.nexus as nxs

import BasicFunctions.generalFunctions as gf


def fileparts(fullFile):
	pathName, file = os.path.split(fullFile)
	return pathName, file


def fileparts2(fullFile):
    filePath, ext = os.path.splitext(fullFile)
    return filePath, ext


# check if file or path exists
def exists(fileOrPath):
	return os.path.exists(fileOrPath)


# check if file exists
def existsFile(file):
	return os.path.isfile(file)


# check if directory exists
def existsDir(dir):
	return os.path.isdir(dir)


# Request files by user selection
def requestFiles(fileType=(("All files", "*.*"),), dialogTitle="", multiselect="on", folder="/"):
	if len(dialogTitle) == 0:
		if multiselect == "on":
			dialogTitle = "Select files"
		else:
			dialogTitle = "Select file"
	# request file(s)
	root = tkinter.Tk()
	root.withdraw()
	if multiselect == "on":
		return filedialog.askopenfilenames(initialdir=folder, title=dialogTitle, filetypes=fileType)
	else:
		return (filedialog.askopenfilename(initialdir=folder, title=dialogTitle, filetypes=fileType),)  # put single file into tupel


def requestFiles2(fileType=(("All files", "*.*"),), dialogTitle="", multiselect="on", folder="/"):
	files = requestFiles(fileType, dialogTitle, multiselect, folder)
	fileNames = list(files)
	for i in range(len(files)):
		pathName, fileName = fileparts(files[i])
		fileNames[i] = fileName
	return fileNames, pathName


def requestSaveFile(fileType=(("All files", "*.*"),), dialogTitle="", folder="/"):
	if len(dialogTitle) == 0:
		dialogTitle = "Select file to save"
	# request file(s)
	root = tkinter.Tk()
	root.withdraw()
	return filedialog.asksaveasfilename(initialdir=folder, title=dialogTitle, filetypes=fileType)


def requestDirectory(parentWindow=None, dialogTitle="Select a folder:", folder=os.getcwd()):
	return filedialog.askdirectory(parent=parentWindow,initialdir=folder,title=dialogTitle)


def writeLine(fileId, data):
	fileId.write(data + gf.newlineChar())


def writeData(fileId, formatStr, dataArray):
	fileId.write(fileId, (formatStr % tuple(dataArray)))


def writeDataLine(fileId, formatStr, dataArray):
	writeLine(fileId, (formatStr % tuple(dataArray)))


def dlmread(file, delim="\t", skipRows=0, usedCols=None, maxRows=None, commentStr="#", convert=None, dtypeVal=float):
	return np.loadtxt(file, dtype=dtypeVal, comments=commentStr, delimiter=delim, converters=convert, skiprows=skipRows,
		usecols=usedCols, max_rows=maxRows)


def dlmwrite(file, data, delim="\t", newln=gf.newlineChar(), head="", foot="", commentStr="", append="off", format="%.18e"):
	if append == "on":
		fid = open(file,'ab')
		np.savetxt(fid, data, fmt=format, delimiter=delim, newline=newln, header=head, footer=foot, comments=commentStr)
		fid.close()
	else:
		np.savetxt(file, data, fmt=format, delimiter=delim, newline=newln, header=head, footer=foot, comments=commentStr)


def readLines(fileName, lines=None, startLine=0):
	file = open(fileName, 'r')
	data = file.readlines()
	file.close()
	if lines is not None:
		relData = data[startLine:lines+startLine]
	else:
		relData = data
	# remove newlines or spaces
	relData = [relData[i].strip() for i in range(len(relData))]
	return relData


def writeList(file, listData, delim='\t'):
	# write list data to file
	fid = open(file, 'w')
	for i in range(gf.size(listData, 0)):
		if gf.size(listData[i], 0) > 1:
			fid.write(str(listData[i][0]))
			for j in range(1, gf.size(listData[i], 0)):
				fid.write(delim + str(listData[i][j]))
		elif gf.size(listData[i], 0) == 1:
			fid.write(str(listData[i]))
		else:
			fid.write('')
		writeLine(fid, '')
	fid.close()


def export(text, resFile=None):
	if resFile is None:
		resFile = requestSaveFile((("Text file", "*.txt"),), "Specify result file")
	fid = open(resFile, 'w')
	fid.write(text)
	fid.close()


def saveSettings(settingsDict, fileName=None):
	if fileName is None:
		fileName = requestSaveFile(fileType=(("Text files", "*.txt"), ("All files", "*.*")),
			dialogTitle="Select settings file", folder="/")
	fid = open(fileName, 'w')
	for item in settingsDict.items():
		parstr = str(item[1])
		if len(parstr) == 0:
			parstr = "''"
		if parstr.isidentifier():
			parstr = "'" + parstr + "'"
		elif parstr.startswith('[['):
			parstr = gf.matrix2str(item[1])
		elif parstr.startswith('['):
			parstr = '[' + gf.array2str(item[1], ', ') + ']'
		writeLine(fid, item[0] + '\t' + parstr)
	fid.close()
	return fileName


def readnxs(filename, show=False):
	nxsdata = nxs.nxload(filename)
	if show:
		print(nxsdata.tree)
	return nxsdata
