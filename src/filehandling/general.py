# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 22:53:13 2019

@author: Sebastian
"""

import tkinter
from tkinter import filedialog
import numpy as np
import os
import sys
import errno

import basics.functions as bf


CUR_FOLDER = os.getcwd()  # this is the folder used as initial for file and folder dialogs


def fileparts(fullFile):
	pathName, file = os.path.split(fullFile)
	return pathName, file


def fileparts2(fullFile):
	filePath, ext = os.path.splitext(fullFile)
	return filePath, ext


def fileparts3(fullFile):
	filePath, ext = os.path.splitext(fullFile)
	pathName, file = os.path.split(filePath)
	return pathName, file, ext


# check if file or path exists
def exists(fileOrPath):
	return os.path.exists(fileOrPath)


# check if file exists
def existsFile(file):
	return os.path.isfile(file)


# check if directory exists
def existsDir(dir):
	return os.path.isdir(dir)


# request files by user selection
def requestFiles(fileType=(("All files", "*.*"),), dialogTitle="", multiselect="on", folder=None):
	global CUR_FOLDER
	if folder is None:
		folder = CUR_FOLDER
	if len(dialogTitle) == 0:
		if multiselect == "on":
			dialogTitle = "Select files"
		else:
			dialogTitle = "Select file"
	# request file(s)
	root = tkinter.Tk()
	root.wm_attributes('-topmost', True)  # show dialog in front of all
	root.withdraw()
	if multiselect == "on":
		files = filedialog.askopenfilenames(initialdir=folder, title=dialogTitle, filetypes=fileType, parent=root)
	else:
		# put single file into tuple
		files = (filedialog.askopenfilename(initialdir=folder, title=dialogTitle, filetypes=fileType, parent=root),)
	# get path of first selected file to set as new working directory
	if len(files) > 0 and len(files[0]) > 0:
		CUR_FOLDER, file = fileparts(files[0])
	return files


def requestFiles2(fileType=(("All files", "*.*"),), dialogTitle="", multiselect="on", folder=None):
	files = requestFiles(fileType, dialogTitle, multiselect, folder)
	fileNames = list(files)
	for i in range(len(files)):
		pathName, fileName = fileparts(files[i])
		fileNames[i] = fileName
	return fileNames, pathName


def requestSaveFile(fileType=(("All files", "*.*"),), dialogTitle="", folder=None):
	global CUR_FOLDER
	if folder is None:
		folder = CUR_FOLDER
	if len(dialogTitle) == 0:
		dialogTitle = "Select file to save"
	# request file(s)
	root = tkinter.Tk()
	root.withdraw()
	file = filedialog.asksaveasfilename(initialdir=folder, title=dialogTitle, filetypes=fileType)
	# get path of specified file to set as new working directory
	if len(file) > 0:
		CUR_FOLDER, _ = fileparts(file)
	return file


def requestDirectory(parentWindow=None, dialogTitle="Select a folder:", folder=None):
	global CUR_FOLDER
	if folder is None:
		folder = CUR_FOLDER
	pathName = filedialog.askdirectory(parent=parentWindow, initialdir=folder, title=dialogTitle)
	# set selected path as new working directory
	if len(pathName) > 0:
		CUR_FOLDER = pathName
	return pathName


def openFile(file, mode='r', buffering=-1, encoding=None, errors=None, newline=None, closefd=True, opener=None):
	if not os.path.exists(os.path.dirname(file)):
		try:
			os.makedirs(os.path.dirname(file))
		except OSError as exc:  # guard against race condition
			if exc.errno != errno.EEXIST:
				raise
	# with open(file, "w") as f:
	# 	f.write("TEST")
	return open(file, mode, buffering, encoding, errors, newline, closefd, opener)


def writeLine(fileId, data):
	fileId.write(data + bf.newlineChar())


def writeData(fileId, formatStr, dataArray):
	fileId.write(fileId, (formatStr % tuple(dataArray)))


def writeDataLine(fileId, formatStr, dataArray):
	writeLine(fileId, (formatStr % tuple(dataArray)))


def dlmread(file, delim="\t", skipRows=0, usedCols=None, maxRows=None, commentStr="#", convert=None, dtypeVal=float):
	return np.loadtxt(file, dtype=dtypeVal, comments=commentStr, delimiter=delim, converters=convert, skiprows=skipRows,
		usecols=usedCols, max_rows=maxRows)


def dlmwrite(file, data, delim="\t", newln=bf.newlineChar(), head="", foot="", commentStr="", append="off", format="%.5g"):
	if append == "on":
		fid = open(file, 'ab')
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
	for i in range(bf.size(listData, 0)):
		if bf.size(listData[i], 0) > 1:
			fid.write(str(listData[i][0]))
			for j in range(1, bf.size(listData[i], 0)):
				fid.write(delim + str(listData[i][j]))
		elif bf.size(listData[i], 0) == 1:
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
			parstr = bf.matrix2str(item[1])
		elif parstr.startswith('['):
			parstr = '[' + bf.array2str(item[1], ', ') + ']'
		writeLine(fid, item[0] + '\t' + parstr)
	fid.close()
	return fileName
