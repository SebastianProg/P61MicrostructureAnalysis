import numpy as np
import os
import sys
import re
from time import sleep, localtime, strftime
from importlib import reload


def init():
	# set amount of elements to print
	np.set_printoptions(threshold=sys.maxsize)


def update(lib):
	reload(lib)


def pythonPath():
	return os.path.dirname(sys.executable)


def newlineChar():
	return '\n'


def pause(secs):
	sleep(secs)


def loctime():
	return localtime()


def datetimestr(formatStr='%d.%m.%Y %H:%M:%S'):
	return strftime(formatStr, localtime())


def datestr(sep='.'):
	return datetimestr('%d' + sep + '%m' + sep + '%Y')


def timestr(sep=':'):
	return datetimestr('%H' + sep + '%M' + sep + '%S')


def datetimearr(formatArr=('%Y', '%m', '%d', '%H', '%M', '%S')):
	return [strftime(item, localtime()) for item in formatArr]


def size(x, dim=None):
	if dim is not None:
		shapeVals = np.shape(x)
		if len(shapeVals) > dim:
			return np.shape(x)[dim]
		else:
			return 1
	else:
		return np.array(np.shape(x))


def lenArray(vals):
	return np.array(vals).size


def length(x):
	return max(size(x, 0), size(x, 1))[0]


def createRange(start=0, step=0.1, stop=1):
	return np.arange(start, stop + 0.001 * step, step)


def flat(vals, lines=0):
	if lines == 0:
		flatVals = np.reshape(vals, lenArray(vals))
	else:
		flatVals = np.reshape(vals, (1, lenArray(vals)))
	return flatVals


def meshdata(xvals, yvals, zvals, cvals=None, mvals=None):
	x_uni = np.unique(xvals)
	y_uni = np.unique(yvals)
	# xout, yout = np.meshgrid(x_uni, y_uni)
	xout = np.reshape(xvals, (len(x_uni), len(y_uni)))
	yout = np.reshape(yvals, (len(x_uni), len(y_uni)))
	zout = np.reshape(zvals, (len(x_uni), len(y_uni)))
	if cvals is None and mvals is None:
		return xout, yout, zout
	else:
		if cvals is not None:
			cout = np.reshape(cvals, (len(x_uni), len(y_uni)))
		if mvals is not None:
			mout = np.reshape(mvals, (len(x_uni), len(y_uni)))
		if cvals is not None and mvals is not None:
			return xout, yout, zout, cout, mout
		elif cvals is not None:
			return xout, yout, zout, cout
		else:
			return xout, yout, zout, mout


def vectordata(m, xm=None, ym=None, sel=None, ind=False):
	if xm is None:
		if ind:
			xx = np.arange(m.shape[0])
		else:
			xx = np.arange(1, m.shape[0] + 1)
		xm = np.array(np.transpose(np.mat(xx)) * np.mat(np.ones(m.shape[1])))
	if ym is None:
		if ind:
			yy = np.arange(m.shape[1])
		else:
			yy = np.arange(1, m.shape[1] + 1)
		ym = np.array(np.transpose(np.mat(np.ones(m.shape[0]))) * np.mat(yy))
	if sel is not None:
		# use only selected values
		xm = xm[sel]
		ym = ym[sel]
		m = m[sel]
	# create result values
	return flat(xm), flat(ym), flat(m)


def vectordata2(m, xm=None, ym=None, sel=None, ind=False):
	xout, yout, zout = vectordata(m, xm, ym, sel, ind)
	return np.array([xout, yout, zout])


def strrepAll(origStr, oldSubstr, newSubstr):
	strVal = origStr
	for i in range(len(oldSubstr)):
		return strVal.replace(oldSubstr[i], newSubstr)


def replace(source, old, new=""):
	return source.replace(old, new)


def replaceAll(source, old, new=""):
	for c in old:
		source = replace(source, c, new)
	return source


def replaceEach(source, oldItems, newItems):
	for i in range(len(oldItems)):
		if len(newItems) > i:
			source = replace(source, oldItems[i], newItems[i])
		else:
			source = replace(source, oldItems[i])
	return source


def replaceInList(source, old, new='', copy=False):
	if copy:
		resList = source.copy()
	else:
		resList = source
	for i in range(len(resList)):
		resList[i] = replace(resList[i], old, new)
	return resList


def replaceAllInList(source, old, new='', copy=False):
	if copy:
		resList = source.copy()
	else:
		resList = source
	for i in range(len(resList)):
		resList[i] = replaceAll(resList[i], old, new)
	return resList


def replaceEachInList(source, oldItems, newItems, copy=False):
	if copy:
		resList = source.copy()
	else:
		resList = source
	for i in range(len(resList)):
		resList[i] = replaceEach(resList[i], oldItems, newItems)
	return resList


def num2strNoSpecialChar(number):
	numText = str(number)
	return replace(replace(numText, '.', 'p'), '-', 'm')


def numsInString(strVal, numlen=0, selPos=0, negAllowed=False):
	if negAllowed:
		pattern = r"[-]?\d"
	else:
		pattern = r"\d"
	if numlen < 1:
		pattern += "+"
	else:
		pattern += "{" + str(numlen) + "}"
	numsStr = re.findall(pattern, strVal)
	# for m in re.finditer(r"\d+", strVal):
	# 	print('%02d-%02d: %s' % (m.start(), m.end(), m.group(0)))
	if len(numsStr) > 0:
		nums = np.array([int(x) for x in numsStr])
		selNum = nums[min(selPos, len(nums)-1)[0]]
	else:
		nums = np.array([])
		selNum = np.nan
	return nums, numsStr, selNum


def fnumsInString(strVal, numlen=0, selPos=0, negAllowed=True):
	patternEnd = r"(\.\d*)?|\.\d+)([eE][-+]?\d+)?"
	if negAllowed:
		patternStart = r"[-]?(\d"  # r"[-+]?(\d"
	else:
		patternStart = r"(\d"
	if numlen < 1:
		patternLen = "+"
	else:
		patternLen = "{" + str(numlen) + "}"
	numsStr = re.findall(patternStart + patternLen + patternEnd, strVal)
	if len(numsStr) > 0:
		# numsStr = np.array(numsStr)
		# numsStr = numsStr[:,0]
		nums = np.array([float(x[0]) for x in numsStr])
		selNum = nums[min(selPos, len(nums)-1)[0]]
	else:
		nums = np.array([])
		selNum = np.nan
	return nums, numsStr, selNum


def split(source, sep=' ', remEmpty=False, max=-1):
	vals = source.split(sep,max)
	if remEmpty:
		return [x for x in vals if x]
	else:
		return vals


def isempty(value):
	return len(value) == 0


def min(values, secondArg=None):
	if secondArg is not None:
		values = [values, secondArg]
	return np.min(values), np.argmin(values)


def max(values, secondArg=None):
	if secondArg is not None:
		values = [values, secondArg]
	return np.max(values), np.argmax(values)


def zeros(dim0, dim1=None, dim2=None, dim3=None, dim4=None, dtypeVal=np.float64):
	if dim1 is None and dim2 is None and dim3 is None and dim4 is None:
		return np.zeros((dim0,), dtype=dtypeVal)
	elif dim2 is None and dim3 is None and dim4 is None:
		return np.zeros((dim0, dim1), dtype=dtypeVal)
	elif dim3 is None and dim4 is None:
		return np.zeros((dim0, dim1, dim2), dtype=dtypeVal)
	elif dim4 is None:
		return np.zeros((dim0, dim1, dim2, dim3), dtype=dtypeVal)
	else:
		return np.zeros((dim0, dim1, dim2, dim3, dim4), dtype=dtypeVal)


def ones(dim0, dim1=None, dim2=None, dim3=None, dim4=None):
	if dim1 is None and dim2 is None and dim3 is None and dim4 is None:
		return np.ones((dim0,))
	elif dim2 is None and dim3 is None and dim4 is None:
		return np.ones((dim0, dim1))
	elif dim3 is None and dim4 is None:
		return np.ones((dim0, dim1, dim2))
	elif dim4 is None:
		return np.ones((dim0, dim1, dim2, dim3))
	else:
		return np.ones((dim0, dim1, dim2, dim3, dim4))


def containsItems(list, members):
	result = np.isin(list, members)
	return result, np.where(result)[0]


def hasKey(dict, key):
	return key in dict


def getDictValOrDef(dict, key, defVal=None):
	if hasKey(dict, key):
		return dict[key]
	else:
		return defVal


def getKeyList(dict):
	return list(dict.keys())


def getValueList(dict):
	return list(dict.values())


def combineListsOrTuples(items1, items2):
	return items1 + items2


def setdict2(dict, key1, key2, value):
	if key1 in dict.keys():
		dict[key1][key2] = value
	else:
		dict[key1] = {key2: value}


def combineDictionaries(dict1, dict2):
	extendDictionary(dict1, dict2)
	return dict1


def extendDictionary(dict1, dict2, selectedKeys=None):
	if selectedKeys is None:
		dict1.update(dict2)
	else:
		for selKey in selectedKeys:
			addCheckedKeyValue(dict1, selKey, getDictValOrDef(dict2, selKey))


def addCheckedKeyValue(dict, key, value, excludeVal=None):
	if value is not excludeVal:
		dict[key] = value
	return dict


def dict2matrix(dict, keyList=None):
	if keyList is None:
		return np.array(getValueList(dict))
	else:
		matrix = np.zeros((len(dict[keyList[0]]), len(keyList)))
		for i in range(len(keyList)):
			matrix[:, i] = dict[keyList[i]]
		return matrix


def matrix2dict(matrix, keyList):
	dictMatrix = dict()
	if size(matrix, 0) == len(keyList):
		for i in range(len(keyList)):
			dictMatrix[keyList[i]] = matrix[i, :]
	elif size(matrix, 1) == len(keyList):
		for i in range(len(keyList)):
			dictMatrix[keyList[i]] = matrix[:, i]
	return dictMatrix


def stringList2string(strList, delim='', inplace=False):
	if inplace:
		resList = strList
	else:
		resList = strList.copy()
	if len(delim) > 0:
		[resList.insert(i*2, delim) for i in range(len(resList))]
		return ''.join(resList[1:])
	else:
		return ''.join(resList)


def array2str(array, sep=' '):
	res = ''
	for item in array:
		res = res + str(item) + sep
	res = res[0:-len(sep)]
	return res


def matrix2str(matrix):
	matrix = np.array(matrix)
	matstr = '['
	for i in range(matrix.shape[0]):
		linestr = '['
		for j in range(matrix.shape[1]):
			linestr = linestr + str(matrix[i, j]) + ', '
		matstr = matstr + linestr[0:-2] + '], '
	matstr = matstr[0:-2] + ']'
	return matstr


def strArray2npArray(strArray):
	return np.array(eval(strArray))


def toInt(stringOrValList):
	return [int(val) for val in stringOrValList]


def toFloat(stringOrValList):
	return [float(val) for val in stringOrValList]


def toString(numList):
	return [str(numVal) for numVal in numList]


def getTextPosition(texts, searchText):
	return texts.index(searchText)


def getTextPositions(texts, searchTexts):
	textPositions = np.zeros(size(searchTexts))
	for i in range(len(searchTexts)):
		textPositions[i] = getTextPosition(texts, searchTexts[i])
	return textPositions


def indexOf(condition):
	usedIndex = np.array(range(len(condition)))
	return usedIndex[condition]


def first(data):
	return data[0]


def last(data):
	return data[-1]
