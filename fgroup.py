import os
from os import path, listdir, access, getcwd
from os.path import isfile, isdir, normpath, basename, splitext, commonprefix
import sys
import numpy as np

#def getargs():
#	if len(sys.argv)>1:
#		return sys.argv[1:]
#	else:
#		return [getcwd()]
def hasriffhead(path):
	f=open(path,"rt")
	ascii4c = f.read(4)
	return ascii4c=="RIFF"
	
def fname(path):
	return splitext(basename(path))[0]
	
############################################################
# Group wave audio files by common prefix
#
# - returns two lists: outPaths and outNames
# - outPaths is a list of groups of paths to files with
#   similar prefixes
# - outNames is a list of the prefixes for file groups in
#   outPaths
#
# outPaths = [["C:\l1.wav", "C:\l2.wav"], ["C:\r1.wav", "C:\r2.wav"]]
# outNames = ["l", "r"]
#
############################################################
def fgrp(pth=None, debug=False):
	if pth==None:
		pth=getcwd()
	pth=normpath(pth)
	if not isdir(pth):
		print "Not a path!"
		return
	files = listdir(pth)
	if debug:
		print "name: isfile, isdir"
	
	wavs = []
	fnames = []
	for arg in files:
		fpath=path.join(pth,arg)
		if debug:	
			print arg,": ",isfile(fpath),", ",isdir(fpath)
		
		if isfile(fpath) and access(fpath, os.R_OK) and hasriffhead(fpath):
			wavs.append((fpath,fname(fpath)))
			fnames.append(fname(fpath))
	
	outNames= rec_pickgroups(rec_mergedown(rec_divide(fnames,0),0))
	outPaths=outNames
	for pth,fn in wavs:
		print "Including wave file ",fn," (",pth,")"
		outPaths = rec_replace(outPaths,fn,pth)
	
	
	
	for idx,arr in enumerate(outNames):
		outNames[idx]=commonprefix(arr)
	
	return outPaths, outNames
	
	
	return wavs
def unique(l):
	return list(set(l))
def rec_divide(inputArr,depth=0):
	searchArr = []
	skipped=0
	for input in inputArr:
		if(len(input)>depth):
			searchArr.append(input[depth])
		else:
			skipped+=1		
	searchArr=unique(searchArr)
	if len(searchArr)==len(inputArr)-skipped:
		return inputArr #everything at this depth is unique already, return
	outArr=[]
	for idx,uq in enumerate(searchArr):
		dArr = []
		for input in inputArr:
			if(len(input)>depth) and input[depth]==uq:
				dArr.append(input)
		outArr.append(dArr)
	for idx,arr in enumerate(outArr):
		outArr[idx]=rec_divide(arr,depth+1)
	return outArr

def rec_replace(inputArr, original, replacement):
	if isinstance(inputArr, basestring):
		if inputArr==original:
			return replacement
		return original
	outputArr=[]
	for arr in inputArr:
		if isinstance(arr, basestring):
			if arr==original:
				outputArr.append(replacement)
			else:
				outputArr.append(arr)
		else:
			outputArr.append(rec_replace(arr,original,replacement))
	return outputArr

def rec_mergeup(inArr, levels):
	if not isinstance(inArr, basestring):
		itArr=inArr
		for i in range(levels):
			outArr=[]
			if len(itArr)>0:
				iArr=[]
				for arr in itArr:
					if isinstance(arr,basestring):
						outArr.append(arr)
					else:
						if hasStrs(arr):
							iArr.extend(arr)
						else:
							iArr.append(rec_mergeup(arr,1))
				outArr.extend(iArr)
			itArr=outArr
		return itArr
	else:
		return inArr
		
def rec_mergedown(inArr, levels):
	if not isinstance(inArr, basestring):
		itArr=inArr
		for i in range(levels):
			outArr=[]
			if len(itArr)>0:
				iArr=[]
				for arr in itArr:
					if isinstance(arr,basestring):
						iArr.append(arr)
					else:
						if onlyStrs(arr):
							iArr.extend(arr)
						else:
							iArr.append(rec_mergedown(arr,1))
				outArr.extend(iArr)
			itArr=outArr
		return itArr
	else:
		return inArr
def rec_pickgroups(inArr):
	outArr=[]
	if len(inArr)>0 and not isinstance(inArr,basestring):
		strArr = []
		for arr in inArr:
			if not isinstance(arr,basestring):
				outArr.extend(rec_pickGroups(arr))
			else:
				#for item in arr:
				#	if isinstance(item,basestring):
				strArr.append(arr)
		if len(strArr)>0:
			outArr.append(strArr)
	else:
		return inArr
	return outArr
def hasStrs(arr):
	if len(arr)>0:
		for a in arr:
			if len(a)>0 and isinstance(a,basestring):
				return True
	return False
def onlyStrs(arr):
	if len(arr)>0:
		for a in arr:
			if len(a)>0 and not isinstance(a,basestring):
				return False
	return True
	
	
#BEGIN OBSOLETE----------------------------------
def rec_split(inArr):
	vlen = 0
	outArr = []
	if len(inArr)>1: #two or more subarrays
		vlen=len(inArr[0])
		for arr in inArr:
			if not len(arr)==vlen or isinstance(arr,basestring):
				outArr.append(arr)
			else:
				outArr.extend(rec_split(arr))
	elif len(inArr)==1 and not isinstance(inArr,basestring):
		return rec_split(inArr)
	else:
		return inArr
	return outArr
	
#flatten to a set depth
def rec_dflatten(inArr, maxdepth=1):
	outArr=[]
	if(maxdepth==0):
		return rec_flatten(inArr)
	for arr in inArr:
		if not isinstance(arr,basestring):
			outArr.append(rec_dflatten(arr,maxdepth-1))
		else:
			outArr.append(rec_flatten(arr))
	return outArr

#flatten
def rec_flatten(inArr):
	outArr = []
	if len(inArr)>0 and not isinstance(inArr,basestring):
		for arr in inArr:
			if isinstance(arr, basestring):
				outArr.append(arr)
			else:
				outArr.extend(rec_flatten(arr))
		return outArr
	return inArr
#END OBSOLETE-----------------------------------
	
#def main():
#	fpath=getargs()[0]
#	print fgrp()
	

#if __name__=="__main__":
#	main()

