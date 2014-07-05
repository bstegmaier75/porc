#!/usr/bin/python -OO
#
# Python Open Room Correction (PORC)
# Copyright (c) 2012 Mason A. Green
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
#   More details about the parallel filter can be found in the papers
#
#   Balazs Bank, "Perceptually Motivated Audio Equalization Using Fixed-Pole Parallel
#   Second-Order Filters", IEEE Signal Processing Letters, 2008.
#   http://www.acoustics.hut.fi/go/spl08-parfilt
#
#   Balazs Bank, "Direct Design of Parallel Second-order Filters for
#   Instrument Body Modeling", International Computer Music Conference,
#   Copenhagen, Denmark, Aug. 2007.
#   http://www.acoustics.hut.fi/go/icmc07-parfilt
#
#   For Mixed-Phase Compensation, see:
#   "Mixed Time-Frequency approach for Mulitpoint Room Response Equalization," by
#   Alberto Carini, et al.

# Python libs
import sys
import textwrap
import wave
from contextlib import closing
import struct

# Scipy, Numpy, and matplotlibs
import numpy as np
import scipy as sp
import scipy.io as sio
import scipy.signal as sig
from scipy.fftpack import ifft, fft, rfft, irfft
from scipy.interpolate import pchip
from scipy.io import wavfile
from scipy.signal import convolve as conv
from scipy.stats import kurtosis, nanstd, linregress
from scipy.stats import norm as Gaussian
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os.path as path
from os.path import isfile, isdir, abspath, commonprefix

# PORC source files
from parfiltid import parfiltid
from tfplot import tfplot, tfplots, debug_log_plot
from freqpoles import freqpoles
from fgroup import fgrp

# Ignore warnings
import warnings; warnings.filterwarnings('ignore')

# MiniDSP's OpenDRC box likes 6144 taps

def rceps(x): 
	y = sp.real(ifft(sp.log(sp.absolute(fft(x)))))
	n = len(x) 
	if (n%2) == 1:
		ym = np.hstack((y[0], 2*y[1:n/2], np.zeros(n/2-1)))
	else:
		ym = np.hstack((y[0], 2*y[1:n/2], y[n/2+1], np.zeros(n/2-1)))
	ym = sp.real(ifft(sp.exp(fft(ym)))) 
	return (y, ym)
	
def parfilt(Bm, Am, FIR, x):
	y = np.zeros(x.size)
	for k in range(Am.shape[1]):
		y += np.ravel(sig.lfilter(Bm[:,k], Am[:,k], x))
	y += np.ravel(sig.lfilter(np.hstack([FIR]), np.hstack([1]), x))
	return y
	
# Normalize signal
def norm(y): return y/np.fabs(y).max()

def dB2Mag(dB):
	return 10**((dB)/20.)
	
def Mag2dB(mag):
	return 20.*np.log10(mag)
 
# The Median Absolute Deviation along given axis of an array
# From statsmodels lib
def mad(a, c=Gaussian.ppf(3/4.), axis=0):  # c \approx .6745
	a = np.asarray(a)
	return np.median((np.fabs(a-np.median(a,0)))/c, axis=axis)
	
def mad_sdev(data):
	return mad(data)/np.std(data)

	
#check if all values in an array are the same
def identical(arr):
	return arr[1:] == arr[:-1]

#return the next odd value less than or equal to val
def makeodd(val):
	if val%2==0:
		return val-1
	else:
		return val
		
def nrange(stop, num, start=0):
	if start>stop:
		tmp=start
		start=stop
		stop=tmp
	return np.unique(np.rint(np.linspace(start,stop,num)))
	

	

	
	
#round down
def rdn(val,step):
	return int(float(val)/step)*step
	
#calculate RT60 and Schroeder frequency based on data
def schroeder(data, fs, dbattl=-20, dbattu=-5, samples=20, rvol = 50, rounddown=1, debug=False):
	rdata2 = data[::-1]**2
	linsig = logtolin(norm(cumtrapz(rdata2)[::-1]))#schroeder integral of the impulse response
	if debug:
		plt.plot(linsig)
		plt.show()
		
	varr=np.linspace(dbattu, dbattl, samples) #attenuation values (in db)
	tarr=[-1]*samples #times corresponding to values in varr

	#calculate decay times
	for i,v in enumerate(linsig):
		for idx in range(samples):
			if v<varr[idx] and tarr[idx]<0:
				tarr[idx]=i/float(fs)
	
	#perform linear regression analysis on the data
	slope,ic,r,p,err = linregress(varr,tarr)
	
	if debug:
		print "y:",tarr
		print "x:",varr
	
	#estimate decay time to -60db
	rt60 = slope*-60+ic
	#schroeder frequency for a room of volume rvol (m^3)
	schroeder=2000.*np.sqrt(rt60/rvol)
	
	if debug:
		xx=np.linspace(-5,-60)
		yy=slope*xx+ic
		plt.plot(varr,tarr,"go",xx,yy,"r-")
		plt.show()
		print "rt60=",rt60
		print "schroeder=",schroeder
		print "rounded down to nearest 50Hz: ",rdn(schroeder,50)
	
	return rdn(schroeder, rounddown)
	
#Davis critical frequency (Don Davis, Eugene Patronis: Sound System Engineering)
def daviscf(minroomdim):
	return (3*344)/minroomdim
	
def peakoffset(data):
	max=0
	offset=0
	for idx, val in enumerate(data):
		c = np.fabs(val)
		if c>max:
			max=c
			offset=idx
	return offset
	#return 48000
	
# align impulse data peaks
def peakalign(dataArr):
	offsets = np.zeros(len(dataArr))
	maxdatalen = 0
	for idx, data in enumerate(dataArr):
		offsets[idx] = peakoffset(data)
		if len(data)>maxdatalen:
			maxdatalen=len(data)
	lpad = np.max(offsets)-offsets
	print lpad
	print offsets
	maxdatalen+=max(lpad)
	print len(dataArr)
	pdataArr = np.zeros((len(dataArr),maxdatalen))
	print pdataArr
	for idx,data in enumerate(dataArr):
		print "aa",maxdatalen-(len(data)+lpad[idx])
		lpadz = np.zeros(lpad[idx])
		rpadz = np.zeros(maxdatalen-(len(data)+lpad[idx]))
		pdataArr[idx]=np.append(lpadz,np.append(data,rpadz))
	return pdataArr
		
# root mean square deviation
def rmsd(data):
	return np.sqrt(np.mean(np.power(data,2)))
	
#recursive trimming function based on the kurtosis of the data
# - The function attempts to trim the impulse so that the kurtosis of the first bin is over the threshold)
# - Bin width should be set to the hop size used in mixed phase correction calculation
# - This method may leave a small number of zero or near-zero samples at the beginning of the impulse, while
#   still ensuring the non-negativity of the kurtosis of the first bin of width binw
def ktrim(data, threshold, maxdepth, binw):
	data = norm(np.real(data))
	bins=np.int(len(data)/binw)
	firstpositive=-1
	globalmax=[-1.,-1] #max. kurtosis, bin no.
	kvals = np.zeros(bins)
	for b in range(bins):
		k = kurtosis(data[(b*binw):((b+1)*binw)],None,True,False)
		if k>threshold:
			if firstpositive<0:
				firstpositive=b
		if k>globalmax[0]: #save bin with maximum kurtosis
			globalmax[0]=k 
			globalmax[1]=b
			kvals[b]=k
			
	if firstpositive==0:
		print "First bin is already positive, no trimming necessary."
		return data
	offset=0
	if globalmax[1]>0:
		print "Locating global maximum"
		maxb = globalmax[1]
		while maxb>0 and kvals[maxb-1]<kvals[maxb]: #
			maxb-=1
		print "Starting recursive search at sample ", binw*maxb
		offset=maxb*binw
		
	return ktrimrec(data[offset:],threshold,maxdepth,binw,binw,binw)
	

def ktrimrec(data,threshold,maxdepth, binw, offset, stepsize):
	if maxdepth>0:
		data = norm(np.real(data))
		bins=np.int(len(data)/binw)
		firstpositive=-1
		for b in range(bins):
			k = kurtosis(data[(b*binw)+offset:((b+1)*binw)+offset],None,True,False)
			if k>threshold:
				if firstpositive<0:
					firstpositive=b

		if firstpositive==0:
			print "-: ", offset
			return ktrimrec(data,threshold,maxdepth-1,binw,offset-stepsize/2, stepsize/2) #start converging on the exact solution
		elif firstpositive==1:
			print "+: ", offset
			return ktrimrec(data,threshold,maxdepth-1,binw,offset+stepsize/2, stepsize/2) #start converging on the exact solution
		else:
			print "Something is not right.. Returning data as-is"
			return data;
	else:
		print "Max depth reached, returning at offset ", offset
		return data[offset:]
		

# smarter trimming algorithm
# works backwards from the impulse
# binw - bin width, smaller values increase accuracy but might cause false positives
# step - iteration step size, in samples. larger values are faster but less accurate
# threshold - 0.01 seems to work well
def smarttrim(data, binw, threshold, step, Fs=48000., debug=False):
	sp = -1
	gmaxv = -1 #absolute value of global maximum, we assume this is a part of the actual impulse
	for pos, sample in enumerate(data):
		if abs(sample)>=gmaxv:
			sp = pos #start search at global max
			gmaxv = sample		
	sp+=np.int(binw/2)
	idv = rmsd(data[sp-binw:sp])
	cdv = 1
	#print "InitSp=",sp
	while np.fabs(cdv)>np.fabs(threshold*idv) and (sp-step)>0:
		sp-=step
		cdv=rmsd(data[sp-binw:sp])
	if debug:
		print "Trimming ", sp/Fs
	return data[sp:]

def minphase(data):
	cp, mp = rceps(data)
	return mp;

def cnorm(c):
	maxmag = max(np.abs(c))
	return c/maxmag
	
	
#linear to db
def lintolog(x):
	return 10**(x/10.)
#db to linear
def logtolin(x):
	return 10.*np.log10(x)
		
def roomcomp(impresps, filter, target, ntaps, mixed_phase, opformat, trim, trimthreshold, noplot, strim, tfreq, lfpoles, hfpoles, debug):
	data = []
	Fs = 0
	if len(impresps)>1:
		if debug:
			print "Averaging ", len(impresps), " impulses"
		fsArr = [None] * len(impresps)
		dataArr = [None] * len(impresps)
		#print len(dataArr)
		for index, impulse in enumerate(impresps):
			fsArr[index], dataArr[index] = wavfile.read(impulse)
		if not identical(fsArr):
			print 'Sample rate mismatch'
			return
		Fs = fsArr[0]
		###########################
		# Average impulse responses
		###########################
		maxlen = -1

		for impulse in dataArr:
			maxlen=max(maxlen,len(impulse))
		ftimpulses = [None] * len(dataArr)
		for idx,impulse in enumerate(dataArr):
			ftimpulses[idx] = logtolin(fft(dataArr[idx],maxlen))
		ftavg = np.mean(ftimpulses,0)
		data=norm(minphase(ifft(lintolog(ftavg))))
	else:
		print "Loading impulse response"
		# Read impulse response
		Fs, data = wavfile.read(impresps[0])
		data = norm(np.hstack(data))

	if trim:
		print "Removing leading silence"
		for spos,sval in enumerate(data):
			if abs(sval)>trimthreshold:
				lzs=max(spos-1,0)
				ld =len(data)
				print 'Impulse starts at position ', spos, '/', len(data)
				print 'Trimming ', float(lzs)/float(Fs), ' seconds of silence'
				data=data[lzs:len(data)] #remove everything before sample at spos
				break
	elif strim:
		data=smarttrim(data,Fs*0.024,trimthreshold,1)
	
	if(debug):
		print "Schroeder freq (assuming 30m^3 room): ",schroeder(data,Fs,rvol=30, debug=debug)  
		print "\nSample rate = ", Fs
  
	print "\nGenerating correction filter"

	###
	## Logarithmic pole positioning
	###

	fplog = np.hstack((sp.logspace(sp.log10(20.), sp.log10(float(tfreq-25)), float(lfpoles)), sp.logspace(sp.log10(float(tfreq+25)), 
			sp.log10(20000.), float(hfpoles))))
	plog = freqpoles(fplog, Fs)

	###
	## Preparing data
	###

	# making the measured response minumum-phase
	cp, minresp = rceps(data)

	# Impulse response
	imp = np.zeros(len(data), dtype=np.float64)
	imp[0]=1.0

	# Target
	outf = []
	db = []

	if target is 'flat':
	
		# Make the target output a bandpass filter
		Bf, Af = sig.butter(4, 30/(Fs/2), 'high')
		outf = sig.lfilter(Bf, Af, imp) 
	
	else:
	
		# load target file
		t = np.loadtxt(target)
		frq = t[:,0]; pwr = t[:,1]
		
		if not frq[0]==0:
			frq = np.append(np.zeros(1),frq)
			pwr = np.append(pwr[0],pwr)
			
		if debug:
			print "Target params:"
			print frq
			print pwr
		
		# calculate the FIR filter via windowing method
		
		#anlin: added antisymmetric flag to ensure type 1 filter,
		#       raised numtaps to odd number closest but smaller than len(minresp)
		fir = sig.firwin2(makeodd(len(minresp)), frq, np.power(10, pwr/20.0), nyq = frq[-1], antisymmetric = False)	
		# Minimum phase, zero padding	
		cp, outf = rceps(np.append(fir, np.zeros(len(minresp) - len(fir))))
		
			
	###
	## Filter design
	###

	#Parallel filter design
	(Bm, Am, FIR) = parfiltid(minresp, outf, plog)

	# equalized loudspeaker response - filtering the 
	# measured transfer function by the parallel filter
	equalizedresp = parfilt(Bm, Am, FIR, data)

	# Equalizer impulse response - filtering a unit pulse
	equalizer = norm(parfilt(Bm, Am, FIR, imp))
	
	# Pad equalizer with zeroes to length ntaps
	equalizer = np.append(equalizer, np.zeros(max(0,ntaps-len(equalizer))))

	# Windowing with a half hanning window in time domain
	han = np.hanning(ntaps*2)[-ntaps:]
	equalizer = han * equalizer[:ntaps]

	###
	## Mixed-phase compensation
	## Based on the paper "Mixed Time-Frequency approach for Multipoint
	## Room Rosponse Equalization," by A. Carini et al.
	## To use this feature, your Room Impulse Response should have all
	## the leading zeros removed.
	###
	if mixed_phase is True:
	
		# prototype function
		hp = norm(np.real(equalizedresp))

		# time integration of the human ear is ~24ms
		# See "Measuring the mixing time in auditoria," by Defrance & Polack
		hop_size = 0.024
		samples = hop_size * Fs

		bins = np.int(np.ceil(len(hp) / samples))
		
		tmix = 0

		# Kurtosis method
		for b in range(bins):
			
			start = np.int(b * samples)
			end = np.int((b+1) * samples)
			k = kurtosis(hp[start:end], None, True, False) #anlin: disabled scipy's default bias correction
			if k <= 0:
				tmix = b * hop_size
				break
		# truncate the prototype function
		taps = np.int(tmix*Fs)
		
		print "\nmixing time(secs) = ", tmix, "; taps = ", taps
	
		if taps > 0:
			# Time reverse the array
			h = hp[:taps][::-1]
			# create all pass filter
			phase = np.unwrap(np.angle(h))
			H = np.exp(1j*phase)
			# convert from db to linear
			mixed = np.power(10, np.real(H)/20.0)
			# create filter's impulse response
			mixed = np.real(ifft(mixed))
			
			# convolve and window to desired length
			equalizer = conv(equalizer, mixed)
			equalizer = han * equalizer[:ntaps]

			
			#data = han * data[:ntaps]
			#eqresp = np.real(conv(equalizer, data))
		else:
			print "zero taps; skipping mixed-phase computation"
			
		
	if opformat in ('wav', 'wav24'):
		# Write data
		wavwrite_24(filter, Fs, norm(np.real(equalizer)))
		print '\nOutput format is wav24'
		print 'Output filter length =', len(equalizer), 'taps'
		print 'Output filter written to ' + filter
		print "\nUse sox to convert output .wav to raw 32 bit IEEE floating point if necessary,"
		print "or to merge left and right channels into a stereo .wav"
		print "\nExample: sox leq48.wav -t f32 leq48.bin"
		print "         sox -M le148.wav req48.wav output.wav\n"

	elif opformat == 'wav32':
		wavwrite_32(filter, Fs, norm(np.real(equalizer)))
		print '\nOutput format is wav32'
		print 'Output filter length =', len(equalizer), 'taps'
		print 'Output filter written to ' + filter
		print "\nUse sox to convert output .wav to raw 32 bit IEEE floating point if necessary,"
		print "or to merge left and right channels into a stereo .wav"
		print "\nExample: sox leq48.wav -t f32 leq48.bin"
		print "         sox -M le148.wav req48.wav output.wav\n"
	elif opformat == 'bin':
		# direct output to bin avoids float64->pcm16->float32 conversion by going direct 
		#float64->float32
		f = open(filter, 'wb')
		norm(np.real(equalizer)).astype('float32').tofile(f)
		f.close()
		print '\nOutput filter length =', len(equalizer), 'taps'
		print 'Output filter written to ' + filter
	elif opformat == 'txt':
		txtwrite_f32(filter,norm(np.real(equalizer)))
	else:
		print 'Output format not recognized, no file generated.'


	###
	## Plots
	###
	if not noplot:
		data *= 500
		# original loudspeaker-room response
		tfplot(data, Fs, avg = 'abs')
		# 1/3 Octave smoothed
		tfplots(data, Fs, 'r')

		#tfplot(mixed, Fs, 'r')

		# equalizer transfer function
		tfplot(0.75*equalizer, Fs, 'g')
		# indicating pole frequencies
		plt.vlines(fplog, -2, 2, color='k', linestyles='solid')

		# equalized loudspeaker-room response
		tfplot(equalizedresp*0.01, Fs, avg = 'abs')
		# 1/3 Octave smoothed
		tfplots(equalizedresp*0.01, Fs, 'r')

		# Add labels
		# May need to reposition these based on input data
		plt.text(325,30,'Unequalized loudspeaker-room response')
		plt.text(100,-15,'Equalizer transfer function')
		plt.text(100,-21,'(Black lines: pole locations)')
		plt.text(130,-70,'Equalized loudspeaker-room response')

		a = plt.gca()
		a.set_xlim([20, 20000])
		a.set_ylim([-80, 80])
		plt.ylabel('Amplitude (dB)', color='b')
		plt.xlabel('Frequency (Hz)')
		plt.grid()
		plt.legend()
		plt.show()

def wavwrite_24(fname, fs, data):
	data_as_bytes = (struct.pack('<i', int(samp*(2**23-1))) for samp in data)
	with closing(wave.open(fname, 'wb')) as wavwriter:
		wavwriter.setnchannels(1)
		wavwriter.setsampwidth(3)
		wavwriter.setframerate(fs)
		for data_bytes in data_as_bytes:
			wavwriter.writeframes(data_bytes[0:3])
			
def wavwrite_32(fname, fs, data):
	data_as_bytes = (struct.pack('<i', int(samp*(2**31-1))) for samp in data)
	with closing(wave.open(fname, 'wb')) as wavwriter:
		wavwriter.setnchannels(1)
		wavwriter.setsampwidth(4)
		wavwriter.setframerate(fs)
		for data_bytes in data_as_bytes:
			wavwriter.writeframes(data_bytes[0:4])
			
def txtwrite_f32(fname,data):
	data_f32 = data.astype('string')
	for index, strv in enumerate(data_f32):
		data_f32[index] = 'b' + str(index) + ' = ' + strv
	fstr = ',\r\n'.join(data_f32)
	file = open(fname,'wb')
	file.write(fstr)
	file.close()
		
	
	
			
def main():
	
	print

	mtxt = textwrap.dedent('''\
	Python Open Room Correction (PORC), version 0.1
	Copyright (c) 2012 Mason A. Green
	Based on the work of Dr. Balazs Bank
	''')

	bye = textwrap.dedent('''
	Example:
	./porc -t b&k.txt -n 8000 l48.wav leq48.bin
		
	See the README for detailed instructions
	''')

	import argparse
	from argparse import RawTextHelpFormatter

	parser = argparse.ArgumentParser(description = mtxt, epilog=bye, formatter_class=RawTextHelpFormatter)

	# Positionals
	parser.add_argument('impresp', metavar='I', type=str, nargs="+", help='Mesaured impulse response. Can be a file, a list of files or a folder')
	parser.add_argument('filter', metavar='F', type=str, help='Output filter path and filename (or common prefix) without extension')

	# Options
	parser.add_argument("-t", dest="target", default='flat',
						help="target curve", metavar="FILE")
	parser.add_argument("-n", dest="ntaps", default = 6144,
						help="filter length, in taps. Default = len(input)", type=int) 
	parser.add_argument('--mixed', action='store_true', default = False,
						help="Implement mixed-phase compensation. see README for details") 
	parser.add_argument("-o", dest="opformat", default = 'bin',
						help="Output file type, default bin optional wav, txt", type=str) 
	parser.add_argument("--lfpoles", dest="lfpoles", default = 14,
						help="Number of low frequency band poles (20Hz->tfreq). Default=14", type=int) 
	parser.add_argument("--hfpoles", dest="hfpoles", default = 13,
						help="Number of high frequency band poles (tfreq->20kHz). Default=13, zero to disable HF correction", type=int) 
	parser.add_argument("--tfreq", dest="tfreq", default = 200,
						help="LF to HF band transition freq.", type=int) 
	parser.add_argument("-s", dest="trimthreshold", default = 0.05,
						help="Normalized silence threshold. Default = 0.05", type=float) 
	parser.add_argument('--trim', action='store_true', default = False,
						help="Trim leading silence")
	parser.add_argument('--strim', action='store_true', default = False,
						help="Smart trim (experimental)")
	parser.add_argument('--noplot', action='store_true', default = False,
						help="Do not open a plot window")
	parser.add_argument('--debug', action='store_true', default = False,
						help="Print debug information")

	args = parser.parse_args()
	
	
	

	
	
	ext = args.opformat[:3] # this won't work with extensions over 3 characters, the if-elif block below is for redundancy
	print ext
	if args.opformat=="wav" or args.opformat=="wav32" or args.opformat=="wav24":
		ext="wav"
	elif args.opformat=="bin":
		ext="bin"
	
	if len(args.impresp)>1 and not isinstance(args.impresp,basestring):
		for impulse in args.impresp:
			if not isfile(impulse):
				print "Invalid input"
				print "Valid inputs are folders, individual files and space separated lists of files"
	elif isdir(args.impresp[0]):
		lfile,lname = fgrp(abspath(args.impresp[0]))
		for idx in range(len(lfile)):
			print lfile[idx]
			roomcomp(lfile[idx], args.filter+lname[idx]+"."+ext, args.target, args.ntaps, args.mixed, args.opformat, args.trim, args.trimthreshold, args.noplot, args.strim, args.tfreq, args.lfpoles, args.hfpoles, args.debug)
		return
	

	roomcomp(args.impresp, args.filter+"."+ext, args.target, args.ntaps, args.mixed, args.opformat, args.trim, args.trimthreshold, args.noplot, args.strim, args.tfreq, args.lfpoles, args.hfpoles, args.debug)

if __name__=="__main__":
	main()  
