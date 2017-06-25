#!/usr/bin/env python
import numpy as np
import scipy as sp
import os, time, glob
import matplotlib.pyplot as plt
import pdb

def dost(time_series):
	rows_time_series = time_series.shape[0]
	dost_coefficients = np.zeros(time_series.shape,dtype = complex)
	# partition of the frequency space
	number_of_freq_bands, band_widths = band_width_partitioning(rows_time_series)
	Fx = Fourier(time_series)
	counter = 0
	for ll in range(number_of_freq_bands):
		frequency_width_cutter = band_widths[ll]
		cut_Fx = Fx[counter:counter+frequency_width_cutter,]
		#pdb.set_trace()
		if frequency_width_cutter == 1:
			dost_coefficients[counter:counter+frequency_width_cutter,] = cut_Fx;
		else:
			dost_coefficients[counter:counter+frequency_width_cutter,] = iFourier(cut_Fx);
		
		counter += frequency_width_cutter
	return dost_coefficients

def Fourier(x):
	y = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x)))/np.sqrt(len(x))
	return y

def iFourier(x):
	y = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(x)))*np.sqrt(len(x))
	return y

def band_width_partitioning(samples):
	max_band_width_exponent = int(np.log2(samples)-2)
	#pdb.set_trace()
	positive_band_width_exponents = np.array(range(max_band_width_exponent+1))
	band_width_powers = np.concatenate(([0],positive_band_width_exponents[::-1],[0],positive_band_width_exponents))
	band_width2 = np.power(2,band_width_powers)
	number_of_freq_bands = len(band_width2)
	return number_of_freq_bands, band_width2


#x = np.array([88.9073,25.0013,38.2623,98.1975,83.0257,84.7297,9.5494,38.7562,87.1301,61.9834,11.3675,97.8701,9.5176,20.4848,18.5307,63.2538]).reshape(4,4)
#print(dost(x))
