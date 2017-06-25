#!/usr/bin/env python
import numpy as np
import scipy as sp
import os, time, glob
import matplotlib.pyplot as plt
import pdb

def idost(dost_coefficients):
	rows_dost_coefficients = dost_coefficients.shape[0]
	reconstructed_signal = np.zeros(dost_coefficients.shape)
	tmp_signal = np.zeros(dost_coefficients.shape,dtype = complex)
	number_of_freq_bands, band_widths = band_width_partitioning(rows_dost_coefficients)
	
	counter = 0
	for ll in range(number_of_freq_bands):
		frequency_width_cutter = band_widths[ll];
		cut_x = dost_coefficients[counter:counter+frequency_width_cutter,]
		#pdb.set_trace()
		if frequency_width_cutter ==1 :
			tmp_signal[counter:counter+frequency_width_cutter,] = cut_x;
		else:
			tmp_signal[counter:counter+frequency_width_cutter,] = Fourier(cut_x);

		counter += frequency_width_cutter

	reconstructed_signal = np.real(iFourier(tmp_signal))
	return reconstructed_signal

def Fourier(x):
	y = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(x)))/np.sqrt(len(x))
	return y

def iFourier(x):
	y = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(x)))*np.sqrt(len(x))
	return y


def band_width_partitioning(samples):
	max_band_width_exponent = int(np.log2(samples)-2)
	positive_band_width_exponents = np.array(range(max_band_width_exponent+1))
	band_width_powers = np.concatenate(([0],positive_band_width_exponents[::-1],[0],positive_band_width_exponents))
	band_width2 = np.power(2,band_width_powers)
	number_of_freq_bands = len(band_width2)
	return number_of_freq_bands, band_width2


#x = np.array([88.9073,25.0013,38.2623,98.1975,83.0257,84.7297,9.5494,38.7562,87.1301,61.9834,11.3675,97.8701,9.5176,20.4848,18.5307,63.2538])
#print(idost(x))
