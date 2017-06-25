#!/usr/bin/env python
import numpy as np
import scipy as sp
import os, time, glob
import matplotlib.pyplot as plt
import pdb

def rearrange_dost(dost_coefficients):
	N = len(dost_coefficients)
	_,freq_band_widths = band_width_partitioning(N)
	time_band_widths = N/freq_band_widths
	matrix_dost_coefficients = np.zeros((N,N))
	kk,ii,counter = 0,0,0
	matrix_of_indices = np.zeros((N,N))
	for hh in freq_band_widths:
		for jj in range(N,0,-time_band_widths[kk]):
			row_index = N-ii-hh
			column_index = abs(N-jj)
			#size1,size2 = matrix_of_indices[row_index:row_index+hh,column_index:column_index+time_band_widths[kk]].shape
			#matrix_of_indices[row_index:row_index+hh,column_index:column_index+time_band_widths[kk]]
			size1,size2 = hh, time_band_widths[kk]
			matrix_of_indices[row_index:row_index+hh,column_index:column_index+time_band_widths[kk]] = counter*np.ones((size1,size2))
			#pdb.set_trace()
			counter += 1
		kk += 1
		ii += hh
	matrix_dost_coefficients = dost_coefficients[matrix_of_indices.astype(int)]

	return matrix_dost_coefficients



def rearrange_dost_back(matrix_dost_coefficients):
	N = matrix_dost_coefficients.shape[0]
	_,freq_band_widths = band_width_partitioning(N)
	time_band_widths = N/freq_band_widths
	dost_coefficients = np.zeros(N)
	ii,kk,counter = 0,0,0
	for hh in freq_band_widths:
		for jj in range(N,0,-time_band_widths[kk]):
			row_index = N-ii-hh
			column_index = abs(N-jj)
			dost_coefficients[counter] = np.mean(matrix_dost_coefficients[row_index:row_index+hh,column_index:column_index+time_band_widths[kk]])
			counter +=1
		kk += 1
		ii += hh

	return dost_coefficients





def band_width_partitioning(samples):
	max_band_width_exponent = int(np.log2(samples)-2)
	#pdb.set_trace()
	positive_band_width_exponents = np.array(range(max_band_width_exponent+1))
	band_width_powers = np.concatenate(([0],positive_band_width_exponents[::-1],[0],positive_band_width_exponents))
	band_width2 = np.power(2,band_width_powers)
	number_of_freq_bands = len(band_width2)
	return number_of_freq_bands, band_width2

x = np.array([1,4,2,3,5,2,3,1])
print(x)
print(rearrange_dost(x))
print(rearrange_dost_back(rearrange_dost(x)))
