import numpy as np 
#import pyfits
import astropy.io.fits as pyfits
import os
import datetime
import time
import sys
from array import array
from pylab import *

##############################################################
# Adapted from Downsamp_FASTpsrfits_freq_time_4pol.py
# 
# Output data for a selected time & freq. range and frequency downsampling rate in Fits format.
# (Output 4 pols Fits format as input datafile)
#
# Miaocc 2019/10/16
##############################################################

filename=sys.argv[1]

starttime=datetime.datetime.now()


hdulist = pyfits.open(filename)
hdu0 = hdulist[0]
data0 = hdu0.data
header0 = hdu0.header
hdu1 = hdulist[1]
data1 = hdu1.data
header1 = hdu1.header

float_data=np.array(data1['DATA'])

print("+++++++++++++++++++++++++")

temp_float_dat_scl=np.array(data1['DAT_SCL'])
nline=header1['NAXIS2']
nsblk=header1['NSBLK']
tbin=header1['TBIN']
npol=header1['NPOL']
nsuboffs=header1['NSUBOFFS']
chan_bw=header1['CHAN_BW']
freq=header0['OBSFREQ']
nchan=header0['OBSNCHAN']
print('                Filename =', filename)
print('Input  Number of subints =', size(temp_float_dat_scl)/npol/nchan, nline)
print('       Num polarisations =', npol)
print('       Ntsamp each subint=', nsblk, size(float_data)/nline/npol/nchan, size(float_data)/nline/nchan/npol)
print('       Central freq(MHz) =', freq)
print('       Freq. bandwidth   =', header0['OBSBW'])
print('       Channel number    =', nchan)
print('       Channel width(MHz)=', chan_bw)
print('       data1[\'DATA\']     =', float_data.shape)
print('       Ntsamp each subint=', 'Same as input')
print('       Central freq(MHz) =', hdu0.header['OBSFREQ'])
print('       Freq. bandwidth   =', hdu0.header['OBSBW'])
print('       Channel number    =', hdu0.header['OBSNCHAN'])

pol_AA = float_data[:,:,0,:,0].reshape(-1,nchan).mean(axis=1)
pol_BB = float_data[:,:,1,:,0].reshape(-1,nchan).mean(axis=1)

outName_AA = filename.split('/')[-1] + "_AA.txt"
outName_BB = filename.split('/')[-1] + "_BB.txt"

np.savetxt(outName_AA, pol_AA)
np.savetxt(outName_BB, pol_BB)
