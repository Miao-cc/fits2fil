import sys
import pywt
import math
import ephem
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
from itertools import groupby
from scipy import interpolate
from itertools import chain

from presto import filterbank



secperday = 3600 * 24
name = sys.argv[1]

def wavelet(sig,threshold = 3, level=3, wavelet='db8'):
    sigma = sig.std()
    dwtmatr = pywt.wavedec(data=sig, wavelet=wavelet, level=level)
    denoised = dwtmatr[:]
    denoised[1:] = [pywt.threshold(i, value=threshold*sigma, mode='soft') for i in dwtmatr[1:]]
    smoothed_sig = pywt.waverec(denoised, wavelet, mode='sp1')[:sig.size]
    noises = sig - smoothed_sig
    return smoothed_sig,noises

def distance_c(bandpass,min_int):
    res=np.diff(bandpass)
    res=np.insert(res,0,1)
    breaks=np.array(bandpass)[np.array(res)>min_int]
    bandpass=list(bandpass)
    [bandpass.insert(list(bandpass).index(i),'div') for i in breaks]
    result=[list(g) for k,g in groupby(bandpass,lambda x:x=='div') if not k]
    return result
def find_wide(cluster):
    cluster_list = np.arange(len(cluster))
    cluster_length_list = [len(cluster[i]) for i in cluster_list]
    wide_band_index =np.where(np.array(cluster_length_list)>25)[0]
    for k in wide_band_index:
        print(k)
        cluster[k]= np.arange(min(cluster[k]),max(cluster[k])+1)
        print(cluster[k])
    return cluster   
def find_minnimun_good(i,bandpass,idxgood):
    min_good = (idxgood - i).sort()
    min_good_thre = min_good[9]
    near_good = idxgood[idxgood - i <= min_good_thre]
    bandpass[i] = bandpass[i]/(bandpass[i]/bandpass[near_good].sum())


fo = filterbank.FilterbankFile(name)

check_len = 5*60
nsample = fo.nspec
nchan = fo.nchan
tsample = fo.tsamp

nsblk = 1024
nsubint = 512
num = 5

subint_samp = nsblk*nsubint*5
print(nsample, nchan, tsample, subint_samp, nsubint, nsample/subint_samp)

loop = 0
checkNum = 0
bandpass = np.zeros(nchan)
for num in range(int(nsample/subint_samp)):
    data = fo.get_spectra(num*subint_samp,subint_samp).data
    bandpass_part = data.mean(axis=1)
    print(num, loop, loop*subint_samp*tsample, data.shape, bandpass_part.shape)
    if loop*subint_samp*tsample<check_len:
        loop = loop+1
        bandpass = bandpass_part + bandpass
    else:
        bandpass = bandpass/(loop+1)
        sig, nos=wavelet(bandpass,level=4)
        idxarr = np.arange(len(nos))
        thremin = np.mean(nos)-0.3*np.std(nos) ##threshold for RFI detecting
        thremax = np.mean(nos)+0.3*np.std(nos)
        idxbad = idxarr[(nos<thremin) | (nos>thremax)]
        cluster = distance_c(idxbad,20) ##threshold for RFI detecting
        #print(cluster)
        idxbad_cluster = find_wide(cluster)
        #print(idxbad_cluster)
        idxbad_wideband = list(chain(*idxbad_cluster))
        #print(idxbad_wideband)
        idxgood = list(set(idxarr)-set(idxbad_wideband))
        np.savetxt('%s_%s-rfi-channel-list.txt' %(name, str(num).rjust(3,'0')),idxbad_wideband, fmt='%i')
        print('RFI Ratio:%.3f%%'%(len(idxbad_wideband)*100.0/len(bandpass)))
        loop = 0
        checkNum = checkNum +1
        bandpass = np.zeros(nchan)
fo.close()

"""
c=min(idxgood)
d=max(idxgood)
xnew=np.arange(c,d,1)

f=interpolate.pchip_interpolate(idxgood,bandpass[idxgood],xnew,axis=0)
ahead=[]
tail=[]
if c>0:
   ahead=list(np.zeros(c-0)+bandpass[c])
if d<len(bandpass):
   tail=list(np.zeros(len(bandpass)-d)+bandpass[d])
f=list(f)
f1=ahead+f
f2=f1+tail
f2=np.array(f2)
####new data###
idxbad_wideband=np.array(idxbad_wideband)
data=data.astype(float)
data[:,idxbad_wideband]=f2[idxbad_wideband]*1.0
plt.plot(bandpass,color='black',label='bandpass')
plt.plot(f2)
plt.plot(data.sum(axis=0)*1.0/l,color='red',label='bandpass with rfi zap')
plt.scatter(idxbad_wideband,bandpass[idxbad_wideband],label='RFI channel')
plt.legend()
plt.show()
"""
