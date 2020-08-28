#http://nbviewer.jupyter.org/github/ridlo/exoplanet_notebook/blob/master/bls_test02.ipynb

import matplotlib.pyplot as plt
import numpy as np
import math

import astropy.units as u
from astropy.constants import G, R_sun, M_sun, R_jup, M_jup, R_earth, M_earth

#%matplotlib inline


#import numpy
#from astropy import units as u
import scipy
from astropy.coordinates import SkyCoord
import glob
import sys
import pylab
#import math
import os
import platform
#import matplotlib.pyplot as plt


def bls(t, x, qmi, qma, fmin, df, nf, nb):
    """Frist trial, BLS algorithm, only minor modification from author's code"""
    
    n = len(t); rn = len(x)
    #! use try
    if n != rn:
        print ("Different size of array, t and x")
        return 0

    rn = float(rn) # float of n

    minbin = 5
    nbmax = 2000
    if nb > nbmax:
        print ("Error: NB > NBMAX!")
        return 0

    tot = t[-1] - t[0] # total time span

    if fmin < 1.0/tot:
        print ("Error: fmin < 1/T")
        return 0

    # parameters in binning (after folding)
    kmi = int(qmi*nb) # nb is number of bin -> a single period
    if kmi < 1: 
        kmi = 1
    kma = int(qma*nb) + 1
    kkmi = rn*qmi # to check the bin size
    if kkmi < minbin: 
        kkmi = minbin

    # For the extension of arrays (edge effect: transit happen at the edge of data set)
    nb1 = nb + 1
    nbkma = nb + kma
        
    # Data centering
    t1 = t[0]
    u = t - t1
    s = np.median(x) # ! Modified
    v = x - s

    bpow = 0.0
    p = np.zeros(nf)
    # setup array for power vs period plot
    powerPeriod=[]
    # Start period search
    for jf in range(nf):
        #f0 = fmin + df*jf # iteration in frequency not period
        #p0 = 1.0/f0

        # Actually iterate in period
        p0 = startPeriod + dp*jf
        f0 = 1.0/p0

        # Compute folded time series with p0 period
        ibi = np.zeros(nbkma)
        y = np.zeros(nbkma)
        # Median version
        yMedian = np.zeros(shape=(nf,n))
        yMedian.fill(np.nan)
        for i in range(n):
            ph = u[i]*f0 # instead of t mod P, he use t*f then calculate the phase (less computation)
            ph = ph - int(ph)
            j = int(nb*ph) # data to a bin 
            ibi[j] = ibi[j] + 1 # number of data in a bin
            y[j] = y[j] + v[i] # sum of light in a bin
            yMedian[j][i]=v[i]

        #print (ibi)
        #print (y)
        #print (yMedian)
        #print (np.nanmedian(yMedian[2,:]))

        # Repopulate y[j] and ibi[j] with the median value
        for i in range(nb+1):
            #print (i)
            ibi[i]=1
            y[i]=np.nanmedian(yMedian[i,:])

        #print (ibi)
        #print (y)
        
        
        
        # Extend the arrays  ibi()  and  y() beyond nb by wrapping
        for j in range(nb1, nbkma):
            jnb = j - nb
            ibi[j] = ibi[jnb]
            y[j] = y[jnb]

        #print (ibi)
        #print (y)
        #sys.exit()

        # Compute BLS statictics for this trial period
        power = 0.0

        for i in range(nb): # shift the test period
            s = 0.0
            k = 0
            kk = 0
            nb2 = i + kma
            # change the size of test period (from kmi to kma)
            for j in range(i, nb2): 
                k = k + 1
                kk = kk + ibi[j]
                s = s + y[j]
                if k < kmi: continue # only calculate SR for test period > kmi
                if kk < kkmi: continue # 
                rn1 = float(kk)
                powo = s*s/(rn1*(rn - rn1))
                if powo > power: # save maximum SR in a test period
                    power = powo # SR value
                    jn1 = i # 
                    jn2 = j
                    rn3 = rn1
                    s3 = s

        power = math.sqrt(power)
        p[jf] = power
        powerPeriod.append([p0,power])
        

        if power > bpow:
            # If it isn't an resonance of a day
            if not ((p0 > 0.95 and p0 < 1.05) or (p0 > 1.95 and p0 < 2.05) or (p0 > 2.98 and p0 < 3.02) or (p0 > 6.65 and p0 < 6.67) or (p0 > 3.32 and p0 < 3.34) or (p0 > 3.64 and p0 < 3.68)):            
                bpow = power # Save the absolute maximum of SR
                in1 = jn1
                in2 = jn2
                qtran = rn3/rn
                # depth = -s3*rn/(rn3*(rn - rn3))
                # ! Modified
                high = -s3/(rn - rn3)
                low = s3/rn3
                depth = high - low
                bper = p0

    # print powerperiod

    plt.figure(figsize=(15,6))

    #print (powerPeriod)
    powerPeriod=np.asarray(powerPeriod)
    plt.subplot(1, 2, 1)
    plt.plot(powerPeriod[:,0], powerPeriod[:,1], 'r.')
    #fite = np.zeros(nbin) + res[8] # H
    #fite[res[1]:res[2]+1] = res[9] # L
        
    #plt.plot(phase, fite)
    #plt.gca().invert_yaxis()
    plt.title("EELBS Period Trials")
    plt.xlabel(r"Trialled Period")
    plt.ylabel(r"Likelihood")
##        plt.text(0.01,0.05,"Best SR: " +str(res[0]), alpha=0.7)
##        plt.text(0.01,0.07,"\nIngress: " + str(res[1]), alpha=0.7)
##        plt.text(0.01,0.09,"\nEgress: "+ str(res[2]), alpha=0.7)
##        plt.text(0.01,0.11,"\nq: "+ str(res[3]), alpha=0.7)
##        plt.text(0.01,0.13,"\nDepth: "+ str(-res[4]), alpha=0.7)
##        plt.text(0.01,0.15,"\nPeriod: "+ str(res[5]), alpha=0.7)
##        plt.text(0.01,0.17,"\nSDE: "+ str(res[6]), alpha=0.7)
    plt.savefig(os.path.join(eelbsPath,str(file).split("/")[-1].split("\\")[-1].replace(".csv","").replace("_calibExcel","")+'_EELBS_Plot.png'))
    #plt.savefig(os.path.join(eelbsPath,str(file).split("\")[-1]+'_EELBS.eps'))
    print (str(file).split("/")[-1].split("\\")[-1].replace(".csv","").replace("_calibExcel","")+'_EELBS_Plot.png')
    #plt.show()
    
    
    # ! add
    sde = (bpow - np.mean(p))/np.std(p) # signal detection efficiency

    return bpow, in1, in2, qtran, depth, bper, sde, p, high, low

# t is times, f is fluxes


##c
##c     Output parameters:
##c     ~~~~~~~~~~~~~~~~~~
##c
##c     p    = array {p(i)}, containing the values of the BLS spectrum
##c            at the i-th frequency value -- the frequency values are 
##c            computed as  f = fmin + (i-1)*df
##c     bper = period at the highest peak in the frequency spectrum
##c     bpow = value of {p(i)} at the highest peak
##c     depth= depth of the transit at   *bper*
##c     qtran= fractional transit length  [ T_transit/bper ]
##c     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
##c     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]
##c
##c
##c     Remarks:
##c     ~~~~~~~~ 
##c
##c     -- *fmin* MUST be greater than  *1/total time span* 
##c     -- *nb*   MUST be lower than  *nbmax* 
##c     -- Dimensions of arrays {y(i)} and {ibi(i)} MUST be greater than 
##c        or equal to  *nbmax*. 
##c     -- The lowest number of points allowed in a single bin is equal 
##c        to   MAX(minbin,qmi*N),  where   *qmi*  is the minimum transit 
##c        length/trial period,   *N*  is the total number of data points,  
##c        *minbin*  is the preset minimum number of the data points per 
##c        bin.

# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")
    trimPath = os.path.join(parentPath,"trimcats\\")
    eelbsPath = os.path.join(parentPath,"eelbs\\")

if "Mac" or "Darwin" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
    trimPath = os.path.join(parentPath,"trimcats//")
    eelbsPath = os.path.join(parentPath,"eelbs//")
    
#fileList=glob.glob('*.ps*')
#fileList = glob.glob(inputPath)
#print (fileList)
#create directory structure
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

if not os.path.exists(outcatPath):
    os.makedirs(outcatPath)

if not os.path.exists(checkPath):
    os.makedirs(checkPath)

if not os.path.exists(trimPath):
    os.makedirs(trimPath)

if not os.path.exists(eelbsPath):
    os.makedirs(eelbsPath)


fileList = glob.glob(outcatPath+'*diffExcel*csv')
r=0
#print (fileList)

##c     Input parameters:
##c     ~~~~~~~~~~~~~~~~~
##c
##c     n    = number of data points
##c     t    = array {t(i)}, containing the time values of the time series
##c     x    = array {x(i)}, containing the data values of the time series
##c     u    = temporal/work/dummy array, must be dimensioned in the 
##c            calling program in the same way as  {t(i)}
##c     v    = the same as  {u(i)}
##c     nf   = number of frequency points in which the spectrum is computed
##c     fmin = minimum frequency (MUST be > 0)
##c     df   = frequency step
##c     nb   = number of bins in the folded time series at any test period       
##c     qmi  = minimum fractional transit length to be tested
##c     qma  = maximum fractional transit length to be tested

qmi = 0.01
qma = 0.1
startPeriod=0.1
endPeriod=3.0
nf=1000


# calculate period range
fmin = 1/endPeriod
fmax = 1/startPeriod
#fmin = 0.333
df = (fmax-fmin)/nf
dp = (endPeriod-startPeriod)/nf
#print (df)
#nf = 100
nb = 200

for file in fileList:
    photFile = np.genfromtxt(file, dtype=float, delimiter=',')
    print ('**********************')
    print ('Testing: ' + str(file))
    #print (photFile[:,0])
    #print (photFile[:,1])
    t = photFile[:,0]
    f = photFile[:,1]
    res = bls(t, f, qmi, qma, fmin, df, nf, nb)
    if (res != 0): # If it did not fail, then do the rest.
        print ("Best SR: ", res[0], "\nIngress: ", res[1], "\nEgress: ", res[2], "\nq: ", res[3], \
    "\nDepth: ", res[4], "\nPeriod: ", res[5], "\nSDE: ", res[6])

        t1 = t[0]
        u = t - t1
        s = np.mean(f)
        v = f - s

        f0 = 1.0/res[5] #  freq = 1/T
        nbin = nb # number of bin
        n = len(t)
        ibi = np.zeros(nbin)
        y = np.zeros(nbin)
        phase = np.linspace(0.0, 1.0, nbin)

        for i in range(n):
            ph = u[i]*f0 
            ph = ph - int(ph)
            j = int(nbin*ph) # data to a bin 
            ibi[j] = ibi[j] + 1.0 # number of data in a bin
            y[j] = y[j] + v[i] # sum of light in a bin


        plt.subplot(1, 2, 2)
        plt.plot(phase, y/ibi, 'r.')
        fite = np.zeros(nbin) + res[8] # H
        fite[res[1]:res[2]+1] = res[9] # L
    
        plt.plot(phase, fite)
        plt.gca().invert_yaxis()
        plt.title("\nDepth: "+ str(-res[4]) + "     " + "Period: {0} d  bin: {1}".format(1/f0, nbin))
        plt.xlabel(r"Phase ($\phi$)")
        plt.ylabel(r"Mean value of $x(\phi)$ in a bin")
##        plt.text(0.01,0.05,"Best SR: " +str(res[0]), alpha=0.7)
##        plt.text(0.01,0.07,"\nIngress: " + str(res[1]), alpha=0.7)
##        plt.text(0.01,0.09,"\nEgress: "+ str(res[2]), alpha=0.7)
##        plt.text(0.01,0.11,"\nnq: "+ str(res[3]), alpha=0.7)
##        plt.text(0.01,0.13,"\nDepth: "+ str(-res[4]), alpha=0.7)
##        plt.text(0.01,0.15,"\nPeriod: "+ str(res[5]), alpha=0.7)
##        plt.text(0.01,0.17,"\nSDE: "+ str(res[6]), alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(eelbsPath,str(file).split("/")[-1].split("\\")[-1].replace(".csv","").replace("_calibExcel","")+'_EELBS_Plot.png'))
        #plt.savefig(os.path.join(eelbsPath,str(file).split("\")[-1]+'_EELBS.eps'))
        print (str(file).split("/")[-1].split("\\")[-1].replace(".csv","").replace("_calibExcel","")+'_EELBS_Plot.png')
        #plt.show()
        plt.clf()

        # Write text file
        texFileName=os.path.join(eelbsPath,str(file).split("/")[-1].split("\\")[-1].replace(".csv","").replace("_calibExcel","")+'_EELBS_Statistics.txt')
        print (texFileName)
        with open(texFileName, "w") as f:
            f.write("Best SR: " +str(res[0])+"\n")
            f.write("Ingress: " + str(res[1])+"\n")
            f.write("Egress: "+ str(res[2])+"\n")
            f.write("nq: "+ str(res[3])+"\n")
            f.write("Depth: "+ str(-res[4])+"\n")
            f.write("Period: "+ str(res[5])+"\n")
            f.write("SDE: "+ str(res[6])+"\n")


sys.exit()



