import numpy
from astropy import units as u
import scipy
from astropy.coordinates import SkyCoord
import glob
import sys
import pylab
import math
import os
import platform
import matplotlib.pyplot as plt

# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")
    #trimPath = os.path.join(parentPath,"trimcats\\")

if "Mac" or "Darwin" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
    #trimPath = os.path.join(parentPath,"trimcats//")
    
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

#if not os.path.exists(trimPath):
    #os.makedirs(trimPath)

# Load in list of used files
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)

fileList = glob.glob(outcatPath+'*diffExcel*csv')
r=0
#print (fileList)
for file in fileList:
    photFile = numpy.genfromtxt(file, dtype=float, delimiter=',')
    exists=os.path.isfile(file.replace('diff','calib'))
    if exists:
        calibFile = numpy.genfromtxt(file.replace('diff','calib'), dtype=float, delimiter=',')
        print ("Calibration difference")
        print (-(photFile[:,1]-calibFile[:,1])[0])
        calibDiff=-((photFile[:,1]-calibFile[:,1])[0])
    #print (photFile[:,1])
    #print (photFile[:,0])
    print (file)
    print (photFile[:,1])
    
    ax=plt.plot(photFile[:,0],photFile[:,1],'ro')
    #print (ax.lines)
    plt.show()
    del ax
    plt.clf()
    plt.cla()
    plt.close()
    plt.close("all")
    print ("Enter Dim (Upper) Cutoff Value:")
    brightV = input()
    print ("Enter Bright (Lower) Cutoff Value:")
    brightD = input()
    measureReject=[]
    for i in range(photFile.shape[0]):
        if photFile[i,1] < float(brightD) or photFile[i,1] > float(brightV) :
            measureReject.append(i)
            print (photFile[i,1])
            print ("REJECT")
    print (photFile.shape[0])
    photFile=numpy.delete(photFile, measureReject, axis=0)
    print (photFile.shape[0])

    # Output trimmed files
    numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"diffPeranso.txt"), photFile, delimiter=" ", fmt='%0.8f')
    numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"diffExcel.csv"), photFile, delimiter=",", fmt='%0.8f')
        
    # Output astroImageJ file
    outputPeransoCalib=[]
    #i=0
    for i in range(numpy.asarray(photFile).shape[0]):
        outputPeransoCalib.append([photFile[i][0]-2450000.0,photFile[i][1],photFile[i][2]])
        #i=i+1

    numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"diffAIJ.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
    numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"diffAIJ.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')

    # Output replot
    pylab.cla()
    outplotx=photFile[:,0]
    outploty=photFile[:,1]
    pylab.xlabel('BJD')
    pylab.ylabel('Differential ' +filterCode+' Mag')
    pylab.plot(outplotx,outploty,'bo')
    pylab.ylim(max(outploty)+0.02,min(outploty)-0.02,'k-')
    pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
    pylab.grid(True)
    pylab.savefig(os.path.join(outputPath,'V' +str(r+1)+'_'+'EnsembleVarDiffMag.png'))
    pylab.savefig(os.path.join(outputPath,'V' +str(r+1)+'_'+'EnsembleVarDiffMag.eps'))
    pylab.cla()
    pylab.clf()
    pylab.close()
    pylab.close("all")


    if exists:
        measureReject=[]
        for i in range(calibFile.shape[0]):
            if calibFile[i,1] < float(brightD)+calibDiff or calibFile[i,1] > float(brightV)+calibDiff :
                measureReject.append(i)
                print (calibFile[i,1])
                print ("REJECT")
        print (calibFile.shape[0])
        calibFile=numpy.delete(calibFile, measureReject, axis=0)
        print (calibFile.shape[0])

        # Output trimmed files
        numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"calibPeranso.txt"), calibFile, delimiter=" ", fmt='%0.8f')
        numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"calibExcel.csv"), calibFile, delimiter=",", fmt='%0.8f')
            
        # Output astroImageJ file
        outputPeransoCalib=[]
        #i=0
        for i in range(numpy.asarray(calibFile).shape[0]):
            outputPeransoCalib.append([calibFile[i][0]-2450000.0,calibFile[i][1],calibFile[i][2]])
            #i=i+1

        numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"calibAIJ.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
        numpy.savetxt(os.path.join(outcatPath,'V' +str(r+1)+'_'+"calibAIJ.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')
        #r=r+1

    r=r+1
    
sys.exit()
