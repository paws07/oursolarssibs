import numpy
from astropy import units as u
import scipy
from astropy.coordinates import SkyCoord
import glob
import sys
import pylab
import matplotlib
import math
import os
import platform



# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")

if "Mac" or "Darwin" in platform.platform():
    inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
#fileList=glob.glob('*.ps*')


# Load in list of used files
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)


fileList = glob.glob(inputPath)
#print (fileList)
#create directory structure
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

if not os.path.exists(outcatPath):
    os.makedirs(outcatPath)

if not os.path.exists(checkPath):
    os.makedirs(checkPath)

with open("LightcurveStats.txt", "w") as f:
  f.write("Lightcurve Statistics \n\n")

targetFile = numpy.genfromtxt('targetstars.csv', dtype=float, delimiter=',')
# Remove any nan rows from targetFile
targetRejecter=[]
if not (targetFile.shape[0] == 4 and targetFile.size ==4):
    for z in range(targetFile.shape[0]):
      if numpy.isnan(targetFile[z][0]):
        targetRejecter.append(z)
    targetFile=numpy.delete(targetFile, targetRejecter, axis=0)

if targetFile.size == 4 and targetFile.shape[0] ==4:
    loopLength=1
else:
    loopLength=targetFile.shape[0]
for q in range(loopLength):
    exists=os.path.isfile(os.path.join(outcatPath,'V'+str(q+1)+'_calibExcel.csv'))
    if exists:
        calibFile = numpy.genfromtxt(os.path.join(outcatPath,'V'+str(q+1)+'_calibExcel.csv'), dtype=float, delimiter=',')
    else:
        calibFile = numpy.genfromtxt(os.path.join(outcatPath,'V'+str(q+1)+'_diffExcel.csv'), dtype=float, delimiter=',')

    print (targetFile.size)
    print (targetFile.shape[0])

    if targetFile.size == 4 and targetFile.shape[0] ==4:
        period = targetFile[2]
        phaseShift = targetFile[3]
    else:
        period = targetFile[q][2]
        phaseShift = targetFile[q][3]

    #print (calibFile[:,6])
    #print (calibFile[:,10])
    #print (calibFile[:,11])

    # Variable lightcurve

    pylab.cla()
    outplotx=calibFile[:,0]
    outploty=calibFile[:,1]
    print (outplotx)
    print (outploty)
    pylab.xlabel('BJD')
    pylab.ylabel('Apparent ' + str(filterCode) + ' Magnitude')
    pylab.plot(outplotx,outploty,'bo')
    #pylab.plot(linex,liney)
    pylab.ylim(max(outploty)-0.04,min(outploty)+0.04,'k-')
    pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
    pylab.grid(True)
    pylab.savefig(os.path.join(outputPath,'Variable'+str(q+1)+'_' + str(filterCode) +'_Lightcurve.png'))
    pylab.savefig(os.path.join(outputPath,'Variable'+str(q+1)+'_' + str(filterCode) +'_Lightcurve.eps'))

    # Phased lightcurve

    pylab.cla()
    fig = matplotlib.pyplot.gcf()
    outplotx=((calibFile[:,0]/period)+phaseShift)%1
    outploty=calibFile[:,1]
    outplotxrepeat=outplotx+1
    print (outplotx)
    print (outploty)
    pylab.xlabel('Phase')
    pylab.ylabel('Apparent ' + str(filterCode) + ' Magnitude')
    pylab.plot(outplotx,outploty,'bo')
    pylab.plot(outplotxrepeat,outploty,'ro')
    #pylab.plot(linex,liney)
    pylab.ylim(max(outploty)+0.04,min(outploty)-0.04,'k-')
    pylab.xlim(-0.01,2.01)
    pylab.errorbar(outplotx, outploty, yerr=3*calibFile[:,2], fmt='-o', linestyle='None')
    pylab.errorbar(outplotxrepeat, outploty, yerr=3*calibFile[:,2], fmt='-o', linestyle='None')
    pylab.grid(True)
    pylab.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.17, wspace=0.3, hspace=0.4)
    fig.set_size_inches(6,3)
    pylab.savefig(os.path.join(outputPath,'Variable'+str(q+1)+'_' + str(filterCode) +'_PhasedLightcurve.png'))
    pylab.savefig(os.path.join(outputPath,'Variable'+str(q+1)+'_' + str(filterCode) +'_PhasedLightcurve.eps'))

    print ("Variable V"+str(q+1))
    print ("Max Magnitude: "+ str(numpy.max(calibFile[:,1])))
    print ("Min Magnitude: "+ str(numpy.min(calibFile[:,1])))
    print ("Amplitude    : "+ str(numpy.max(calibFile[:,1])-numpy.min(calibFile[:,1])))
    print ("Mid Magnitude: "+ str((numpy.max(calibFile[:,1])+numpy.min(calibFile[:,1]))/2))

    with open("LightcurveStats.txt", "a+") as f:
        f.write("Variable V"+str(q+1)+"\n")
        f.write("Max Magnitude: "+ str(numpy.max(calibFile[:,1]))+"\n")
        f.write("Min Magnitude: "+ str(numpy.min(calibFile[:,1]))+"\n")
        f.write("Amplitude    : "+ str(numpy.max(calibFile[:,1])-numpy.min(calibFile[:,1]))+"\n")
        f.write("Mid Magnitude: "+ str((numpy.max(calibFile[:,1])+numpy.min(calibFile[:,1]))/2)+"\n\n")
    
    

sys.exit()
