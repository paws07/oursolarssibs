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

#filterCode = 3 # u=0, g=1, r=2, i=3, z=4
#calibFlag = 0 # 0 = no calibration attempted, 1 = calibration attempted.
#thresholdCounts=1000000 # This is the target countrate for the ensemble comparison... the lowest variability stars will be added until this countrate is reached.
#variabilityMax=0.025 # This will stop adding ensemble comparisons if it starts using stars higher than this variability
minimumVariableCounts=10000  # Do not try to detect variables dimmer than this.
acceptDistance=1.0 # Furtherest distance in arcseconds for matches
minimumNoOfObs=30 # Minimum number of observations to count as a potential variable.
targetIDonly=1 # Do this just to produce a list of all possible stars to analyse
outputToTargetStars=0 # makes this output the list to targetstars.csv ... not implemented yet.

# Initialise Directories
parentPath = os.getcwd()
if "Windows" in platform.platform():
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")

if "Mac" or "Darwin" in platform.platform():
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")

#create directory structure
if not os.path.exists(outputPath):
    os.makedirs(outputPath)

if not os.path.exists(outcatPath):
    os.makedirs(outcatPath)

if not os.path.exists(checkPath):
    os.makedirs(checkPath)

# Load in list of used files
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())
    
# LOAD Phot FILES INTO LIST

photFileArray=[]
for file in fileList:
    photFileArray.append(numpy.genfromtxt(file, dtype=float, delimiter=','))

photFileArray=numpy.asarray(photFileArray)

# LOAD IN COMPARISON FILE
preFile = numpy.genfromtxt('stdComps.csv', dtype=float, delimiter=',')
if preFile.shape[0] !=13:
    preFile=(preFile[preFile[:,2].argsort()])# GET REFERENCE IMAGE
# Sort through and find the largest file and use that as the reference file
fileSizer=0
print ("Finding image with most stars detected")
for imgs in range(photFileArray.shape[0]):
    photFile = photFileArray[imgs]
    if photFile.size > fileSizer:
        referenceFrame=photFile
        print (photFile.size)
        fileSizer=photFile.size  

compFile=numpy.genfromtxt(os.path.join(parentPath,"compsUsed.csv"), dtype=float, delimiter=',')
print ("Stable Comparison Candidates below variability threshold")
outputPhot=[]

# Get total counts for each file

fileCount=[]
compArray=[]
allCountsArray=[]

if (targetIDonly != 1): 
    for imgs in range(photFileArray.shape[0]):
        allCounts=0.0
        allCountsErr=0.0
        photFile = photFileArray[imgs]
        fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)
        #Array of comp measurements

        print ("***************************************")
        print ("Calculating total Comparison counts for")
        print (fileList[imgs])

        for j in range(compFile.shape[0]):        

            if compFile.size == 2:
                matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
            else:
                matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
            idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)

            allCounts=allCounts+photFile[idx][4]
            allCountsErr=allCountsErr+photFile[idx][5]

        allCountsArray.append([allCounts,allCountsErr])

    print (allCountsArray)


# Define targetlist as every star in referenceImage above a count threshold
print ("Setting up Variable Search List")
targetFile=referenceFrame
# Although remove stars that are below the variable countrate
starReject=[]
for q in range(targetFile.shape[0]):
    if targetFile[q][4] < minimumVariableCounts:
        starReject.append(q)
print ("Total number of stars in reference Frame")
print (targetFile.shape[0])
targetFile=numpy.delete(targetFile, starReject, axis=0)
print ("Total number of stars with sufficient counts")
print (targetFile.shape[0])

## NEED TO REMOVE COMPARISON STARS FROM TARGETLIST
if (targetIDonly != 1):

    allcountscount=0
    # For each variable calculate the variability
    outputVariableHolder=[]
    for q in range(targetFile.shape[0]):
        print ("***********************************************************************\nProcessing Variable " +str(q+1)+ "\nRA " +str(targetFile[q][0]) + " DEC " + str(targetFile[q][1]))
        varCoord = SkyCoord(targetFile[q][0],(targetFile[q][1]), frame='icrs', unit=u.deg) # Need to remove target stars from consideration
        outputPhot=[] # new
        compArray=[]
        compList=[]


        diffMagHolder=[]
        
        allcountscount=0

        for imgs in range(photFileArray.shape[0]):
            compList=[]
            photFile = photFileArray[imgs]
        
            fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)

            idx, d2d, d3d = varCoord.match_to_catalog_sky(fileRaDec)
            if (numpy.less(d2d.arcsecond, acceptDistance) and ((numpy.multiply(-2.5,numpy.log10(numpy.divide(photFile[idx][4],allCountsArray[allcountscount][0])))) != numpy.inf )):
                
                diffMagHolder=numpy.append(diffMagHolder,(numpy.multiply(-2.5,numpy.log10(numpy.divide(photFile[idx][4],allCountsArray[allcountscount][0])))))
                
            allcountscount=numpy.add(allcountscount,1)
            
        ## REMOVE MAJOR OUTLIERS FROM CONSIDERATION
        while True:
            stdVar=numpy.std(diffMagHolder)
            avgVar=numpy.average(diffMagHolder)
            starReject=[]
            z=0
            for j in range(numpy.asarray(diffMagHolder).shape[0]):
                if diffMagHolder[j] > avgVar+(4*stdVar) or diffMagHolder[j] < avgVar-(4*stdVar) :
                    starReject.append(j)
                    print ("REJECT " + str(diffMagHolder[j]))
                    z=1
            diffMagHolder=numpy.delete(diffMagHolder, starReject, axis=0)
            if z==0:
                break
            

        print ("Standard Deviation in mag")
        print (numpy.std(diffMagHolder))
        print ("Median Magnitude")
        print (numpy.median(diffMagHolder))
        print ("Number of Observations")
        print (numpy.asarray(diffMagHolder).shape[0])
        if (  numpy.asarray(diffMagHolder).shape[0] > minimumNoOfObs):
            outputVariableHolder.append( [targetFile[q][0],targetFile[q][1],numpy.median(diffMagHolder), numpy.std(diffMagHolder), numpy.asarray(diffMagHolder).shape[0]])

    print ("Outputting list of the variability of all possible target stars to starVariability.csv")
    print ("Please transfer the RA and Decs of your target stars into targetstars.csv")
    numpy.savetxt(os.path.join(parentPath,"starVariability.csv"), outputVariableHolder, delimiter=",", fmt='%0.8f')

else:
    print ("Outputting list of all possible target stars to possibleTargets.csv")
    print ("Please transfer the RA and Decs of your target stars into targetstars.csv")
    print (targetFile)
    numpy.savetxt(os.path.join(parentPath,"possibleTargets.csv"), targetFile, delimiter=",", fmt='%0.8f')
  
sys.exit()
