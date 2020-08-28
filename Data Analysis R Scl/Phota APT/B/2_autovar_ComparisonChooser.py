import numpy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vo_conesearch import conesearch
from astroquery.vo_conesearch import ConeSearch
#from astroquery.sdss import SDSS
from astroquery.vizier import Vizier
import glob
import sys
import os
import platform

# This is how many standard deviations above the mean to cut off the top.
# The cycle will continue until there are no stars this many std.dev above the mean
stdMultiplier=3
thresholdCounts=1000000000 # UNNECESSARY I THINK! This is the target countrate for the ensemble comparison... the lowest variability stars will be added until this countrate is reached.
variabilityMultiplier=2.5 # The script will keep adding stars until it reaches stars with variability of this value times the minimum variability in the data
removeTargets=1 # Set this to 1 to remove targets from consideration for comparison stars
acceptDistance=1.0 # Furtherest distance in arcseconds for matches
max_sep=acceptDistance * u.arcsec

# Get list of phot files
parentPath = os.getcwd()
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())

# LOAD Phot FILES INTO LIST
photFileArray=[]
for file in fileList:
    photFileArray.append(numpy.genfromtxt(file, dtype=float, delimiter=','))
photFileArray=numpy.asarray(photFileArray)

#Grab the candidate comparison stars
compFile = numpy.genfromtxt('screenedComps.csv', dtype=float, delimiter=',')

# Remove targets from consideration for target stars (for obvious reasons!
if removeTargets == 1:
    print ("Removing Target Stars from potential Comparisons")
    targetFile = numpy.genfromtxt('targetstars.csv', dtype=float, delimiter=',')
    fileRaDec = SkyCoord(ra=compFile[:,0]*u.degree, dec=compFile[:,1]*u.degree)
    #print (len(targetFile))
    #print (targetFile)

    # Remove any nan rows from targetFile
    targetRejecter=[]
    if not (targetFile.shape[0] == 4 and targetFile.size ==4):
        for z in range(targetFile.shape[0]):
          if numpy.isnan(targetFile[z][0]):
            targetRejecter.append(z)
        targetFile=numpy.delete(targetFile, targetRejecter, axis=0)

    # Remove targets from consideration    
    if len(targetFile)== 4:
        loopLength=1
    else:
        loopLength=targetFile.shape[0]
    targetRejects=[]
    for q in range(loopLength):
        print (targetFile.shape[0])
        
        if int(len(targetFile)) == 4:
          varCoord = SkyCoord(targetFile[0],(targetFile[1]), frame='icrs', unit=u.deg) # Need to remove target stars from consideration
        else:
          varCoord = SkyCoord(targetFile[q][0],(targetFile[q][1]), frame='icrs', unit=u.deg) # Need to remove target stars from consideration
        #print (varCoord)
        #print (fileRaDec)
        idx, d2d, _ = varCoord.match_to_catalog_sky(fileRaDec)
        print (d2d.arcsecond)
        if d2d.arcsecond < acceptDistance:
          targetRejects.append(idx)

    #Remove target and restore skycoord list
    compFile=numpy.delete(compFile, idx, axis=0)
    fileRaDec = SkyCoord(ra=compFile[:,0]*u.degree, dec=compFile[:,1]*u.degree)

print (compFile.shape[0])


# Get Average RA and Dec from file
if compFile.shape[0] == 13:
    print (compFile[0])
    print (compFile[1])
    avgCoord=SkyCoord(ra=(compFile[0])*u.degree, dec=(compFile[1]*u.degree))

else:
    print (numpy.average(compFile[:,0]))
    print (numpy.average(compFile[:,1]))
    avgCoord=SkyCoord(ra=(numpy.average(compFile[:,0]))*u.degree, dec=(numpy.average(compFile[:,1]))*u.degree)


# Check VSX for any known variable stars and remove them from the list
variableResult=Vizier.query_region(avgCoord, '0.33 deg', catalog='VSX')['B/vsx/vsx']

print (variableResult)

variableResult=variableResult.to_pandas()

print (variableResult.keys())

variableSearchResult=variableResult[['RAJ2000','DEJ2000']].as_matrix()


raCat=variableSearchResult[:,0]
print (raCat)
decCat=variableSearchResult[:,1]
print (decCat)

varStarReject=[]
for t in range(raCat.size):
  print (raCat[t])
  compCoord=SkyCoord(ra=raCat[t]*u.degree, dec=decCat[t]*u.degree)
  print (compCoord)
  catCoords=SkyCoord(ra=compFile[:,0]*u.degree, dec=compFile[:,1]*u.degree)
  idxcomp,d2dcomp,d3dcomp=compCoord.match_to_catalog_sky(catCoords)
  print (d2dcomp)
  if d2dcomp *u.arcsecond < max_sep*u.arcsecond:
    print ("match!")
    varStarReject.append(t)
  else:
    print ("no match!")
    

print ("Number of stars prior to VSX reject")
print (compFile.shape[0])
compFile=numpy.delete(compFile, varStarReject, axis=0)
print ("Number of stars post to VSX reject")
print (compFile.shape[0])


if (compFile.shape[0] ==1):
    print ("Looks like you have a single comparison star!")
    compFile=[[compFile[0][0],compFile[0][1],0.01]]
    compFile=numpy.asarray(compFile)
    numpy.savetxt(os.path.join(parentPath,"compsUsed.csv"), compFile, delimiter=",", fmt='%0.8f')
    sortStars=[[compFile[0][0],compFile[0][1],0.01,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]]
    sortStars=numpy.asarray(sortStars)
    numpy.savetxt("stdComps.csv", sortStars, delimiter=",", fmt='%0.8f')
    sys.exit()
    

while True:
    # First half of Loop: Add up all of the counts of all of the comparison stars
    # To create a gigantic comparison star.
    fileCount=[]

    #print ("Please wait... calculating ensemble comparison star for each image")
    for imgs in range(photFileArray.shape[0]):
        allCounts=0.0
        photFile = photFileArray[imgs]
        fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)
        for j in range(compFile.shape[0]):        
            matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
            idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)
            allCounts=numpy.add(allCounts,photFile[idx][4])                    
        print (fileList[imgs] + "\nTotal Counts in Image: " + str(allCounts) +"\n*")
        fileCount=numpy.append(fileCount, allCounts)

    # Second half of Loop: Calculate the variation in each candidate comparison star in brightness
    # compared to this gigantic comparison star.
    rejectStar=[]
    stdCompStar=[]
    sortStars=[]
    for j in range(compFile.shape[0]):
        compDiffMags=[]
        q=0
        print ("*************************")
        print ("RA : " + str(compFile[j][0]))
        print ("DEC: " + str(compFile[j][1]))
        for imgs in range(photFileArray.shape[0]):
            #print file
            photFile = photFileArray[imgs]
            fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)                      
            matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
            idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)
            compDiffMags=numpy.append(compDiffMags,2.5 * numpy.log10(photFile[idx][4]/fileCount[q]))
            q=numpy.add(q,1)

        print ("VAR: " +str(numpy.std(compDiffMags)))
        stdCompStar.append(numpy.std(compDiffMags))
        sortStars.append([compFile[j][0],compFile[j][1],numpy.std(compDiffMags),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

    # Calculate and present the sample statistics
    print (fileCount)
    print (stdCompStar)
    stdCompMed=numpy.median(stdCompStar)
    print (numpy.median(stdCompStar))
    stdCompStd=numpy.std(stdCompStar)
    print (numpy.std(stdCompStar))

    # Delete comparisons that have too high a variability
    starRejecter=[]
    for j in range(len(stdCompStar)):
        print (stdCompStar[j])
        if ( stdCompStar[j] > (stdCompMed + (stdMultiplier*stdCompStd)) ):
            print ("Star Rejected, Variability too high!")
            starRejecter.append(j)
            #potVar.append([compFile[j][0],compFile[j][1],stdCompStar[j]])
        if ( numpy.isnan(stdCompStar[j]) ) :
            print ("Star Rejected, Invalid Entry!")
            starRejecter.append(j)
    compFile=numpy.delete(compFile, starRejecter, axis=0)
    sortStars=numpy.delete(sortStars, starRejecter, axis=0)

    # Calculate and present statistics of sample of candidate comparison stars.
    print ("Median variability")
    print (numpy.median(stdCompStar))
    print ("Std variability")
    print (numpy.std(stdCompStar))
    print ("Min variability")
    print (numpy.min(stdCompStar))
    print ("Max variability")
    print (numpy.max(stdCompStar))
    print ("Number of Stable Comparison Candidates")
    print ((compFile.shape[0]))

    # Once we have stopped rejecting stars, this is our final candidate catalogue
    if (starRejecter == []):

        print (' ')
        print ('Statistical stability reached.')
        print (' ')
        print ('List of stable comparison candidates output to stdComps.csv')
        
        numpy.savetxt("stdComps.csv", sortStars, delimiter=",", fmt='%0.8f')

        # GET REFERENCE IMAGE
        # Sort through and find the largest file and use that as the reference file
        fileSizer=0
        print ("Finding image with most stars detected")
        for imgs in range(photFileArray.shape[0]):
            photFile = photFileArray[imgs]
            if photFile.size > fileSizer:
                referenceFrame=photFile
                print (photFile.size)
                fileSizer=photFile.size  
        print ("Setting up reference Frame")
        fileRaDec = SkyCoord(ra=referenceFrame[:,0]*u.degree, dec=referenceFrame[:,1]*u.degree)
        
        # SORT THE COMP FILE such that least variable comparison is first
        sortStars=(sortStars[sortStars[:,2].argsort()])

        # PICK COMPS UNTIL OVER THE THRESHOLD OF COUNTS OR VRAIABILITY ACCORDING TO REFERENCE IMAGE
        print ("PICK COMPS UNTIL OVER THE THRESHOLD ACCORDING TO REFERENCE IMAGE")

        variabilityMax=(numpy.min(stdCompStar)*variabilityMultiplier)
        
        compFile=[]
        tempCountCounter=0.0
        finalCountCounter=0.0
        for j in range(sortStars.shape[0]):        
            matchCoord=SkyCoord(ra=sortStars[j][0]*u.degree, dec=sortStars[j][1]*u.degree)        
            idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)
            tempCountCounter=numpy.add(tempCountCounter,referenceFrame[idx][4])

            
            if tempCountCounter < thresholdCounts:
                if sortStars[j][2] < variabilityMax: 
                    compFile.append([sortStars[j][0],sortStars[j][1],sortStars[j][2]])
                    print ("Comp " + str(j+1) + " std: " + str(sortStars[j][2]))
                    print ("Cumulative Counts thus far: " + str(tempCountCounter))
                    finalCountCounter=numpy.add(finalCountCounter,referenceFrame[idx][4])

              

        print ("Selected stars listed below:")
        print (compFile)
        
        print ("Final Ensemble Counts: " + str(finalCountCounter))
        compFile=numpy.asarray(compFile)

        print (str(compFile.shape[0]) + " Stable Comparison Candidates below variability threshold output to compsUsed.csv")
        print ("Variabilities of these stars")
        print (compFile[:,2])

        numpy.savetxt(os.path.join(parentPath,"compsUsed.csv"), compFile, delimiter=",", fmt='%0.8f')

        print ("\nStars rejected due to being identified variables in VSX: " + str(len(varStarReject)))
        sys.exit()


sys.exit()
