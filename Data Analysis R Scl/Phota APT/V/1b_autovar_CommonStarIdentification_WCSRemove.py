import numpy
from astropy import units as u
from astropy.coordinates import SkyCoord
import glob
import sys
import os
import platform

#Initialisation values
acceptDistance=1.0 # Furtherest distance in arcseconds for matches
minimumCounts=2000 # look for comparisons brighter than this
maximumCounts=1000000 # look for comparisons dimmer than this
imageFracReject=0.5 # This is a value which will reject images based on number of stars detected
starFracReject=0.1 # This ia a value which will reject images that reject this fraction of available stars after....
rejectStart=7 # This many initial images (lots of stars are expected to be rejected in the early images)
minCompStars=1 # This is the minimum number of comp stars required
usedImages=[] # This holds the list of acceptable files to use. 

# Get list of phot files
parentPath = os.getcwd()
print (platform.platform())
if "Windows" in platform.platform():
    inputPath = os.path.join(parentPath,"inputs\\*.*")
if "Mac" or "Darwin" in platform.platform():
    inputPath = os.path.join(parentPath,"inputs//*.*")
fileList = glob.glob(inputPath)

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)

# Check that all of the files have the same filter
for file in fileList:
    filterCheck=(file.replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
    if filterCheck != filterCode:
        
        print ("Check your images, the script detected multiple filters in your file list. Autovar currently only does one filter at a time.")
        sys.exit()


# Sort through and find the largest file and use that as the reference file
fileSizer=0
print ("Finding image with most stars detected")
for file in fileList:
    print (file)
    photFile = numpy.genfromtxt(file, dtype=float, delimiter=',', encoding="utf8")
    if photFile.size > fileSizer:
        if (( numpy.asarray(photFile[:,0]) > 360).sum() == 0) and ( numpy.asarray(photFile[0][0]) != 'null') and ( numpy.asarray(photFile[0][0]) != 0.0) :
            referenceFrame=photFile
            print (photFile.size)
            fileSizer=photFile.size
            print (file)
print ("Setting up reference Frame")
fileRaDec = SkyCoord(ra=referenceFrame[:,0]*u.degree, dec=referenceFrame[:,1]*u.degree)
print ("Removing stars with low or high counts")
rejectStars=[]
# Check star has adequate counts
for j in range(referenceFrame.shape[0]):
    if ( referenceFrame[j][4] < minimumCounts or referenceFrame[j][4] > maximumCounts ):
        rejectStars.append(int(j))
print ("Number of stars prior")
print (referenceFrame.shape[0])
referenceFrame=numpy.delete(referenceFrame, rejectStars, axis=0)
print ("Number of stars post")
print (referenceFrame.shape[0])


imgsize=imageFracReject * fileSizer # set threshold size
rejStartCounter=0
imgReject=0 # Number of images rejected due to high rejection rate
loFileReject=0
wcsFileReject=0
for file in fileList:
    rejStartCounter=rejStartCounter +1
    photFile = numpy.genfromtxt(file, dtype=float, delimiter=',')

    if (photFile.size > 7):
        if (( numpy.asarray(photFile[:,0]) > 360).sum() == 0) and ( numpy.asarray(photFile[0][0]) != 'null') and ( numpy.asarray(photFile[0][0]) != 0.0) :
            fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)

            print (' ')
            print ('Image Number: ' + str(rejStartCounter))
            print (file)
            print ("Image treshold size: "+str(imgsize))
            print ("Image catalogue size: "+str(photFile.size))
            if photFile.size > imgsize:
            
                # Checking existance of stars in all photometry files
                rejectStars=[] # A list to hold what stars are to be rejected

                # Find whether star in reference list is in this phot file, if not, reject star.
                for j in range(referenceFrame.shape[0]):
                    photRAandDec=SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)
                    testStar=SkyCoord(ra=referenceFrame[j][0]*u.degree, dec=referenceFrame[j][1]*u.degree)
                    idx, d2d, d3d = testStar.match_to_catalog_sky(photRAandDec)
                    if (d2d.arcsecond > acceptDistance):
                        #"No Match! Nothing within range."
                        rejectStars.append(int(j))

                
                # if the rejectstar list is not empty, remove the stars from the reference List
                if rejectStars != []:
                    
                    if not (((len(rejectStars) / referenceFrame.shape[0]) > starFracReject) and rejStartCounter > rejectStart):
                        referenceFrame=numpy.delete(referenceFrame, rejectStars, axis=0)
                        print ('**********************') 
                        print ('Stars Removed  : ' +str(len(rejectStars))) 
                        print ('Remaining Stars: ' +str(referenceFrame.shape[0]))
                        print ('**********************')
                        usedImages.append(file)
                    else:
                        print ('**********************')
                        print ('Image Rejected due to too high a fraction of rejected stars')
                        print (len(rejectStars) / referenceFrame.shape[0])
                        print ('**********************')
                        imgReject=imgReject+1
                else:
                    print ('**********************') 
                    print ('All Stars Present')
                    print ('**********************')
                    usedImages.append(file)

                # If we have removed all stars, we have failed!
                if (referenceFrame.shape[0]==0):
                    print ("All Stars Removed. Try removing problematic files or raising the imageFracReject")
                    sys.exit()

                if (referenceFrame.shape[0]< minCompStars):
                    print ("There are less than the requested number of Comp Stars. Try removing problematic files or raising the imageFracReject")
                    sys.exit()
        else:
            print ('**********************') 
            print ("WCS Coordinates broken")
            print ('**********************')
            wcsFileReject=wcsFileReject+1
    else:
        print ('**********************') 
        print ("CONTAINS TOO FEW STARS")
        print ('**********************')
        loFileReject=loFileReject+1
      


outputComps=[]
for j in range (referenceFrame.shape[0]):
    outputComps.append([referenceFrame[j][0],referenceFrame[j][1]])

print ("These are the identified common stars of sufficient brightness that are in every image")
print (outputComps)

print  (' ')
print ('Images Rejected due to high star rejection rate: ' + str(imgReject))
print ('Images Rejected due to low file size: ' + str(loFileReject))
print ('Images Rejected due to broken WCS: ' + str(wcsFileReject))
print ('Out of this many original images: ' + str(len(fileList)))

print  (' ')

print ("Number of candidate Comparison Stars Detected: " + str(len(outputComps)))
print  (' ')
print  ('Output sent to screenedComps.csv ready for use in CompDeviation')
numpy.savetxt("screenedComps.csv", outputComps, delimiter=",", fmt='%0.8f')
print  ('UsedImages ready for use in 2_CompChooser')
with open("usedImages.txt", "w") as f:
    for s in usedImages:
        f.write(str(s) +"\n")
print  (' ')

print  ('If you are going to enter targets into targetstars.csv, do this now')

sys.exit()
