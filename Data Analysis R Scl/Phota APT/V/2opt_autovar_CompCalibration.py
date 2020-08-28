import numpy
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vo_conesearch import conesearch
from astroquery.vo_conesearch import ConeSearch
from astroquery.sdss import SDSS
from astroquery.vizier import Vizier
import glob
import pandas
import sys
import os
import platform
import math
import pylab

Vizier.ROW_LIMIT = -1
max_sep=1.0 * u.arcsec
max_magerr=0.05

variabilityMultiplier=2.0 # The script will keep adding stars until it reaches stars with variability of this value times the minimum variability in the data

stdMultiplier=2

panStarrsInstead=0 # If set to 1 this will use Panstarrs rather than SDSS for northern targets. SDSS tends to have a terrible z filter with terrible fringing, so this was incorporated to use panstarss for z band images.

# Get list of phot files
parentPath = os.getcwd()

if "Windows" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")
    calibPath = os.path.join(parentPath,"calibcats\\")

if "Mac" or "Darwin" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
    calibPath = os.path.join(parentPath,"calibcats//")


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

if not os.path.exists(calibPath):
    os.makedirs(calibPath)

# Get List of Files Used
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)

# Load compsused
compFile = numpy.genfromtxt('stdComps.csv', dtype=float, delimiter=',')
print (compFile.shape[0])

if compFile.shape[0] == 13:
    compCoords=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
else:
    compCoords=SkyCoord(ra=compFile[:,0]*u.degree, dec=compFile[:,1]*u.degree)

# Get Average RA and Dec from file
if compFile.shape[0] == 13:
    print (compFile[0])
    print (compFile[1])
    avgCoord=SkyCoord(ra=(compFile[0])*u.degree, dec=(compFile[1]*u.degree))

else:
    print (numpy.average(compFile[:,0]))
    print (numpy.average(compFile[:,1]))
    avgCoord=SkyCoord(ra=(numpy.average(compFile[:,0]))*u.degree, dec=(numpy.average(compFile[:,1]))*u.degree)

# get results from internetz

if filterCode=='B' or filterCode=='V':
    

    #collect APASS results
    apassResult=Vizier.query_region(avgCoord, '0.33 deg', catalog='APASS')['II/336/apass9']
    
    print (apassResult)
    apassResult=apassResult.to_pandas()
    #print (apassResult.to_pandas())

    print (apassResult.keys())

    apassSearchResult=apassResult[['RAJ2000','DEJ2000','Bmag','e_Bmag','Vmag','e_Vmag']].as_matrix()

    print (apassSearchResult)

    #apassSearchResult=apassSearchResult.values

    raCat=apassSearchResult[:,0]
    print (raCat)
    decCat=apassSearchResult[:,1]
    print (decCat)
    if filterCode=='B':
        magCat=apassSearchResult[:,2]
        print (magCat)
        emagCat=apassSearchResult[:,3]
        print (emagCat)
        
    if filterCode=='V':
        magCat=apassSearchResult[:,4]
        print (magCat)
        emagCat=apassSearchResult[:,5]
        print (emagCat)


if filterCode=='up' or filterCode=='gp' or filterCode=='rp' or filterCode=='ip' or filterCode=='zs':
    # Are there entries in SDSS?
    sdssResult=SDSS.query_region(avgCoord, '0.33 deg')
    #print(sdssResult)
    sdssFind=1

    # If not in SDSS, try Skymapper
    if sdssResult==None:
        sdssFind=0
        print ("Not found in SDSS, must be in the South.")
        #print (ConeSearch.URL)
        ConeSearch.URL='http://skymapper.anu.edu.au/sm-cone/aus/query?'
        sdssResult=ConeSearch.query_region(avgCoord, '0.33 deg')

        print (sdssResult)

        searchResult=sdssResult.to_table().to_pandas()

        print (searchResult)

        sdssSearchResult=searchResult[['raj2000','dej2000','u_psf','e_u_psf','g_psf','e_g_psf','r_psf','e_r_psf','i_psf','e_i_psf','z_psf','e_z_psf']].as_matrix()

        print (sdssSearchResult[:,0])

        raCat=sdssSearchResult[:,0]
        print (raCat)
        decCat=sdssSearchResult[:,1]
        print (decCat)
        if filterCode=='up':
            magCat=sdssSearchResult[:,2]
            print (magCat)
            emagCat=sdssSearchResult[:,3]
            print (emagCat)
        if filterCode=='gp':
            magCat=sdssSearchResult[:,4]
            print (magCat)
            emagCat=sdssSearchResult[:,5]
            print (emagCat)
        if filterCode=='rp':
            magCat=sdssSearchResult[:,6]
            print (magCat)
            emagCat=sdssSearchResult[:,7]
            print (emagCat)
        if filterCode=='ip':
            magCat=sdssSearchResult[:,8]
            print (magCat)
            emagCat=sdssSearchResult[:,9]
            print (emagCat)
        if filterCode=='zs':
            magCat=sdssSearchResult[:,10]
            print (magCat)
            emagCat=sdssSearchResult[:,11]
            print (emagCat)

    elif panStarrsInstead ==1 and filterCode!='up':
        print ("Panstarrs!")
        
        sdssResult=Vizier.query_region(avgCoord, '0.33 deg', catalog='PanStarrs')['II/349/ps1']
        print (sdssResult)
        print (sdssResult.keys())
        

        searchResult=sdssResult.to_pandas()

        print (searchResult)

        sdssSearchResult=searchResult[['RAJ2000','DEJ2000','gmag','e_gmag','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag']].as_matrix()

        print (sdssSearchResult[:,0])

        raCat=sdssSearchResult[:,0]
        print (raCat)
        decCat=sdssSearchResult[:,1]
        print (decCat)
        if filterCode=='up':
            magCat=sdssSearchResult[:,2]
            print (magCat)
            emagCat=sdssSearchResult[:,3]
            print (emagCat)
        if filterCode=='gp':
            magCat=sdssSearchResult[:,4]
            print (magCat)
            emagCat=sdssSearchResult[:,5]
            print (emagCat)
        if filterCode=='rp':
            magCat=sdssSearchResult[:,6]
            print (magCat)
            emagCat=sdssSearchResult[:,7]
            print (emagCat)
        if filterCode=='ip':
            magCat=sdssSearchResult[:,8]
            print (magCat)
            emagCat=sdssSearchResult[:,9]
            print (emagCat)
        if filterCode=='zs':
            magCat=sdssSearchResult[:,10]
            print (magCat)
            emagCat=sdssSearchResult[:,11]
            print (emagCat)


    else:
        print ("goodo lets do the sdss stuff then.")
        sdssResult=Vizier.query_region(avgCoord, '0.33 deg', catalog='SDSS')['V/147/sdss12']
        print (sdssResult)
        print (sdssResult.keys())

        searchResult=sdssResult.to_pandas()

        print (searchResult)

        sdssSearchResult=searchResult[['RA_ICRS','DE_ICRS','umag','e_umag','gmag','e_gmag','rmag','e_rmag','imag','e_imag','zmag','e_zmag']].as_matrix()

        print (sdssSearchResult[:,0])

        raCat=sdssSearchResult[:,0]
        print (raCat)
        decCat=sdssSearchResult[:,1]
        print (decCat)
        if filterCode=='up':
            magCat=sdssSearchResult[:,2]
            print (magCat)
            emagCat=sdssSearchResult[:,3]
            print (emagCat)
        if filterCode=='gp':
            magCat=sdssSearchResult[:,4]
            print (magCat)
            emagCat=sdssSearchResult[:,5]
            print (emagCat)
        if filterCode=='rp':
            magCat=sdssSearchResult[:,6]
            print (magCat)
            emagCat=sdssSearchResult[:,7]
            print (emagCat)
        if filterCode=='ip':
            magCat=sdssSearchResult[:,8]
            print (magCat)
            emagCat=sdssSearchResult[:,9]
            print (emagCat)
        if filterCode=='zs':
            magCat=sdssSearchResult[:,10]
            print (magCat)
            emagCat=sdssSearchResult[:,11]
            print (emagCat)

#Setup standard catalogue coordinates
catCoords=SkyCoord(ra=raCat*u.degree, dec=decCat*u.degree)


#Get calib mags for least variable IDENTIFIED stars.... not the actual stars in compUsed!! Brighter, less variable stars may be too bright for calibration!
#So the stars that will be used to calibrate the frames to get the OTHER stars.
calibStands=[]
if compFile.shape[0] ==13:
    lenloop=1
else:
    lenloop=len(compFile[:,0])
for q in range(lenloop):
    if compFile.shape[0] ==13:
        compCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
    else:
        compCoord=SkyCoord(ra=compFile[q][0]*u.degree, dec=compFile[q][1]*u.degree)
    idxcomp,d2dcomp,d3dcomp=compCoord.match_to_catalog_sky(catCoords)

    if d2dcomp *u.arcsecond < max_sep*u.arcsecond:
        if not numpy.isnan(magCat[idxcomp]):
            #print (idxcomp)
            #print (d2dcomp)
            #print (magCat[idxcomp])
            #print (emagCat[idxcomp])

            if compFile.shape[0] ==13:
                calibStands.append([compFile[0],compFile[1],compFile[2],magCat[idxcomp],emagCat[idxcomp]])
            else:
                calibStands.append([compFile[q][0],compFile[q][1],compFile[q][2],magCat[idxcomp],emagCat[idxcomp]])

# Get the set of least variable stars to use as a comparison to calibrate the files (to eventually get the *ACTUAL* standards
#print (numpy.asarray(calibStands).shape[0])
if numpy.asarray(calibStands).shape[0] == 0:
    print ("We could not find a suitable match between any of your stars and the calibration catalogue")
    print ("You might need to reduce the low value (usually 10000) to get some dimmer stars in script 1")
    print ("Perhaps try 5000 then 1000. You are trying to find dim stars to calibrate to.")
    sys.exit()

varimin=(numpy.min(numpy.asarray(calibStands)[:,2])) * variabilityMultiplier

print ("varimin")
#print (numpy.asarray(calibStands)[:,2])
print (varimin)

calibStandsReject=[]
for q in range(len(numpy.asarray(calibStands)[:,0])):
    if calibStands[q][2] > varimin:
        calibStandsReject.append(q)
        #print (calibStands[q][2])

calibStands=numpy.delete(calibStands, calibStandsReject, axis=0) 

calibStand=numpy.asarray(calibStands)  

numpy.savetxt("calibStands.csv", calibStands , delimiter=",", fmt='%0.8f')






# Lets use this set to calibrate each datafile and pull out the calibrated compsused magnitudes

compUsedFile = numpy.genfromtxt('compsUsed.csv', dtype=float, delimiter=',')

calibCompUsed=[]

print ("CALIBRATING EACH FILE")
for file in fileList:

    print (file)

    #Get the phot file into memory
    photFile = numpy.genfromtxt(file, dtype=float, delimiter=',')
    photCoords=SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)

    

    #Convert the phot file into instrumental magnitudes

    for r in range(len(photFile[:,0])):
        #print (photFile[r,5])
        #print (photFile[r,4])
        if photFile[r,4] == 0.0:
            photFile[r,4]=numpy.nan
            photFile[r,5]=numpy.nan
        else:
            photFile[r,5]=1.0857 * (photFile[r,5]/photFile[r,4])
            photFile[r,4]=-2.5*math.log10(photFile[r,4])        

    #Pull out the CalibStands out of each file
    tempDiff=[]
    for q in range(len(calibStands[:,0])):
        calibCoord=SkyCoord(ra=calibStand[q][0]*u.degree,dec=calibStand[q][1]*u.degree)
        idx,d2d,d3d=calibCoord.match_to_catalog_sky(photCoords)
        tempDiff.append(calibStand[q,3]-photFile[idx,4])

    #print (tempDiff)
    tempZP= (numpy.median(tempDiff))
    #print (numpy.std(tempDiff))

    #Shift the magnitudes in the phot file by the zeropoint
    for r in range(len(photFile[:,0])):
        photFile[r,4]=photFile[r,4]+tempZP                    

    #Save the calibrated photfiles to the calib directory
    googFile = (file.replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',""))
    numpy.savetxt(os.path.join(calibPath ,googFile.replace('.','calibrated.')), photFile, delimiter=",", fmt='%0.8f')
    #numpy.savetxt(os.path.join(calibPath ,googFile.replace('.','CALIB.')), calFile, delimiter=",", fmt='%0.8f')

    #Look within photfile for ACTUAL usedcomps.csv and pull them out
    lineCompUsed=[]
    if compUsedFile.shape[0] ==3 and compUsedFile.size == 3:
        lenloop=1
    else:
        lenloop=len(compUsedFile[:,0])

    #print (compUsedFile.size)
    for r in range(lenloop):
        if compUsedFile.shape[0] ==3 and compUsedFile.size ==3:
            compUsedCoord=SkyCoord(ra=compUsedFile[0]*u.degree,dec=compUsedFile[1]*u.degree)
        else:
            compUsedCoord=SkyCoord(ra=compUsedFile[r][0]*u.degree,dec=compUsedFile[r][1]*u.degree)
        idx,d2d,d3d=compUsedCoord.match_to_catalog_sky(photCoords)
        lineCompUsed.append(photFile[idx,4])

    #print (lineCompUsed)
    calibCompUsed.append(lineCompUsed)


# Finalise calibcompsusedfile
#print (calibCompUsed)

calibCompUsed=numpy.asarray(calibCompUsed)
#print (calibCompUsed[0,:])

finalCompUsedFile=[]
sumStd=[]
for r in range(len(calibCompUsed[0,:])):
    #Calculate magnitude and stdev
    #print (calibCompUsed[:,r])
    #print (numpy.median(calibCompUsed[:,r]))
    #print (numpy.std(calibCompUsed[:,r]))
    #sumStd=sumStd+numpy.std(calibCompUsed[:,r])
    sumStd.append(numpy.std(calibCompUsed[:,r]))
    #print (calibCompUsed[:,r])
    #print (numpy.std(calibCompUsed[:,r]))
    if compUsedFile.shape[0] ==3  and compUsedFile.size ==3:
        finalCompUsedFile.append([compUsedFile[0],compUsedFile[1],compUsedFile[2],numpy.median(calibCompUsed[:,r]),numpy.asarray(calibStands[0])[4]])
    else:
        finalCompUsedFile.append([compUsedFile[r][0],compUsedFile[r][1],compUsedFile[r][2],numpy.median(calibCompUsed[:,r]),numpy.std(calibCompUsed[:,r])])

#print (finalCompUsedFile)
print (" ")
sumStd=numpy.asarray(sumStd)

errCalib = numpy.median(sumStd) / pow((len(calibCompUsed[0,:])), 0.5)

#print (len(calibCompUsed[0,:]))
if len(calibCompUsed[0,:]) == 1:
    print ("As you only have one comparison, the uncertainty in the calibration is unclear")
    print ("But we can take the catalogue value, although we should say this is a lower uncertainty")
    print ("Error/Uncertainty in Calibration: " +str(numpy.asarray(calibStands[0])[4]))
else:
    print ("Median Standard Deviation of any one star: " + str(numpy.median(sumStd)))
    print ("Standard Error/Uncertainty in Calibration: " +str(errCalib))

with open("calibrationErrors.txt", "w") as f:
    f.write("Median Standard Deviation of any one star: " + str(numpy.median(sumStd)) +"\n")
    f.write("Standard Error/Uncertainty in Calibration: " +str(errCalib))
    
#print (finalCompUsedFile)
numpy.savetxt("calibCompsUsed.csv", numpy.asarray(finalCompUsedFile), delimiter=",", fmt='%0.8f')    

sys.exit()
