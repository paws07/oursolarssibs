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

errorReject=0.5 # reject measurements with instrumental errors larger than this (this is not total error, just the estimated error in the single measurement of the variable)
acceptDistance=5.0 # Furtherest distance in arcseconds for matches

# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")
    doerPath = os.path.join(parentPath,"outputcats\\doer*.csv")

if "Mac" or "Darwin" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
    doerPath = os.path.join(parentPath,"outputcats//doer*.csv")

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

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)

# LOAD Phot FILES INTO LIST
photFileArray=[]
for file in fileList:
    loadPhot=numpy.genfromtxt(file, dtype=float, delimiter=',')
    if loadPhot.shape[1] > 6:
        loadPhot=numpy.delete(loadPhot,6,1)
        loadPhot=numpy.delete(loadPhot,6,1)
    photFileArray.append(loadPhot)
photFileArray=numpy.asarray(photFileArray)

#Load in the targets
targetFile = numpy.genfromtxt('targetstars.csv', dtype=float, delimiter=',')
# Remove any nan rows from targetFile
targetRejecter=[]
if not (targetFile.shape[0] == 4 and targetFile.size ==4):
    for z in range(targetFile.shape[0]):
      if numpy.isnan(targetFile[z][0]):
        targetRejecter.append(z)
    targetFile=numpy.delete(targetFile, targetRejecter, axis=0)

# Detect whether Differential or Calibrated photometry
exists=os.path.isfile('calibCompsUsed.csv')
if exists:
    print ("Calibrated")
    compFile=numpy.genfromtxt('calibCompsUsed.csv', dtype=float, delimiter=',')
    calibFlag=1    
else:
    print ("Differential")
    compFile=numpy.genfromtxt('compsUsed.csv', dtype=float, delimiter=',')
    calibFlag=0

# Get total counts for each file
fileCount=[]
compArray=[]
allCountsArray=[]
for imgs in range(photFileArray.shape[0]):
    allCounts=0.0
    allCountsErr=0.0
    photFile = photFileArray[imgs]
    fileRaDec = SkyCoord(ra=photFile[:,0]*u.degree, dec=photFile[:,1]*u.degree)
    print ("Calculating total Comparison counts for " + str(fileList[imgs]))
    #print (compFile.shape[0])


    if (compFile.shape[0]== 5 and compFile.size ==5) or (compFile.shape[0]== 3 and compFile.size ==3):
        loopLength=1
    else:
        loopLength=compFile.shape[0]
    #print (compFile.size)
                #sys.exit()
    for j in range(loopLength):        
        if compFile.size == 2 or (compFile.shape[0]== 3 and compFile.size ==3) or (compFile.shape[0]== 5 and compFile.size ==5):
##    
##    for j in range(compFile.shape[0]):        
##        if compFile.size == 2 or compFile.shape[0] == 5:
            matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
        else:
            matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
        idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)
        allCounts=numpy.add(allCounts,photFile[idx][4])
        allCountsErr=numpy.add(allCountsErr,photFile[idx][5])

    allCountsArray.append([allCounts,allCountsErr])

print (allCountsArray)

allcountscount=0

if len(targetFile)== 4:
    loopLength=1
else:
    loopLength=targetFile.shape[0]
# For each variable calculate all the things
for q in range(loopLength):
    starErrorRejCount=0
    starDistanceRejCount=0
    print ("***********************************************************************")
    print ("Processing Variable " +str(q+1))
    print ("RA")
    if int(len(targetFile)) == 4:
        print (targetFile[0])
    else:
        print (targetFile[q][0])
    print ("DEC")
    if int(len(targetFile)) == 4:
        print (targetFile[1])
    else:
        print (targetFile[q][1])
    if int(len(targetFile)) == 4:
        varCoord = SkyCoord(targetFile[0],(targetFile[1]), frame='icrs', unit=u.deg) # Need to remove target stars from consideration
    else:
        varCoord = SkyCoord(targetFile[q][0],(targetFile[q][1]), frame='icrs', unit=u.deg) # Need to remove target stars from consideration
     
    # Grabbing variable rows
    print ("Extracting and Measuring Differential Magnitude in each Photometry File")
    outputPhot=[] # new
    compArray=[]
    compList=[]
    allcountscount=0
    for imgs in range(photFileArray.shape[0]):
        compList=[]
        fileRaDec = SkyCoord(ra=photFileArray[imgs][:,0]*u.degree, dec=photFileArray[imgs][:,1]*u.degree)
        idx, d2d, _ = varCoord.match_to_catalog_sky(fileRaDec)
        starRejected=0
        if (numpy.less(d2d.arcsecond, acceptDistance)):
            magErrVar = 1.0857 * (photFileArray[imgs][idx][5]/photFileArray[imgs][idx][4])
            if magErrVar < errorReject:

                magErrEns = 1.0857 * (allCountsErr/allCounts)
                magErrTotal = pow( pow(magErrVar,2) + pow(magErrEns,2),0.5)

                #templist is a temporary holder of the resulting file. 
                tempList=photFileArray[imgs][idx,:]                    
                googFile = (fileList[imgs].replace(parentPath,"").replace('inputs',"").replace('//',""))
                tempList=numpy.append(tempList, float(googFile.split("_")[2].replace("d",".")))
                tempList=numpy.append(tempList, float(googFile.split("_")[4].replace("a",".")))
                tempList=numpy.append(tempList, allCountsArray[allcountscount][0])
                tempList=numpy.append(tempList, allCountsArray[allcountscount][1])

                #Differential Magnitude
                #tempList=numpy.append(tempList,-2.5 * numpy.log10(photFileArray[imgs][idx][4]/allCountsArray[allcountscount][0]))
                tempList=numpy.append(tempList, 2.5 * numpy.log10(allCountsArray[allcountscount][0]/photFileArray[imgs][idx][4]))
                tempList=numpy.append(tempList, magErrTotal)
                tempList=numpy.append(tempList, photFileArray[imgs][idx][4])
                tempList=numpy.append(tempList, photFileArray[imgs][idx][5])

##                for j in range(compFile.shape[0]):        
##                    if compFile.size == 2 or compFile.shape[0] == 5:
##                        matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
##                    else:
##                        matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
##                    idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)        
##                    tempList=numpy.append(tempList, photFileArray[imgs][idx][4])

                if (compFile.shape[0]== 5 and compFile.size ==5) or (compFile.shape[0]== 3 and compFile.size ==3):
                    loopLength=1
                else:
                    loopLength=compFile.shape[0]
                #print (compFile.size)
                #sys.exit()
                for j in range(loopLength):        
                    if compFile.size == 2 or (compFile.shape[0]== 3 and compFile.size ==3) or (compFile.shape[0]== 5 and compFile.size ==5):
                        matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
                    else:
                        matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
                    idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)        
                    tempList=numpy.append(tempList, photFileArray[imgs][idx][4])  
                    
                outputPhot.append(tempList)
                fileCount.append(allCounts)
                allcountscount=allcountscount+1

            else:
                starErrorRejCount=starErrorRejCount+1
                starRejected=1
        else:
            starDistanceRejCount=starDistanceRejCount+1
            starRejected=1

        if ( starRejected == 1):
            
                #templist is a temporary holder of the resulting file. 
                tempList=photFileArray[imgs][idx,:]                    
                googFile = (fileList[imgs].replace(parentPath,"").replace('inputs',"").replace('//',""))
                tempList=numpy.append(tempList, float(googFile.split("_")[2].replace("d",".")))
                tempList=numpy.append(tempList, float(googFile.split("_")[4].replace("a",".")))
                tempList=numpy.append(tempList, allCountsArray[allcountscount][0])
                tempList=numpy.append(tempList, allCountsArray[allcountscount][1])

                #Differential Magnitude
                tempList=numpy.append(tempList,numpy.nan)
                tempList=numpy.append(tempList,numpy.nan)
                tempList=numpy.append(tempList, photFileArray[imgs][idx][4])
                tempList=numpy.append(tempList, photFileArray[imgs][idx][5])
                #print "Extracting individual comps"
##                for j in range(compFile.shape[0]):        
##                    if compFile.size == 2 or compFile.shape[0] == 5:
##                        matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
##                    else:
##                        matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
##                    idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)        
##                    tempList=numpy.append(tempList, photFileArray[imgs][idx][4])

                if (compFile.shape[0]== 5 and compFile.size ==5) or (compFile.shape[0]== 3 and compFile.size ==3):
                    loopLength=1
                else:
                    loopLength=compFile.shape[0]
                #print (compFile.shape[0])
                #sys.exit()
                for j in range(loopLength):        
                    if compFile.size == 2 or (compFile.shape[0]== 3 and compFile.size ==3) or (compFile.shape[0]== 5 and compFile.size ==5):
                        matchCoord=SkyCoord(ra=compFile[0]*u.degree, dec=compFile[1]*u.degree)
                    else:
                        matchCoord=SkyCoord(ra=compFile[j][0]*u.degree, dec=compFile[j][1]*u.degree)        
                    idx, d2d, d3d = matchCoord.match_to_catalog_sky(fileRaDec)        
                    tempList=numpy.append(tempList, photFileArray[imgs][idx][4])                    
                outputPhot.append(tempList)            
                fileCount.append(allCounts)
                allcountscount=allcountscount+1

    # Check for dud images
    imageReject=[]
    for j in range(numpy.asarray(outputPhot).shape[0]):
        if numpy.isnan(outputPhot[j][11]):
            imageReject.append(j)
    outputPhot=numpy.delete(outputPhot, imageReject, axis=0)

    ## REMOVE MAJOR OUTLIERS FROM CONSIDERATION
    stdVar=numpy.nanstd(numpy.asarray(outputPhot)[:,10])
    avgVar=numpy.nanmean(numpy.asarray(outputPhot)[:,10])
    starReject=[]
    stdevReject=0
    for j in range(numpy.asarray(outputPhot).shape[0]):
        if outputPhot[j][10] > avgVar+(4*stdVar) or outputPhot[j][10] < avgVar-(4*stdVar) :
            starReject.append(j)
            stdevReject=stdevReject+1

    print ("Rejected Stdev Measurements: " + str(stdevReject))
    print ("Rejected Error Measurements: " + str(starErrorRejCount))
    print ("Rejected Distance Measurements: " + str(starDistanceRejCount))
    print ("Variability of Comparisons")
    print ("Average : " +str(avgVar))
    print ("Stdev   : "+str(stdVar))
    
    outputPhot=numpy.delete(outputPhot, starReject, axis=0)
    if outputPhot.shape[0] > 2:
        numpy.savetxt(os.path.join(outcatPath,"doerPhot_V" +str(q+1) +".csv"), outputPhot, delimiter=",", fmt='%0.8f')

#Make Normal Plots

###BJD vs Airmass
##pylab.cla()
###print compArray
##outplotx=numpy.asarray(compArray)[:,0]
##outploty=numpy.asarray(compArray)[:,1]
###print outplotx
###print outploty
##pylab.xlabel('BJD')
##pylab.ylabel('Airmass')
##pylab.plot(outplotx,outploty,'bo')
###pylab.plot(linex,liney)
##pylab.ylim(max(outploty)+0.1,min(outploty)-0.1,'k-')
##pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
##pylab.grid(True)
##pylab.savefig(os.path.join(checkPath,'BJDairmass.png'))
##pylab.savefig(os.path.join(checkPath,'BJDairmass.eps'))
##
##pylab.cla()
##outplotx=numpy.asarray(compArray)[:,0]
##outploty=numpy.asarray(compArray)[:,2]
###print outplotx
###print outploty
##pylab.xlabel('BJD')
##pylab.ylabel('Ensemble Counts')
##pylab.plot(outplotx,outploty,'bo')
###pylab.plot(linex,liney)
##pylab.ylim(min(outploty)-1000,max(outploty)+1000,'k-')
##pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
##pylab.grid(True)
##pylab.savefig(os.path.join(checkPath,'BJDensCounts.png'))
##pylab.savefig(os.path.join(checkPath,'BJDensCounts.eps'))
##
##pylab.cla()
##outplotx=numpy.asarray(compArray)[:,1]
##outploty=numpy.asarray(compArray)[:,2]
###print outplotx
###print outploty
##pylab.xlabel('Airmass')
##pylab.ylabel('Ensemble Counts')
##pylab.plot(outplotx,outploty,'bo')
###pylab.plot(linex,liney)
##pylab.ylim(min(outploty)-1000,max(outploty)+1000,'k-')
##pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
##pylab.grid(True)
##pylab.savefig(os.path.join(checkPath,'AirmassEnsCounts.png'))
##pylab.savefig(os.path.join(checkPath,'AirmassEnsCounts.eps'))

fileList = glob.glob(doerPath)

for file in fileList:
    
    outputPhot=numpy.genfromtxt(file, delimiter=",", dtype='float')
    r = file.split("_")[-1].replace(".csv","")
    print ("Making Plots and Catalogues for Variable " + str(r))

    # Output Differential peranso file
    outputPeransoCalib=[]
    for i in range(outputPhot.shape[0]):
        outputPeransoCalib.append([outputPhot[i][6],outputPhot[i][10],outputPhot[i][11]])

    numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"diffPeranso.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
    numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"diffExcel.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')
        
    # Output Differential astroImageJ file
    outputPeransoCalib=[]
    for i in range(numpy.asarray(outputPhot).shape[0]):
        outputPeransoCalib.append([outputPhot[i][6]-2450000.0,outputPhot[i][10],outputPhot[i][11]])
    numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"diffAIJ.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
    numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"diffAIJ.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')

    pylab.cla()
    outplotx=numpy.asarray(outputPhot)[:,6]
    outploty=numpy.asarray(outputPhot)[:,10]
    pylab.xlabel('BJD')
    pylab.ylabel('Differential ' +filterCode+' Mag')
    pylab.plot(outplotx,outploty,'bo')
    pylab.ylim(max(outploty)+0.02,min(outploty)-0.02,'k-')
    pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
    pylab.grid(True)
    pylab.savefig(os.path.join(outputPath,str(r)+'_'+'EnsembleVarDiffMag.png'))
    pylab.savefig(os.path.join(outputPath,str(r)+'_'+'EnsembleVarDiffMag.eps'))

    pylab.cla()
    outplotx=numpy.asarray(outputPhot)[:,7]
    outploty=numpy.asarray(outputPhot)[:,10]
    pylab.xlabel('Airmass')
    pylab.ylabel('Differential ' +filterCode+' Mag')
    pylab.plot(outplotx,outploty,'bo')
    pylab.ylim(min(outploty)-0.02,max(outploty)+0.02,'k-')
    pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
    pylab.grid(True)
    pylab.savefig(os.path.join(checkPath,str(r)+'_'+'AirmassEnsVarDiffMag.png'))
    pylab.savefig(os.path.join(checkPath,str(r)+'_'+'AirmassEnsVarDiffMag.eps'))

    pylab.cla()
    outplotx=numpy.asarray(outputPhot)[:,7]
    outploty=numpy.asarray(outputPhot)[:,8]
    pylab.xlabel('Airmass')
    pylab.ylabel('Variable Counts')
    pylab.plot(outplotx,outploty,'bo')
    pylab.ylim(min(outploty)-1000,max(outploty)+1000,'k-')
    pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
    pylab.grid(True)
    pylab.savefig(os.path.join(checkPath,str(r)+'_'+'AirmassVarCounts.png'))
    pylab.savefig(os.path.join(checkPath,str(r)+'_'+'AirmassVarCounts.eps'))

    # Make a calibrated version
    # Need to shift the shape of the curve against the lowest error in the catalogue. 

    if calibFlag == 1:
        calibCompFile=numpy.genfromtxt('calibCompsUsed.csv', dtype=float, delimiter=',')
        print ("Calibrating Photometry")
        # Load in calibrated magnitudes and add them
        #print (compFile.size)
        if compFile.shape[0] == 5 and compFile.size != 25:
            ensembleMag=calibCompFile[3]
        else:
            ensembleMag=calibCompFile[:,3]
        ensMag=0

        if compFile.shape[0] == 5 and compFile.size != 25:
            lenloop=1
        else:
            lenloop=len(calibCompFile[:,3])
        for q in range(lenloop):
            if compFile.shape[0] == 5 and compFile.size != 25:
                ensMag=pow(10,-ensembleMag*0.4)
            else:
                ensMag=ensMag+(pow(10,-ensembleMag[q]*0.4))
        #print (ensMag)
        ensembleMag=-2.5*math.log10(ensMag)
        print ("Ensemble Magnitude: "+str(ensembleMag))
        

        #calculate error
        if compFile.shape[0] == 5 and compFile.size !=25:
            ensembleMagError=calibCompFile[4]
            #ensembleMagError=numpy.average(ensembleMagError)*1/pow(ensembleMagError.size, 0.5)
        else:
            ensembleMagError=calibCompFile[:,4]
            ensembleMagError=numpy.average(ensembleMagError)*1/pow(ensembleMagError.size, 0.5)

        #for file in fileList:
        for i in range(outputPhot.shape[0]):
            outputPhot[i][10]=outputPhot[i][10]+ensembleMag
            #outputPhot[i][11]=pow((pow(outputPhot[i][11],2)+pow(ensembleMagError,2)),0.5)


        pylab.cla()
        outplotx=numpy.asarray(outputPhot)[:,6]
        outploty=numpy.asarray(outputPhot)[:,10]
        pylab.xlabel('BJD')
        pylab.ylabel('Calibrated ' +filterCode+' Mag')
        pylab.plot(outplotx,outploty,'bo')
        pylab.ylim(max(outploty)+0.02,min(outploty)-0.02,'k-')
        pylab.xlim(min(outplotx)-0.01,max(outplotx)+0.01)
        pylab.grid(True)
        pylab.savefig(os.path.join(outputPath,str(r)+'_'+'EnsembleVarCalibMag.png'))
        pylab.savefig(os.path.join(outputPath,str(r)+'_'+'EnsembleVarCalibMag.eps'))

        
        # Output Calibed peranso file
        outputPeransoCalib=[]
        r = file.split("_")[-1].replace(".csv","")
        for i in range(outputPhot.shape[0]):
            outputPeransoCalib.append([outputPhot[i][6],outputPhot[i][10],outputPhot[i][11]])
        numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"calibPeranso.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
        numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"calibExcel.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')
            
        # Output astroImageJ file
        outputPeransoCalib=[]
        for i in range(numpy.asarray(outputPhot).shape[0]):
            outputPeransoCalib.append([outputPhot[i][6]-2450000.0,outputPhot[i][10],outputPhot[i][11]])
            #i=i+1

        numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"calibAIJ.txt"), outputPeransoCalib, delimiter=" ", fmt='%0.8f')
        numpy.savetxt(os.path.join(outcatPath,str(r)+'_'+"calibAIJ.csv"), outputPeransoCalib, delimiter=",", fmt='%0.8f')

sys.exit()   

