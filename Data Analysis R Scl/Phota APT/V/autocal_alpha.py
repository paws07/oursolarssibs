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

Vizier.ROW_LIMIT = -1
max_sep=1.0 * u.arcsec
max_magerr=0.05

fileList=glob.glob('*.cs*')

print (fileList)

for file in fileList:

    calFile = numpy.genfromtxt(file, dtype=float, delimiter=',')
    
    print (file)
    #print (calFile)

    print (file.split('_')[1])
    filterCode=file.split('_')[1]

    # Get Median RA and Dec from file
    print (numpy.median(calFile[:,0]))
    print (numpy.median(calFile[:,1]))
    avgCoord=SkyCoord(ra=(numpy.median(calFile[:,0]))*u.degree, dec=(numpy.median(calFile[:,1]))*u.degree)


    #avgCoord=SkyCoord(ra=(0.4*u.degree), dec=(20*u.degree)) <---- coordinate to test SDSS northern object

    #Get fileCoordinates
    fileCoords=SkyCoord(ra=calFile[:,0]*u.degree, dec=calFile[:,1]*u.degree)

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

            print (sdsssearchResult[:,0])

            raCat=sdssSearchResult[:,0]
            print (raCat)
            decCat=sdssSearchResult[:,1]
            print (decCat)
            if filterCode=='u':
                magCat=sdssSearchResult[:,2]
                print (magCat)
                emagCat=sdssSearchResult[:,3]
                print (emagCat)
            if filterCode=='g':
                magCat=sdssSearchResult[:,4]
                print (magCat)
                emagCat=sdssSearchResult[:,5]
                print (emagCat)
            if filterCode=='r':
                magCat=sdssSearchResult[:,6]
                print (magCat)
                emagCat=sdssSearchResult[:,7]
                print (emagCat)
            if filterCode=='i':
                magCat=sdssSearchResult[:,8]
                print (magCat)
                emagCat=sdssSearchResult[:,9]
                print (emagCat)
            if filterCode=='z':
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

            print (sdsssearchResult[:,0])

            raCat=sdssSearchResult[:,0]
            print (raCat)
            decCat=sdssSearchResult[:,1]
            print (decCat)
            if filterCode=='u':
                magCat=sdssSearchResult[:,2]
                print (magCat)
                emagCat=sdssSearchResult[:,3]
                print (emagCat)
            if filterCode=='g':
                magCat=sdssSearchResult[:,4]
                print (magCat)
                emagCat=sdssSearchResult[:,5]
                print (emagCat)
            if filterCode=='r':
                magCat=sdssSearchResult[:,6]
                print (magCat)
                emagCat=sdssSearchResult[:,7]
                print (emagCat)
            if filterCode=='i':
                magCat=sdssSearchResult[:,8]
                print (magCat)
                emagCat=sdssSearchResult[:,9]
                print (emagCat)
            if filterCode=='z':
                magCat=sdssSearchResult[:,10]
                print (magCat)
                emagCat=sdssSearchResult[:,11]
                print (emagCat)


    # Match to fileCoords
    catCoords=SkyCoord(ra=raCat*u.degree, dec=decCat*u.degree)
    idx,d2d,d3d=fileCoords.match_to_catalog_sky(catCoords)
    print(idx)
##        sep_constraint= d2d < max_sep
##        c_matches=fileCoords[sep_constraint] # matches in original catalogue
##        print (c_matches[0])
##        catalog_matches=fileCoords[idx[sep_constraint]] # matches in catalogue
##        print (catalog_matches[0])

    print (d2d*u.arcsecond)


    magholder=[]
    calibratorz=[]
    for q in range(len(calFile[:,0])):
        #print (calFile[q,0])
        #print (idx[q])
        #print (d2d[q]*u.arcsecond)
        magTempErr = 1.0857 * (calFile[q,5]/calFile[q,4])
        magTemp=-2.5*math.log10(calFile[q,4])
        #enter in appropriate stars into the calibratorz array
        if d2d[q]*u.arcsecond < max_sep*u.arcsecond:
##                if not numpy.isnan(magCat[idx[q]]):
##                    if not numpy.isnan(emagCat[idx[q]]):
                    
            
            magholder.append([magTemp,magTempErr,magCat[idx[q]],emagCat[idx[q]]])
            calibratorz.append([magTemp,magTempErr,magCat[idx[q]],emagCat[idx[q]]])
        else:
            magholder.append([magTemp,magTempErr,numpy.nan,numpy.nan])

    print (calibratorz)

    calibratorz=numpy.asarray(calibratorz)
    magholder=numpy.asarray(magholder)

    # get rid of Nans
    calibRejects=[]
    for s in range(len(calibratorz[:,0])):
        if (numpy.isnan(calibratorz[s,2])):
            calibRejects.append(s)
    if calibRejects!=[]:
        calibratorz=numpy.delete(calibratorz, calibRejects, axis=0)

    # get rid of high errors
    calibRejects=[]
    for s in range(len(calibratorz[:,0])):
        if (calibratorz[s,3] > max_magerr):
            calibRejects.append(s)
    if calibRejects!=[]:
        calibratorz=numpy.delete(calibratorz, calibRejects, axis=0)

    # Estimate best zeropoint
    while True:    
        
        tempZP=calibratorz[:,2]-calibratorz[:,0]
        #print (tempZP)
        #print (numpy.nanmedian(tempZP))
        #print (numpy.nanstd(tempZP))
        tempmedian=numpy.nanmedian(tempZP)
        tempnanstd=numpy.nanstd(tempZP)
        calibRejects=[]
        for s in range(len(calibratorz[:,0])):
            if (calibratorz[s,2]-calibratorz[s,0]) > (tempmedian + (tempnanstd*2)):
                calibRejects.append(s)
            if (calibratorz[s,2]-calibratorz[s,0]) < (tempmedian - (tempnanstd*2)):
                calibRejects.append(s)
        if calibRejects==[]:
            print (tempZP)
            print (numpy.nanmedian(tempZP))
            print (numpy.nanstd(tempZP))
            break
        else:
            calibratorz=numpy.delete(calibratorz, calibRejects, axis=0)

        
        
            
    #magErrTotal = pow( pow(magErrVar,2) + pow(magErrEns,2),0.5)
    
    numpy.savetxt(file.replace('.csv','COMPARISONS.csv'), calibratorz, delimiter=",", fmt='%0.8f')

    #replace counts by mag
    for t in range(len(magholder[:,0])):
        calFile[t,4]=magholder[t,0]+tempmedian
        calFile[t,5]=magholder[t,1]
            

    numpy.savetxt(file.replace('.csv','CALIBED.csv'), calFile, delimiter=",", fmt='%0.8f')
        


sys.exit()
