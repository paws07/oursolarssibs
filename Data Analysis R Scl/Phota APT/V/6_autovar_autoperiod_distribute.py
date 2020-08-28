#from astropy.coordinates import SkyCoord
#from astropy import units as u
#from astropy.time import Time
import numpy
#import bokeh.plotting as bk
#from bokeh.models import Range1d
#from bokeh.layouts import row
#from bokeh.plotting import figure
#import json
import sys
import os
import platform
import glob
import matplotlib.pyplot as plt

numBins = 10
minperiod=0.2
maxperiod=1.2
periodsteps=10000


# Potentially incorporate -> http://www.astroml.org/gatspy/

#########################################

def sortByPhase (phases, fluxes):
  phaseIndices = numpy.asarray(phases).argsort()
  sortedPhases = []
  sortedFluxes = []
  for i in range(0,len(phases)):
    sortedPhases.append(phases[phaseIndices[i]])
    sortedFluxes.append(fluxes[phaseIndices[i]])

  return (sortedPhases, sortedFluxes)

#########################################

def normalize(fluxes):
  
  normalizedFluxes = []
  
  for flux in fluxes:
    normalizedFlux = (flux - min(fluxes)) / (max(fluxes) - min(fluxes))
    normalizedFluxes.append(normalizedFlux)
  
  return(normalizedFluxes)

######################################### 
  
def getPhases(julian_dates, fluxes, period):

  phases = []
    
  for n in range(0, len(julian_dates)):
    phases.append((julian_dates[n] / (period)) % 1)
    
  (sortedPhases, sortedFluxes) = sortByPhase(phases, fluxes)
    
  return(sortedPhases, sortedFluxes)

#########################################

# Find the value in array2 corresponding to the minimum value in array1.

def find_minimum(array1, array2):
  value_to_return = 0.0
  
  minimum = min(array1)
   
  for i in range(0, len(array1)):
    if (array1[i] == minimum):
      value_to_return = array2[i]

  return (value_to_return, minimum)

#########################################

def sum_distances (sortedPhases, sortedNormalizedFluxes):
  
  distanceSum = 0.0
  
  for i in range(0, (len(sortedPhases) - 1)):
    fluxdiff = sortedNormalizedFluxes[i + 1] - sortedNormalizedFluxes[i]
    phasediff = sortedPhases[i + 1] - sortedPhases[i]
    
    distanceSum = distanceSum + (((fluxdiff ** 2) + (phasediff ** 2)) ** 0.5)

  return(distanceSum)

#########################################

def sum_stdevs (sortedPhases, sortedNormalizedFluxes, numBins):
   
  stdevSum = 0.0
  
  for i in range (0, numBins):
    fluxes_inrange = []
    minIndex = (float(i) / float(numBins)) * float(len(sortedPhases))
    maxIndex = (float(i + 1) / float(numBins)) * float(len(sortedPhases))

    for j in range (0, len(sortedPhases)):
      if (j >= minIndex and j < maxIndex):
        fluxes_inrange.append(sortedNormalizedFluxes[j])
    
    stdev_of_bin_i = numpy.std(fluxes_inrange)
    stdevSum = stdevSum + stdev_of_bin_i
    
  return(stdevSum)
    
#########################################

def phase_dispersion_minimization ():  

  periodguess_array = []
  
  distance_results = []
  stdev_results = []
                        
  (julian_dates, fluxes) = (varData[:,0],varData[:,1])
  normalizedFluxes = normalize(fluxes)

  for r in range(periodsteps): 
    periodguess = minperiod + (r * ((maxperiod-minperiod)/periodsteps))
    (sortedPhases, sortedNormalizedFluxes) = getPhases(julian_dates, normalizedFluxes, periodguess)
    
    distance_sum = sum_distances(sortedPhases, sortedNormalizedFluxes)
    stdev_sum = sum_stdevs(sortedPhases, sortedNormalizedFluxes, numBins)

    periodguess_array.append(periodguess)
    distance_results.append(distance_sum)
    stdev_results.append(stdev_sum)

  periodTrialMatrix=[]
  for r in range(periodsteps):
    periodTrialMatrix.append([periodguess_array[r],distance_results[r],stdev_results[r]])
  periodTrialMatrix=numpy.asarray(periodTrialMatrix)
  numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"Trials.csv"), periodTrialMatrix, delimiter=",", fmt='%0.8f')
        
  (distance_minperiod, distance_min) = find_minimum(distance_results, periodguess_array)
  (stdev_minperiod, stdev_min) = find_minimum(stdev_results, periodguess_array)
  
  pdm = {}
  pdm["periodguess_array"] = periodguess_array
  pdm["distance_results"] = distance_results
  pdm["distance_minperiod"] = distance_minperiod
  pdm["stdev_results"] = stdev_results
  pdm["stdev_minperiod"] = stdev_minperiod

  # Estimating the error
  # stdev method
  #print (numpy.min(stdev_results))
  #print (numpy.max(stdev_results))
  # Get deviation to the left
  totalRange=numpy.max(stdev_results) - numpy.min(stdev_results)
  for q in range(len(periodguess_array)):
    if periodguess_array[q]==pdm["stdev_minperiod"]:
      beginIndex=q
      beginValue=stdev_results[q]
  #print (beginIndex)
  #print (beginValue)
  currentperiod=stdev_minperiod
  stepper=0
  thresholdvalue=beginValue+(0.5*totalRange)
  while True:
    #print (beginIndex-stepper)
    #print (stdev_results[beginIndex-stepper])
    if stdev_results[beginIndex-stepper] > thresholdvalue:
      #print ("LEFTHAND PERIOD!")
      #print (periodguess_array[beginIndex-stepper])
      lefthandP=periodguess_array[beginIndex-stepper]
      #print (distance_results)
      break
    stepper=stepper+1

  stepper=0
  thresholdvalue=beginValue+(0.5*totalRange)
  #print (beginIndex)
  #print (periodsteps)
  
  while True:
    #print (beginIndex+stepper)
    #print (stdev_results[beginIndex+stepper])
    if beginIndex+stepper+1 == periodsteps:
      righthandP=periodguess_array[beginIndex+stepper]
      print ("Warning: Peak period for stdev method too close to top of range")
      break
    if stdev_results[beginIndex+stepper] > thresholdvalue:
      #print ("RIGHTHAND PERIOD!")
      #print (periodguess_array[beginIndex+stepper])
      righthandP=periodguess_array[beginIndex+stepper]
      #print (distance_results)
      break
    stepper=stepper+1


  #print ("Stdev method error: " + str((righthandP - lefthandP)/2))
  pdm["stdev_error"] = (righthandP - lefthandP)/2


  # Estimating the error
  # stdev method
  #print (numpy.min(stdev_results))
  #print (numpy.max(stdev_results))
  # Get deviation to the left
  totalRange=numpy.max(distance_results) - numpy.min(distance_results)
  for q in range(len(periodguess_array)):
    if periodguess_array[q]==pdm["distance_minperiod"]:
      beginIndex=q
      beginValue=distance_results[q]
  #print (beginIndex)
  #print (beginValue)
  currentperiod=distance_minperiod
  stepper=0
  thresholdvalue=beginValue+(0.5*totalRange)
  while True:
    #print (beginIndex-stepper)
    #print (stdev_results[beginIndex-stepper])
    if distance_results[beginIndex-stepper] > thresholdvalue:
      #print ("LEFTHAND PERIOD!")
      #print (periodguess_array[beginIndex-stepper])
      lefthandP=periodguess_array[beginIndex-stepper]
      #print (distance_results)
      break
    stepper=stepper+1

  stepper=0
  thresholdvalue=beginValue+(0.5*totalRange)
  while True:
    if beginIndex+stepper+1 == periodsteps:
      righthandP=periodguess_array[beginIndex+stepper]
      print ("Warning: Peak period for distance method too close to top of range")
      break
    #print (beginIndex+stepper)
    #print (stdev_results[beginIndex+stepper])
    if distance_results[beginIndex+stepper] > thresholdvalue:
      #print ("RIGHTHAND PERIOD!")
      #print (periodguess_array[beginIndex+stepper])
      righthandP=periodguess_array[beginIndex+stepper]
      #print (distance_results)
      break
    stepper=stepper+1

  #print ("Distance method error: " + str((righthandP - lefthandP)/2))
  pdm["distance_error"] = (righthandP - lefthandP)/2
  
  return (pdm)

#########################################



trialRange=[minperiod, maxperiod]

# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs\\*.p*")
    outputPath = os.path.join(parentPath,"outputplots\\")
    outcatPath = os.path.join(parentPath,"outputcats\\")
    checkPath = os.path.join(parentPath,"checkplots\\")
    #doerPath = os.path.join(parentPath,"outputcats\\doer*.csv")
    periodPath = os.path.join(parentPath,"periods\\")

if "Mac" or "Darwin" in platform.platform():
    #inputPath = os.path.join(parentPath,"inputs//*.p*")
    outputPath = os.path.join(parentPath,"outputplots//")
    outcatPath = os.path.join(parentPath,"outputcats//")
    checkPath = os.path.join(parentPath,"checkplots//")
    periodPath = os.path.join(parentPath,"periods//")
    #doerPath = os.path.join(parentPath,"outputcats//doer*.csv")

# Load in list of used files to get filterCode
fileList=[]
with open("usedImages.txt", "r") as f:
  for line in f:
    fileList.append(line.strip())

# Detect Filter Set being used
filterCode = (fileList[0].replace(parentPath,"").replace('inputs',"").replace('//',"").replace('\\',"").split('_')[1])
print ("Filter Set: " + filterCode)



#create directory structure
if not os.path.exists(outcatPath):
    os.makedirs(outcatPath)

if not os.path.exists(periodPath):
    os.makedirs(periodPath)

fileList = glob.glob(outcatPath+'*_diffExcel.csv')
with open("periodEstimates.txt", "w") as f:
  f.write("Period Estimates \n\n")
  
# Load in the files
for file in fileList:
  print (file)
  variableName=str(file).replace(str(outcatPath).replace('//','\\'),"").split('_')[0]
  #print (str(outcatPath).replace('//',''))
  print (variableName)
  varData=numpy.genfromtxt(file, dtype=float, delimiter=',')
  if os.path.isfile(file.replace('diff','calib')):
    calibData=numpy.genfromtxt(file.replace('diff','calib'), dtype=float, delimiter=',')

  #print (minDate)

  pdm_results = {}

  pdm=phase_dispersion_minimization()

  plt.figure(figsize=(15, 5))

  print("Distance Method Estimate (days): " + str(pdm["distance_minperiod"]))
  #print(pdm["distance_minperiod"] )
  print ("Distance method error: " + str(pdm["distance_error"]))
  phaseTest=(varData[:,0] / (pdm["distance_minperiod"])) % 1
  with open("periodEstimates.txt", "a+") as f:
    f.write("Variable : "+str(variableName) +"\n")
    f.write("Distance Method Estimate (days): " + str(pdm["distance_minperiod"])+"\n")
    f.write("Distance method error: " + str(pdm["distance_error"])+"\n")
  
  
  plt.plot(pdm["periodguess_array"], pdm["distance_results"])
  plt.gca().invert_yaxis()
  plt.title("Range {0} d  Steps: {1}".format(trialRange, periodsteps))
  plt.xlabel(r"Trial Period")
  plt.ylabel(r"Likelihood of Period")
  plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"StringLikelihoodPlot.png"))
  plt.clf()

  plt.plot(phaseTest, varData[:,1], 'bo', linestyle='None')
  plt.plot(phaseTest+1, varData[:,1], 'ro', linestyle='None')
  plt.errorbar(phaseTest, varData[:,1], yerr=varData[:,2], linestyle='None')
  plt.errorbar(phaseTest+1, varData[:,1], yerr=varData[:,2], linestyle='None')
  plt.gca().invert_yaxis()
  plt.title("Period: {0} d  Steps: {1}".format(pdm["distance_minperiod"], periodsteps))
  plt.xlabel(r"Phase ($\phi$)")
  plt.ylabel(r"Differential " + str(filterCode) + " Magnitude")
  plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"StringTestPeriodPlot.png"))
  plt.clf()

  if os.path.isfile(file.replace('diff','calib')):
    phaseTestCalib=(calibData[:,0] / (pdm["distance_minperiod"])) % 1
    plt.plot(phaseTestCalib, calibData[:,1], 'bo', linestyle='None')
    plt.plot(phaseTestCalib+1, calibData[:,1], 'ro', linestyle='None')
    plt.errorbar(phaseTestCalib, calibData[:,1], yerr=varData[:,2], linestyle='None')
    plt.errorbar(phaseTestCalib+1, calibData[:,1], yerr=varData[:,2], linestyle='None')
    plt.gca().invert_yaxis()
    plt.title("Period: {0} d  Steps: {1}".format(pdm["distance_minperiod"], periodsteps))
    plt.xlabel(r"Phase ($\phi$)")
    plt.ylabel(r"Calibrated " + str(filterCode) + " Magnitude")
    plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"StringTestPeriodPlot_Calibrated.png"))
    plt.clf()
  
  tempPeriodCatOut=[]
  for g in range(len(phaseTest)):
    tempPeriodCatOut.append([phaseTest[g],varData[g,1]])
  tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
  numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"StringTrial.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')

  if os.path.isfile(file.replace('diff','calib')): 
    tempPeriodCatOut=[]
    for g in range(len(calibData[:,0])):
      tempPeriodCatOut.append([(calibData[g,0]/(pdm["distance_minperiod"]) % 1), calibData[g,1], calibData[g,2]])
    tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
    numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"String_PhasedCalibMags.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')

  tempPeriodCatOut=[]
  for g in range(len(varData[:,0])):
    tempPeriodCatOut.append([(varData[g,0]/(pdm["distance_minperiod"]) % 1), varData[g,1], varData[g,2]])
  tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
  numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"String_PhasedDiffMags.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')
  
  
  print("PDM Method Estimate (days): "+ str(pdm["stdev_minperiod"]))
  #print(pdm["stdev_minperiod"])
  phaseTest=(varData[:,0] / (pdm["stdev_minperiod"])) % 1
  print ("PDM method error: " + str(pdm["stdev_error"]))

  with open("periodEstimates.txt", "a+") as f:
    f.write("PDM Method Estimate (days): "+ str(pdm["stdev_minperiod"])+"\n")
    f.write("PDM method error: " + str(pdm["stdev_error"])+"\n\n")

  
  plt.plot(pdm["periodguess_array"], pdm["stdev_results"])
  plt.gca().invert_yaxis()
  plt.title("Range {0} d  Steps: {1}".format(trialRange, periodsteps))
  plt.xlabel(r"Trial Period")
  plt.ylabel(r"Likelihood of Period")
  plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"PDMLikelihoodPlot.png"))
 
  plt.clf()


  plt.plot(phaseTest, varData[:,1], 'bo', linestyle='None')
  plt.plot(phaseTest+1, varData[:,1], 'ro', linestyle='None')
  plt.errorbar(phaseTest, varData[:,1], yerr=varData[:,2], linestyle='None')
  plt.errorbar(phaseTest+1, varData[:,1], yerr=varData[:,2], linestyle='None')
  plt.gca().invert_yaxis()
  plt.title("Period: {0} d  Steps: {1}".format(pdm["stdev_minperiod"], periodsteps))
  plt.xlabel(r"Phase ($\phi$)")
  plt.ylabel(r"Differential " + str(filterCode) + " Magnitude")
  plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"PDMTestPeriodPlot.png"))
  plt.clf()

  if os.path.isfile(file.replace('diff','calib')):
    phaseTestCalib=(calibData[:,0] / (pdm["stdev_minperiod"])) % 1
    plt.plot(phaseTestCalib, calibData[:,1], 'bo', linestyle='None')
    plt.plot(phaseTestCalib+1, calibData[:,1], 'ro', linestyle='None')
    plt.errorbar(phaseTestCalib, calibData[:,1], yerr=varData[:,2], linestyle='None')
    plt.errorbar(phaseTestCalib+1, calibData[:,1], yerr=varData[:,2], linestyle='None')
    plt.gca().invert_yaxis()
    plt.title("Period: {0} d  Steps: {1}".format(pdm["stdev_minperiod"], periodsteps))
    plt.xlabel(r"Phase ($\phi$)")
    plt.ylabel(r"Calibrated " + str(filterCode) + " Magnitude")
    plt.savefig(os.path.join(periodPath,str(variableName)+'_'+"PDMTestPeriodPlot_Calibrated.png"))
    plt.clf()

  tempPeriodCatOut=[]
  for g in range(len(phaseTest)):
    tempPeriodCatOut.append([phaseTest[g],varData[g,1]])
  tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
  numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"PDMTrial.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')

  if os.path.isfile(file.replace('diff','calib')):
    tempPeriodCatOut=[]
    for g in range(len(calibData[:,0])):
      tempPeriodCatOut.append([(calibData[g,0]/(pdm["stdev_minperiod"])) % 1, calibData[g,1], calibData[g,2]])
    tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
    numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"PDM_PhasedCalibMags.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')

  tempPeriodCatOut=[]
  for g in range(len(varData[:,0])):
    tempPeriodCatOut.append([(varData[g,0]/(pdm["stdev_minperiod"])) % 1, varData[g,1], varData[g,2]])
  tempPeriodCatOut=numpy.asarray(tempPeriodCatOut)
  numpy.savetxt(os.path.join(periodPath,str(variableName)+'_'+"PDM_PhaseddiffMags.csv"), tempPeriodCatOut, delimiter=",", fmt='%0.8f')



sys.exit()


      
