import os
import shutil

if os.path.exists('calibcats'):
    shutil.rmtree('calibcats')
if os.path.exists('periods'):
    shutil.rmtree('periods')
    
if os.path.exists('checkplots'):
    shutil.rmtree('checkplots')
if os.path.exists('eelbs'):
    shutil.rmtree('eelbs')
if os.path.exists('outputcats'):
    shutil.rmtree('outputcats')
if os.path.exists('outputplots'):
    shutil.rmtree('outputplots')
if os.path.exists('trimcats'):
    shutil.rmtree('trimcats')

if os.path.isfile('calibCompsUsed.csv'):
    os.remove('calibCompsUsed.csv')

if os.path.isfile('calibStands.csv'):
    os.remove('calibStands.csv')

if os.path.isfile('compsUsed.csv'):
    os.remove('compsUsed.csv')

if os.path.isfile('screenedComps.csv'):
    os.remove('screenedComps.csv')

if os.path.isfile('starVariability.csv'):
    os.remove('starVariability.csv')

if os.path.isfile('stdComps.csv'):
    os.remove('stdComps.csv')

if os.path.isfile('usedImages.txt'):
    os.remove('usedImages.txt')

if os.path.isfile('LightcurveStats.txt'):
    os.remove('LightcurveStats.txt')

if os.path.isfile('periodEstimates.txt'):
    os.remove('periodEstimates.txt')

if os.path.isfile('calibrationErrors.txt'):
    os.remove('calibrationErrors.txt')

if os.path.isfile('targetstars.csv'):
    os.remove('targetstars.csv')

