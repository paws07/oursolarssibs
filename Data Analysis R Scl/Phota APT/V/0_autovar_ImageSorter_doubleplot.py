import numpy
from astropy.io import fits
import glob
import sys
import os
import platform
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import shutil

photFileProcess=1 # set this to 1 to also reject photfiles

# Get list of phot files
parentPath = os.getcwd()
if "Windows" in platform.platform():
    imagesPath = os.path.join(parentPath,"images\\*.f*")
    inputPath = os.path.join(parentPath,"inputs\\*.*")
    rejectPath = os.path.join(parentPath,"rejects\\")
    imagesdonePath = os.path.join(parentPath,"imagesdone\\")
    inputsdonePath = os.path.join(parentPath,"inputsdone\\")

if "Mac" or "Darwin" in platform.platform():
    imagesPath = os.path.join(parentPath,"images//*.f*")
    inputPath = os.path.join(parentPath,"inputs//*.*")
    rejectPath = os.path.join(parentPath,"rejects//")
    imagesdonePath = os.path.join(parentPath,"imagesdone//")
    inputsdonePath = os.path.join(parentPath,"inputsdone//")

print (platform.platform())

if not os.path.exists(rejectPath):
    os.makedirs(rejectPath)
if not os.path.exists(imagesdonePath):
    os.makedirs(imagesdonePath)
if not os.path.exists(inputsdonePath):
    os.makedirs(inputsdonePath)

fileList = glob.glob(imagesPath)
# fileList = sorted(glob.glob(imagesPath), reverse=True)

for file in fileList:
    fitsFile=fits.open(file)
    image_data=fitsFile[0].data
    fitsFile.close()
    image_data=numpy.arcsinh(image_data)
    plt.subplot(1,2,1)
    plt.imshow(image_data, cmap='Greys', norm=LogNorm())
    plt.colorbar()
    plt.subplot(1,2,2)
    plt.imshow(image_data, cmap='gray', norm=LogNorm())
    plt.colorbar()
    mng = plt.get_current_fig_manager()
    # mng.window.state('zoomed')
    # this command could be
    # mng.frame.Maximize(True)
    # or
    mng.window.showMaximized()
    # on a Mac... let me know if it is!
    plt.show()
    while True:
        print ("Keep that image?:")
        keepQ = input()
        if keepQ == 'y' or keepQ =='Y':
            if "Windows" in platform.platform():
                print (' ')
                print (file.split('\\')[-1])
                print (' ')
                fileName = file.split('\\')[-1]
                print (os.path.join(parentPath,'imagesdone\\'+fileName))
                shutil.move(file, os.path.join(parentPath,"imagesdone\\"+fileName))
                if photFileProcess ==1 :
                    photFileName=fileName.split('.')[0]+'.*'
                    photFileList = glob.glob(os.path.join(parentPath,"inputs\\"+photFileName))
                    print (photFileList)
                    for movephotfile in photFileList:
                        movephotfile=movephotfile.split('\\')[-1]
                        print (movephotfile)
                        shutil.move(os.path.join(parentPath,"inputs\\"+movephotfile), os.path.join(parentPath,"inputsdone\\"+movephotfile))
            if ("Mac" in platform.platform() )or ("Darwin" in platform.platform()):
                print (' ')
                print (file.split('/')[-1])
                print (' ')
                fileName = file.split('/')[-1]
                print (os.path.join(parentPath,"imagesdone/"+fileName))
                shutil.move(file, os.path.join(parentPath,"imagesdone/"+fileName))
                if photFileProcess ==1 :
                    photFileName=fileName.split('.')[0]+'.*'
                    photFileList = glob.glob(os.path.join(parentPath,"inputs/"+photFileName))
                    print (photFileList)
                    for movephotfile in photFileList:
                        movephotfile=movephotfile.split('/')[-1]
                        print (movephotfile)
                        shutil.move(os.path.join(parentPath,"inputs/"+movephotfile), os.path.join(parentPath,"inputsdone/"+movephotfile))
            break
        if keepQ == 'n' or keepQ =='N':
            print ("reject!")
            if "Windows" in platform.platform():
                print (' ')
                print (file.split('\\')[-1])
                print (' ')
                fileName = file.split('\\')[-1]
                print (os.path.join(parentPath,'rejects\\'+fileName))
                shutil.move(file, os.path.join(parentPath,"rejects\\"+fileName))
                if photFileProcess ==1 :
                    photFileName=fileName.split('.')[0]+'.*'
                    photFileList = glob.glob(os.path.join(parentPath,"inputs\\"+photFileName))
                    print (photFileList)
                    for movephotfile in photFileList:
                        movephotfile=movephotfile.split('\\')[-1]
                        print (movephotfile)
                        shutil.move(os.path.join(parentPath,"inputs\\"+movephotfile), os.path.join(parentPath,"rejects\\"+movephotfile))
            if ("Mac" in platform.platform() )or ("Darwin" in platform.platform()):
                print (' ')
                print (file.split('/')[-1])
                print (' ')
                fileName = file.split('/')[-1]
                print (os.path.join(parentPath,"rejects/"+fileName))
                shutil.move(file, os.path.join(parentPath,"rejects/"+fileName))
                if photFileProcess ==1 :
                    photFileName=fileName.split('.')[0]+'.*'
                    photFileList = glob.glob(os.path.join(parentPath,"inputs/"+photFileName))
                    print (photFileList)
                    for movephotfile in photFileList:
                        movephotfile=movephotfile.split('/')[-1]
                        print (movephotfile)
                        shutil.move(os.path.join(parentPath,"inputs/"+movephotfile), os.path.join(parentPath,"rejects/"+movephotfile))



                
            break

    
    
    
sys.exit()



