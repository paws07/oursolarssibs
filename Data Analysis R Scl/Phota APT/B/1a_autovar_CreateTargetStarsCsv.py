import numpy

# Generate a blank targetstars.csv file
targetStars=[(0,0,0,0)]
numpy.savetxt("targetstars.csv", targetStars, delimiter=",", fmt='%0.8f')
