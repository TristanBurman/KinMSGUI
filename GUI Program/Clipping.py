import numpy as np
from scipy import interpolate,ndimage
from astropy.io import fits
import matplotlib.pyplot as plt


#take the channels with emission out as a separate cube, and glue the channels with noise together as another cube
def emisChannels(cube, start, stop):
   emis = cube[:,:,start:stop]
   noise = np.concatenate((cube[:,:,:start],cube[:,:,stop:]), axis=2)
   return emis, noise

def innersquare(cube, N):
    '''
    cuts out the inner square (spatial) of the data cube or image
    '''
    start = int(len(cube[1])*(7./16.))  #If really bad, use 7/16 and 9/16
    stop = int(len(cube[1])*(9./16.))
    if N == 3:
        return cube[start:stop, start:stop, :]
    elif N == 2:
        return cube[start:stop, start:stop]
    else:
        print("Enter a correct value for N")
def clip(emiscube, noisecube, cutoff):
    '''
    "clip" the cube: set all pixels with a value lower than a certain times sigma to zero
    '''
    #Since the PB correction has already been applied, measure the rms in the inner square of the spatial axes (PB correction makes the outer edges more noisy, so we don't want to use those if we measure the rms
    innersq = innersquare(noisecube, 3)
    sigma = np.std(innersq[np.isfinite(innersq)])
    emiscube[emiscube<cutoff*sigma]=0
    return emiscube
def smoothclip(cube, emiscube, cutoff, header, start, stop):
    '''
    Apply a Gaussian smoothing to the data cube.
    cube = the raw data cube
    emiscube = the part of the cube that contains the emission line [bit with data you want]
    cutoff = the level at which you want to clip (in multiples of the rms noise of the cube)
    header = header of the data cube
    start = first channel that contains emission
    stop = last channel that contains emission    
    '''
    #copy the datacube
    cube_copy = cube.copy() #a copy of the data cube will be used to do the smoothing. This smoothed cube will then be used as a mask to clip the orignial cube
    # extract relevant information from header (pixel size, beam size)
    res = hdulist[0].header['CDELT2']  # deg/pixel
    bmaj = BeamData['BMAJ'] # degrees
    bmaj = hdulist[0].header['BMAJ']
    beam = bmaj / res  #beam size in pixels, use the major axis
    
    # convert the FWHM of the beam to the std dev in the Gaussian distribution and use 1.5 times the beam size to convolve with
    sigma = 1.5 * beam / np.sqrt(8.*np.log(2.))
    Test = np.ones(len(sigma)) * 4
    #apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best). The mode 'nearest' seems to give the best results.
    cube_smoothed = ndimage.filters.gaussian_filter(cube_copy, (Test,sigma,sigma), order=0, mode='nearest')
    #clip the copied cube, using the smaller cube where the emission is
    smoothemis, smoothnoise = emisChannels(cube_smoothed, start, stop)  #this separates the channels with emission from the ones without emission and glues the latter bits together in one "noise cube"
    clipped = clip(smoothemis, smoothnoise, cutoff)
    #mask the original cube with the smoothed cube
    emiscube[clipped == 0] = 0
    return emiscube




#hdulist = fits.open('C:/Users/trist/Desktop/KinMSpy-master/test_suite/NGC5806.CO.cube.fits')
hdulist = fits.open('C:/Users/trist/Desktop/KinMSpy-master/test_suite/NGC3607_cube_nopbcorr.fits')
cube = hdulist[0].data.T
#cube = cube[:,:,:,0]
#BeamData = hdulist[1].data
cutoff = 2
header = hdulist[0]
start = 20
stop = 80
emiscube = cube.copy()[:,:,start:stop]


Mask = smoothclip(cube, emiscube, cutoff, header, start, stop)




#TotMask = np.sum(Mask, axis = 2)
#X = np.arange(0, 100, 1)
#Y = np.arange(0, 100, 1)
#X, Y = np.meshgrid(X, Y)
#plt.figure('Mask')
#plt.contourf(X, Y, TotMask.T, linewidth=0, antialiased=False)












