from scipy import ndimage, interpolate
import numpy as np
from makebeam import makebeam
import matplotlib.pyplot as plt
from sauron_colormap import sauron
from astropy.convolution import convolve_fft



#Masking Code
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
        emiscube[emiscube<cutoff*sigma] = 0
        return emiscube
def smoothclip(cube, cutoff, beamsize, cellsize, start, stop, psf):
        '''
        Apply a Gaussian smoothing to the data cube.
        cube = the raw data cube
        emiscube = the part of the cube that contains the emission line
        cutoff = the level at which you want to clip (in multiples of the rms noise of the cube)
        beamsize = Array of [Bmin, Bmaj, Bpa]
        start = first channel that contains emission
        stop = last channel that contains emission    
        '''
        #copy the datacube
        cube_copy = cube.copy() #a copy of the data cube will be used to do the smoothing. This smoothed cube will then be used as a mask to clip the orignial cube
        emiscube = cube.copy()[:,:,start:stop]
        # extract relevant information from header (pixel size, beam size)
        beam = beamsize[1] / cellsize  #beam size in pixels, use the major axis
        # convert the FWHM of the beam to the std dev in the Gaussian distribution and use 1.5 times the beam size to convolve with
        sigma = 1.5 * beam / np.sqrt(8.*np.log(2.))
        #apply a Gaussian blur, using sigma = 4 in the velocity direction (seems to work best). The mode 'nearest' seems to give the best results.
#        cube_smoothed = ndimage.filters.gaussian_filter(cube_copy, (4.,sigma,sigma), order=0, mode='nearest')
        cube_smoothed = cube_copy.copy() * 0
        w2do = np.where(cube.sum(axis=0).sum(axis=0) >0)[0]
        for i in range(0,w2do.size): cube_smoothed[:,:,w2do[i]] = convolve_fft(cube[:,:,w2do[i]], psf)
        
        #clip the copied cube, using the smaller cube where the emission is
        smoothemis, smoothnoise = emisChannels(cube_smoothed, start, stop)  #this separates the channels with emission from the ones without emission and glues the latter bits together in one "noise cube"
        clipped = clip(smoothemis, smoothnoise, cutoff)
        #mask the original cube with the smoothed cube
        emiscube[clipped == 0] = 0
        return emiscube
  
#Plotting 
def makeplots(f,xsize,ysize,vsize,cellsize,dv,beamsize,posang=0,overcube=False,pvdthick=2,nconts=11.,title=False, **kwargs):
    
# ;;;; Create plot data from cube ;;;;
    mom0rot=f.sum(axis=2)
    if np.any(overcube): mom0over=overcube.sum(axis=2)
    x1=np.arange(-xsize/2.,xsize/2.,cellsize)
    y1=np.arange(-ysize/2.,ysize/2.,cellsize)
    v1=np.arange(-vsize/2.,vsize/2.,dv)

    mom1=(mom0rot*0.0)-10000.0
    for i in range(0,int(xsize/cellsize)):
         for j in range(0,int(ysize/cellsize)):
             if mom0rot[i,j] > 0.1*np.max(mom0rot):
                 mom1[i,j]=(v1*f[i,j,:]).sum()/f[i,j,:].sum()
    
    pvdcube=f 
    
    pvdcube=ndimage.interpolation.rotate(f, 90-posang, axes=(1, 0), reshape=False)
    if np.any(overcube): pvdcubeover=ndimage.interpolation.rotate(overcube, 90-posang, axes=(1, 0), reshape=False)
        
    pvd=pvdcube[:,np.int((ysize/2.)-pvdthick):np.int((ysize/2.)+pvdthick),:].sum(axis=1)
    if np.any(overcube): pvdover=pvdcubeover[:,np.int((ysize/2.)-pvdthick):np.int((ysize/2.)+pvdthick),:].sum(axis=1)
    
    if not isinstance(beamsize, (list, tuple, np.ndarray)):
        beamsize=np.array([beamsize,beamsize,0])
    beamtot=(makebeam(xsize,ysize,[beamsize[0]/cellsize,beamsize[1]/cellsize],rot=beamsize[2])).sum()
    #       (makebeam(xs   ,ys   ,[beamsize[0]/cellsize,beamsize[1]/cellsize],rot=beamsize[2])).sum()
    spec=f.sum(axis=0).sum(axis=0)/beamtot
    if np.any(overcube): specover=overcube.sum(axis=0).sum(axis=0)/beamtot    
# ;;;;

# ;;;; Plot the results ;;;;
    levs=v1[np.min(np.where(spec != 0)):np.max(np.where(spec != 0))]
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(221, aspect='equal')
    plt.xlabel('Offset (")')
    plt.ylabel('Offset (")')
    ax1.contourf(x1,y1,mom0rot.T,levels=np.linspace(1,0,num=10,endpoint=False)[::-1]*np.max(mom0rot), cmap="YlOrBr")
    if np.any(overcube): ax1.contour(x1,y1,mom0over.T,colors=('black'),levels=np.arange(0.1, 1.1, 0.1)*np.max(mom0over))
    if 'yrange' in kwargs: ax1.set_ylim(kwargs['yrange'])
    if 'xrange' in kwargs: ax1.set_xlim(kwargs['xrange'])
    ax2 = fig.add_subplot(222, aspect='equal')
    plt.xlabel('Offset (")')
    plt.ylabel('Offset (")')
    ax2.contourf(x1,y1,mom1.T,levels=levs, cmap=sauron)
    if 'yrange' in kwargs: ax2.set_ylim(kwargs['yrange'])
    if 'xrange' in kwargs: ax2.set_xlim(kwargs['xrange'])
    ax3 = fig.add_subplot(223)
    plt.xlabel('Offset (")')
    plt.ylabel(r'Velocity (km s$^{-1}$)')
    ax3.contourf(x1,v1,pvd.T,levels=np.linspace(1,0,num=10,endpoint=False)[::-1]*np.max(pvd), cmap="YlOrBr" ,aspect='auto')
    if np.any(overcube): ax3.contour(x1,v1,pvdover.T,colors=('black'),levels=np.arange(0.1, 1.1, 0.1)*np.max(pvdover))
    if 'vrange' in kwargs: ax3.set_ylim(kwargs['vrange'])
    if 'xrange' in kwargs: ax3.set_xlim(kwargs['xrange'])
    ax4 = fig.add_subplot(224)
    plt.ylabel('Flux')
    plt.xlabel(r'Velocity (km s$^{-1}$)')
    ax4.plot(v1,spec, drawstyle='steps')
    if np.any(overcube): ax4.plot(v1,specover,'r', drawstyle='steps')
    if 'vrange' in kwargs: ax4.set_xlim(kwargs['vrange'])
    if title: plt.suptitle(title)
    plt.show()