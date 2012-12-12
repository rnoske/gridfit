# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.ndimage as spimg
import scipy.signal as spsig
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_errors

def read_fits_nparray(name = 'test.fit', number = 0):
        """ Read .fits file from iStar camera
        
        name (str): file name
        number (int): number of hdulist (usually 0)
        
        Returns:
            _header (pyfits.header.Header): dictionary type something
            _arr (numpy.ndarray): numpy array
        
        """
        import pyfits
        _file = name #self. workspace + name
        _fits = pyfits.open(_file)
        _header = _fits[number].header
        _arr = _fits[number].data
        _arr = _arr[0,:,:] #da nur 2D Array
        return _arr, _header
        
def grid_fit(nparray, params, plotflag = True):
    """ Fit a grid to a 2D data set
    
    """
    def grid(h, w, nv, nh, dv, dh, ov, oh, a):
        """ Return a 2d grid
        
        h = height of image
        w = width of image
        nv = number of vertical lines
        nh = number of horizontal lines
        dv = spacing vertical
        dh = spacing horizontal
        ov = offset vertical
        oh = offset horizontal
        a = amplitude of grid
        
        """
        g = np.zeros((h,w)) #nparray full of zeros
        
        for _h in xrange(h): # for every height
            for _nv in xrange(nv): #for every vertical line
                pos = _nv*dv+ov #position jeder verticalen linie
                try:
                    g[_h, int(pos)] = (pos - int(pos)) *a
                    g[_h, int(pos+1)] = (int(pos)+1 - pos) *a
                except IndexError:
                    pass
                
        for _w in xrange(w): #for every width
            for _nh in xrange(nh): #for every horizontal line
                pos = _nh*dh+oh
                try:
                    g[int(pos), _w] = (pos - int(pos)) *a
                    g[int(pos)+1, _w] = (int(pos)+1 - pos)*a
                except IndexError:
                    pass
                
        return g
    
    def res(params, nparray):
        nv = int(params['nv'].value)
        nh = int(params['nh'].value)
        dv = params['dv'].value
        dh = params['dh'].value
        ov = params['ov'].value
        oh = params['oh'].value
        a = params['a'].value
        h = nparray.shape[0]
        w = nparray.shape[1]
        b = params['b'].value
        barr = np.ones((h,w))*b
        
        model = grid(h, w, nv, nh, dv, dh, ov, oh, a)
        
        err = nparray + model - barr
        err = err.flatten()
        
        return err
    
    #do fit
    minimize(res, params, args=(nparray,)) 
    #komma ist wichtig, da tupel erwartet!
    
    if plotflag == True:
        nv = int(params['nv'].value)
        nh = int(params['nh'].value)
        dv = params['dv'].value
        dh = params['dh'].value
        ov = params['ov'].value
        oh = params['oh'].value
        a = params['a'].value
        h = nparray.shape[0]
        w = nparray.shape[1]
        b = params['b'].value
        barr = np.ones((h,w))*b
        
        g = grid(h, w, nv, nh, dv, dh, ov, oh, a)
        sp.misc.imsave('grid.jpg', g)
        fit = nparray+g-barr
        plt.cla()
        plt.clf()
        plt.imshow(fit)
        plt.colorbar()
        #plt.savefig('fit.jpg')
        sp.misc.imsave('fit.jpg', fit)
        plt.show()
        
    return params
    

if __name__ == "__main__":
    fp = 'D:/Raimund Buero/Python/SpyDev/gridfit/3cm.fits'
    img, header = read_fits_nparray(name=fp)
    #img = spimg.imread(fp, flatten=True)
    img = spimg.interpolation.rotate(img, -0.9, order = 5, reshape=False)
    img = img[401:714, 364:676]
    img = spimg.filters.median_filter(img, size=(3,3))
    sobel_x = [[-1, 0, 1],[-2,0,2],[-1,0,1]]
    sobel_y = [[-1,-2,-1],[0,0,0],[1,2,1]]
    #img = spimg.convolve(img, sobel_x)
    #img = spimg.convolve(img, sobel_y)
    
    #img = spsig.medfilt2d(img, kernel_size=3)
    #img = img * (-1)
    
    #show image
    plt.cla()
    plt.clf()
    plt.imshow(img, origin='lower')
    plt.colorbar()
    #plt.savefig('original.tif')
    sp.misc.imsave('original.jpg', img)
    plt.show()
    
    # create a set of Parameters
    params = Parameters()
    params.add('nv', value=20, vary=False)
    params.add('nh', value=19, vary=False)
    params.add('dv', value=16, vary=True, min=5.0, max=25.0)
    #params.add('dh', value=15, vary=True, min=5.0, max=25.0)
    params.add('dh', expr='dv')
    params.add('ov', value=14, vary=True, min=0., max=50)
    params.add('oh', value=8, vary=True, min=1., max=50)
    params.add('a', value=0, vary=True)
    params.add('b', value=30000, vary=True)
    
    #do fit
    #grid_fit(img, params)
    
    print params['dv'].value, params['a'].value, params['b'].value, params['ov'].value, params['oh'].value