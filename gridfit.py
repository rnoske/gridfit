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
    def grid(ov, oh, s, h, w, a):
        """ Return a 2d grid
        
        h = height of image
        w = width of image
        ov = origin vertical
        oh = origin horizontal
        s = spacing
        a = amplitude
        
        """
        g = np.zeros((h,w)) #nparray full of zeros
        nv = int(w/s +2) #+2 just to make sure
        nh = int(h/s+2)
        
        #calculate grid origin offset
        _hp = ov
        while _hp >= s:
            _hp = _hp - s

        #print _hp
        for _nv in xrange(nv): # for every vertical line
            #calculate horizontal position   
            _pos = _nv * s + _hp
            #print _pos
            for _h in xrange(h): #for every height
                try:
                    g[_h, int(_pos)] = (1-(_pos - int(_pos))) *a
                    g[_h, int(_pos+1)] = (1-(int(_pos+1)-_pos)) *a
                except IndexError:
                    pass
        
        _vp = oh
        while _vp >= s:
            _vp = _vp -s
            
        for _nh in xrange(nh):
            _pos = _nh * s + _vp
            for _v in xrange(w):
                try:
                    g[int(_pos), _v] = (1-(_pos - int(_pos))) *a
                    g[int(_pos+1), _v] = (1-(int(_pos+1)-_pos)) *a
                except IndexError:
                    pass
           
        return g
    
    def res(params, nparray):
        s = params['s'].value
        ov = params['ov'].value
        oh = params['oh'].value
        a = params['a'].value
        h = nparray.shape[0]
        w = nparray.shape[1]
        b = params['b'].value
        barr = np.ones((h,w))*b
        
        model = grid(ov, oh, s, h, w, a)
        
        err = nparray - (barr-model)
        err = err.flatten()
        
        return err
    
    #do fit
    minimize(res, params, args=(nparray,))
    #komma ist wichtig, da tupel erwartet!
    report_errors(params)
    
    if plotflag == True:
        s = params['s'].value
        ov = params['ov'].value
        oh = params['oh'].value
        a = params['a'].value
        h = nparray.shape[0]
        w = nparray.shape[1]
        b = params['b'].value
        barr = np.ones((h,w))*b
        
        g = grid(ov, oh, s, h, w, a)
        #sp.misc.imsave('grid.jpg', g)
        #fit = nparray+g-barr
        plt.cla()
        plt.clf()
        plt.imshow(img)
        #plt.hold(True)
        plt.imshow(g, alpha=0.5)
        #plt.colorbar()
        #plt.savefig('fit.jpg')
        sp.misc.imsave('grid.jpg', g)
        plt.show()
        
        
        
    return params
    

if __name__ == "__main__":
    fp = 'D:/Raimund Buero/Python/SpyDev/gridfit/3cm.fits'
    img, header = read_fits_nparray(name=fp)
    #img = spimg.imread(fp, flatten=True)
    img = spimg.interpolation.rotate(img, -0.9, order = 5, reshape=False)
    #img = img[401:714, 364:676]
    img = img[380:985, :] #cut interesting image part
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
    #plt.show()
    
    # create a set of Parameters
    params = Parameters()    
    params.add('ov', value=520, vary=True)
    params.add('oh', value=28, vary=True)
    params.add('s', value=15.7, vary=True, min=10.0, max=20.0)#15.7
    params.add('a', value=5000, vary=True, min = 0, max=20000)
    params.add('b', value=30000, vary=True, min=10000)
    
    #do fit
    grid_fit(img, params)
    
    #print params['s'].value, params['a'].value, params['b'].value, params['ov'].value, params['oh'].value