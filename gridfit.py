# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.ndimage as spimg
import scipy.signal as spsig
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_errors

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
        
        model = grid(h, w, nv, nh, dv, dh, ov, oh, a)
        
        err = nparray - model
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
        
        g = grid(h, w, nv, nh, dv, dh, ov, oh, a)
        fit = nparray-g
        plt.cla()
        plt.clf()
        plt.imshow(fit)
        plt.colorbar()
        #plt.savefig('fit.jpg')
        sp.misc.imsave('fit.jpg', fit)
        plt.show()
        
    return params
    
if __name__ == "__main__":
    fp = 'C:/Python/SpyDev/gridfit/tar_mitte.jpg'
    img = spimg.imread(fp, flatten=True)
    img = spimg.interpolation.rotate(img, 0.9, order = 5, reshape=False)
    img = img[5:390, 5:320]
    #img = spimg.filters.median_filter(img, size=(3,3))
    img = spsig.medfilt2d(img, kernel_size=3)
    #img = img * (-1)
    
    #show image
    plt.cla()
    plt.clf()
    plt.imshow(img)
    plt.colorbar()
    sp.misc.imsave('original.jpg', img)
    #plt.show()
    
    # create a set of Parameters
    params = Parameters()
    params.add('nv', value=19, vary=True)
    params.add('nh', value=24, vary=True)
    params.add('dv', value=15, vary=True, min=5.0, max=25.0)
    #params.add('dh', value=15, vary=True, min=5.0, max=25.0)
    params.add('dh', expr='dv')
    params.add('ov', value=12, vary=True, min=0., max=50)
    params.add('oh', value=16, vary=True, min=1., max=50)
    params.add('a', value=10, vary=True)
    
    #do fit
    grid_fit(img, params)
    
    print params['dv'].value, params['a'].value, params['ov'].value, params['oh'].value