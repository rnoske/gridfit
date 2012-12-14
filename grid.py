# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import scipy.ndimage as spimg
import scipy.signal as spsig
import matplotlib.pyplot as plt

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
        
def grid(ov, oh, s, h, w):
    """ Return a 2d grid
    
    h = height of image
    w = width of image
    ov = origin vertical
    origin horizontal
    s = spacing
    
    """
    g = np.zeros((h,w)) #nparray full of zeros
    p = np.zeros((h,w)) #nparray full of zeros
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
                g[_h, int(_pos)] = (1-(_pos - int(_pos))) *255
                g[_h, int(_pos+1)] = (1-(int(_pos+1)-_pos)) *255
                p[_h, int(_pos)] += 1
                #p[_h, int(_pos)] += (1-(_pos - int(_pos))) *1
                #p[_h, int(_pos+1)] += (1-(int(_pos+1)-_pos)) *1
            except IndexError:
                pass
    
    _vp = oh
    while _vp >= s:
        _vp = _vp -s
        
    for _nh in xrange(nh):
        _pos = _nh * s + _vp
        for _v in xrange(w):
            try:
                g[int(_pos), _v] = (1-(_pos - int(_pos))) *255
                g[int(_pos+1), _v] = (1-(int(_pos+1)-_pos)) *255
                p[int(_pos), _v] += 1
                #p[int(_pos), _v] += (1-(_pos - int(_pos))) *1
                #p[int(_pos+1), _v] += (1-(int(_pos+1)-_pos)) *1
            except IndexError:
                pass
    """        
    points = []
    for _h in xrange(h):
        for _v in xrange(w):
            if p[_h, _v] <= 1:
                p[_h, _v] = 0
            else:
                p[_h, _v] = 255
                points.append((_h, _v))
    """ 
    
    return g#,p, points #g=grid p=crossingpoints, points=list of crossingpoints
        
def calibrate_image():
    """ Calibrate an image with grid function
    
    """
    fp = 'D:/Raimund Buero/Python/SpyDev/gridfit/3cm.fits'
    name = '3cm'
    img, header = read_fits_nparray(name=fp)
    #img = spimg.imread(fp, flatten=True)
    img = spimg.interpolation.rotate(img, -0.9, order = 5, reshape=False)
    img = img[380:985, :] #cut interesting image part
    h, w = img.shape
    img = spimg.filters.median_filter(img, size=(3,3))
    7
    ov, oh = (520.3, 13.5) #origin of grid (middle bottom of target)13,5
    s = 15.65 #spacing in px for grid
    
    """
    img geschnitten: img[380:985, :]
    3cm:
        ov, oh = (520.3, 13.5)
        s = 15.65
    0cm:
        ov, oh = (518.5, 6.5)
        s = 16.4
    6cm:
        ov, oh = (519.5, 19)
        s = 14.9
    """
    
    g = grid(ov,oh,s,h,w)
    b = np.zeros(img.shape)

    img *= (255/img.max())
    g *= (255/g.max())
    
    rgb = np.dstack((img,g,b))
    #rbb = np.dstack((img,b,b))
    #showing and saving
    plt.cla()
    plt.clf()
    plt.imshow(img)
    plt.colorbar()
    #plt.savefig('fit.jpg')
    #sp.misc.imsave('name-original.jpg', img)
    #plt.show()
    
    plt.cla()
    plt.clf()
    plt.imshow(img)
    #plt.hold(True)
    plt.imshow(g, alpha=0.5)
    #plt.colorbar()
    plt.savefig('fit.jpg')
    #sp.misc.imsave('grid.jpg', g)
    plt.show()
    
    from PIL import Image
    _img = Image.fromarray(np.uint8(img))
    _img.save(name + '-ori.bmp')
    _img = Image.fromarray(np.uint8(rgb))
    _img.save(name + '-ori+grid.bmp')
    
def calc_distance(points1, points2):
    p1 = points1
    p2 = points2
    dist = []
    for _p1 in p1:
        dist_old = 99999999
        h_p1, v_p1 = _p1
        for _p2 in p2:
            h_p2, v_p2 = _p2
            dist_new = np.sqrt((h_p2-h_p1)**2 + (v_p2-v_p1)**2)
            if dist_new < dist_old:
                dist_old = dist_new
                dist_p = h_p2, v_p2
        h_p2, v_p2 = dist_p
        dist.append([h_p1, v_p1, h_p2, v_p2])
    
    return dist
    
def calc_values(s, new_origin, h, w):
    v = np.zeros((h,w))
    new_h, new_w = new_origin
    for _h in xrange(h):
        for _w in xrange(w):
            v[_h, _w] = np.sqrt(((_h-new_h)**2)+((_w-new_w)**2))/s
    
    return v
    
def calc_steigung(dstackarray):
    stack = dstackarray
    h, w, z = stack.shape
    #print stack.shape
    m = np.zeros((h,w))
    x = np.arange(0,z,1)
    from scipy import polyfit
    for _h in xrange(h):
        for _w in xrange(w):
            a, b = polyfit(x, stack[_h,_w,:],1)
            m[_h, _w] = a
    return m
        

if __name__ == "__main__":
    #calibrate_image()
    h, w = (605, 1024) #size of cutted image
    
    new_origin = (512-380, 512) #old origin of image ist 512,512
    
    p_0 = calc_values(16.4, new_origin, h, w)
    p_3 = calc_values(15.65, new_origin, h, w)
    p_6 = calc_values(14.9, new_origin, h, w)
    
    stack = np.dstack((p_0, p_3, p_6))
    m = calc_steigung(stack)
    
    plt.imshow(m)
    plt.colorbar()
    plt.savefig('fsteigung.jpg')
    plt.show()
    
    from PIL import Image
    m /= np.amax(m)
    m *= 255
    _img = Image.fromarray(np.uint8(m))
    _img.save('v5.bmp')
    
    
    """
    g_0, p_0, po_0= grid(518.5, 6.5, 16.4, h, w)
    g_3, p_3, po_3 = grid(520.3, 13.5, 15.65, h, w)
    g_6, p_6, po_6 = grid(519.5, 19, 14.9, h, w)

    from PIL import Image, ImageDraw
    rgb = np.dstack((g_0, g_3, g_6))
    _img = Image.fromarray(np.uint8(rgb))
    _img.save('vergleich1.bmp')

    rgb = np.dstack((p_0, p_3, p_6))
    _img = Image.fromarray(np.uint8(rgb))
    _img.save('vergleich2.bmp')
    
    dist = calc_distance(po_0, po_6)
    draw = ImageDraw.Draw(_img)
    for entry in dist:
        draw.line(entry, fill = 128)
    _img.save('vergleich4.bmp')
    """
        
    
    