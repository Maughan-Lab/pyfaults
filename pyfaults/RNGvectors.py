#########################################################################################
# pyfaults.RNGvectors
# Author: Sinclair R. Combs
#########################################################################################

import random as r
import numpy as np

''' RNGvectors class -- Lists randomly generated stacking vectors '''
#----------------------------------------------------------------------------------------

class RNGvectors(list):
    '''
    Parameters
    ----------
    num : int
        Number of stacking vectors to generate
    xRange : list (float), optional
        Tuple of minimum and maximum x-values. The default is (0,0).
    yRange : list (float), optional
        Tuple of minimum and maximum y-values. The default is (0,0).
    zRange : list (float), optional
        Tuple of minimum and maximum z-values. The default is (0,0).
    numDecimals : int, optional
        Specify number of decimal places in vector components. The default is 2.
    '''
        
    # properties ------------------------------------------------------------------------
    num = property(lambda self: self._num,
                   lambda self, val: self.setParam(num=val))
        
    xRange = property(lambda self: self._xRange,
                      lambda self, val: self.setParam(xRange=val))
        
    yRange = property(lambda self: self._yRange,
                      lambda self, val: self.setParam(yRange=val))
        
    zRange = property(lambda self: self._zRange,
                      lambda self, val: self.setParam(zRange=val)) 
    
    numDecimals = property(lambda self: self._numDecimals,
                           lambda self, val: self.setParam(numDecimals=val))
    
    vecList = property(lambda self: self._vecList)
        
    # creates instance of RNGvectors object ---------------------------------------------
    def __init__(self, num, xRange=None, yRange=None, zRange=None, numDecimals=None):
        
        # initialize parameters
        self._num = None
        self._xRange = None
        self._yRange = None
        self._zRange = None
        self._numDecimals = None
        
        self.setParam(num, xRange, yRange, zRange, numDecimals)
        
        self._vecList = self.setVec()
        
        return
    
    # set parameters --------------------------------------------------------------------
    def setParam(self, num, xRange=None, yRange=None, zRange=None, numDecimals=None):
        
        self._num = num
        
        if xRange is not None:
            self._xRange = xRange
        elif xRange is None:
            self._xRange = (0,0)
        
        if yRange is not None:
            self._yRange = yRange
        elif yRange is None:
            self._yRange = (0,0)
        
        if zRange is not None:
            self._zRange = zRange
        elif zRange is None:
            self._zRange = (0,0)
            
        if numDecimals is not None:
            self._numDecimals = numDecimals
        elif numDecimals is None:
            self._numDecimals = 2
        
        return
          
    def setVec(self):
        vecList = []
        for i in range(self.num):
            if self.xRange != (0,0):
                xval = round(r.uniform(self.xRange[0], self.xRange[1]), self.numDecimals)
            elif self.xRange == (0,0):
                xval = 0
            if self.yRange != (0,0):
                yval = round(r.uniform(self.yRange[0], self.yRange[1]), self.numDecimals)
            elif self.yRange == (0,0):
                yval = 0
            if self.zRange != (0,0):
                zval = round(r.uniform(self.zRange[0], self.zRange[1]), self.numDecimals)
            elif self.zRange == (0,0):
                zval = 0
                
            s = np.array([xval, yval, zval])
            vecList.append(s)
            
        return vecList
    
    
    
    
    
    
    
    
    
    
    
    