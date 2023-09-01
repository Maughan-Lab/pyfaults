#########################################################################################
# pyfaults.vectorGen
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
          
    # generate vectors ------------------------------------------------------------------
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
  

''' vecGridGen class -- Lists combinations of stacking vectors generated stepwise '''
#----------------------------------------------------------------------------------------    

class vecGridGen(list):
    
    # properties ------------------------------------------------------------------------ 
    xRange = property(lambda self: self._xRange,
                      lambda self, val: self.setParam(xRange=val))
        
    yRange = property(lambda self: self._yRange,
                      lambda self, val: self.setParam(yRange=val))
        
    zRange = property(lambda self: self._zRange,
                      lambda self, val: self.setParam(zRange=val)) 
    
    numDecimals = property(lambda self: self._numDecimals,
                           lambda self, val: self.setParam(numDecimals=val))
    
    vecList = property(lambda self: self._vecList)
    
    def __init__(self, xRange=None, yRange=None, zRange=None, numDecimals=None):
        
        # initialize parameters
        self._xRange = None
        self._yRange = None
        self._zRange = None
        self._numDecimals = None
        
        self._xStep = None
        self._yStep = None
        self._zStep = None
        
        self._xNum = None
        self._yNum = None
        self._zNum = None
        
        self.setParam(xRange, yRange, zRange, numDecimals)
        
        self._vecList = self.setVec()
        
        return
    
    # set parameters --------------------------------------------------------------------
    def setParam(self, xRange=None, yRange=None, zRange=None, numDecimals=None):
        
        if xRange is not None:
            self._xRange = [xRange[0], xRange[1]]
            self._xStep = xRange[2]
            self._xNum = int((xRange[1] - xRange[0]) / xRange[2])
        elif xRange is None:
            self._xRange = (0,0)
            self._xStep = 0
            self._xNum = 1
        
        if yRange is not None:
            self._yRange = [yRange[0], yRange[1]]
            self._yStep = yRange[2]
            self._yNum = int((yRange[1] - yRange[0]) / yRange[2])
        elif yRange is None:
            self._yRange = (0,0)
            self._yStep = 0
            self._yNum = 1
        
        if zRange is not None:
            self._zRange = [zRange[0], zRange[1]]
            self._zStep = zRange[2]
            self._zNum = int((zRange[1] - zRange[0]) / zRange[2])
        elif zRange is None:
            self._zRange = (0,0)
            self._zStep = 0
            self._zNum = 1
            
        if numDecimals is not None:
            self._numDecimals = numDecimals
        elif numDecimals is None:
            self._numDecimals = 2
        
        return
    
    # generate vectors ------------------------------------------------------------------
    def setVec(self):
        vecList = []
        
        xList = []
        if self.xNum > 1:
            for i in range(self.xNum):
                addX = self.xRange[0] + (self.xStep * i)
                xList.append(addX)
        elif self.xNum == 1:
            xList.append(self.xRange[0])
            
        yList = []
        if self.yNum > 1:
            for i in range(self.yNum):
                addY = self.yRange[0] + (self.yStep * i)
                yList.append(addY)
        elif self.yNum == 1:
            yList.append(self.yRange[0])
            
        zList = []
        if self.zNum > 1:
            for i in range(self.zNum):
                addZ = self.zRange[0] + (self.zStep * i)
                zList.append(addZ)
        elif self.zNum == 1:
            zList.append(self.zRange[0])
            
        for x in xList:
            for y in yList:
                for z in zList:
                    newVec = np.array([x, y, z])
        vecList.append(newVec)
            
        return vecList
    
    
    
    
    
    
    
    
    
    
    