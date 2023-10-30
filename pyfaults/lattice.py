##################################################################################
# Author: Sinclair R. Combs
##################################################################################

''' LATTICE CLASS '''
#---------------------------------------------------------------------------------
class Lattice(object):
    # properties -----------------------------------------------------------------
    a =\
        property(lambda self: self._a,
                 lambda self, val: self.setParam(a=val),
                 doc='float : unit cell vector a in Angstroms')
    
    b =\
        property(lambda self: self._b,
                 lambda self, val: self.setParam(b=val),
                 doc='float : unit cell vector b in Angstroms')
    
    c =\
        property(lambda self: self._c,
                 lambda self, val: self.setParam(c=val),
                 doc='float : unit cell vector c in Angstroms')
    
    alpha =\
        property(lambda self: self._alpha,
                 lambda self, val: self.setParam(alpha=val),
                 doc='float : unit cell angle alpha in degrees')
    
    beta =\
        property(lambda self: self._beta,
                 lambda self, val: self.setParam(beta=val),
                 doc='float : unit cell angle beta in degrees')
    
    gamma =\
        property(lambda self: self._gamma,
                 lambda self, val: self.setParam(gamma=val),
                 doc='float : unit cell angle gamma in degrees')
    
    # creates instance of Lattice object -----------------------------------------
    def __init__(self, a, b, c, alpha, beta, gamma):
        # initialize parameters
        self._a = None
        self._b = None
        self._c = None
        self._alpha = None
        self._beta = None
        self._gamma = None
        self.setParam(a, b, c, alpha, beta, gamma)
        return
    
    # sets parameters ------------------------------------------------------------
    def setParam(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None):
        if a is not None:
            self._a = a
        if b is not None:
            self._b = b
        if c is not None:
            self._c = c
        if alpha is not None:
            self._alpha = alpha
        if beta is not None:
            self._beta = beta
        if gamma is not None:
            self._gamma = gamma
        return
    
    # prints lattice parameters --------------------------------------------------
    def display(self):
        displayStr = '\n'.join(('a : ' + self.a, 
                         'b : ' + self.b, 
                         'c : ' + self.c,
                         'alpha : ' + self.alpha, 
                         'beta : ' + self.beta,
                         'gamma : ' + self.gamma))
        print(displayStr)
        return displayStr