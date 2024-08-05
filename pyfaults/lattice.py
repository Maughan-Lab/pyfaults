##################################################################################
# Author: Sinclair R. Combs
##################################################################################

# Lattice object class ----------
class Lattice(object):
    '''
    Parameters
    ----------
    a (float) : unit cell vector a
    b (float) : unit cell vector b
    c (float) : unit cell vector c
    alpha (float) : unit cell angle alpha
    beta (float) : unit cell angle beta
    gamma (float) : unit cell angle gamma
    '''
    # properties ----------
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
    
    # initialization, defines Lattice defaults ----------
    def __init__(self, a, b, c, alpha, beta, gamma):
        self._a = None
        self._b = None
        self._c = None
        self._alpha = None
        self._beta = None
        self._gamma = None
        self.setParam(a, b, c, alpha, beta, gamma)
        return
    
    # sets parameters of Lattice ----------
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
    
    # prints lattice parameters ----------
    def display(self):
        displayStr = '\n'.join(('a : ' + str(self.a), 
                         'b : ' + str(self.b), 
                         'c : ' + str(self.c),
                         'alpha : ' + str(self.alpha), 
                         'beta : ' + str(self.beta),
                         'gamma : ' + str(self.gamma)))
        print(displayStr)
        return
