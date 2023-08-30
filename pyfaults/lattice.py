#########################################################################################
# pyfaults.lattice
# Author: Sinclair R. Combs
#########################################################################################

''' Lattice object class -- Stores unit cell lattice parameters '''
#----------------------------------------------------------------------------------------
class Lattice(object):
    '''
    Parameters
    ----------
    a : float
        Unit cell lattice vector a (A)
    b : float
        Unit cell lattice vector b (A)
    c : float
        Unit cell lattice vector c (A)
    alpha : float
        Unit cell lattice angle alpha (deg)
    beta : float
        Unit cell lattice angle beta (deg)
    gamma : float
        Unit cell lattice angle gamma (deg)
    '''
    
    # properties ------------------------------------------------------------------------
    a = property(lambda self: self._a,
                 lambda self, val: self.setParam(a=val),
                 doc="float : unit cell vector a")
    
    b = property(lambda self: self._b,
                 lambda self, val: self.setParam(b=val),
                 doc="float : unit cell vector b")
    
    c = property(lambda self: self._c,
                 lambda self, val: self.setParam(c=val),
                 doc="float : unit cell vector c")
    
    alpha = property(lambda self: self._alpha,
                 lambda self, val: self.setParam(alpha=val),
                 doc="float : unit cell angle alpha")
    
    beta = property(lambda self: self._beta,
                 lambda self, val: self.setParam(beta=val),
                 doc="float : unit cell angle beta")
    
    gamma = property(lambda self: self._gamma,
                 lambda self, val: self.setParam(gamma=val),
                 doc="float : unit cell angle gamma")
    
    # creates instance of Lattice object ------------------------------------------------
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
    
    # sets lattice parameters -----------------------------------------------------------
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
    
    # prints lattice parameters ---------------------------------------------------------
    def display(self):
        latt = "\n".join(("a : " + self.a, 
                       "b : " + self.b,
                       "c : " + self.c,
                       "alpha : " + self.alpha,
                       "beta : " + self.beta,
                       "gamma : " + self.gamma))
        print(latt)