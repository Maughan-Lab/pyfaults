import numpy as np

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' creates a single stacking vectory (sx, sy) '''

class S_vector(list):
    def __init__(self, s_x, s_y):
        '''
        s_x (float)
        s_y (float)
        '''
        
        # assign s_x and s_y stacking vectors
        self.s_x = s_x
        self.s_y = s_y
        
        self.append(s_x, s_y)

        return
#------------------------------------------------------------------------------  
    ''' S_vector methods '''   
    
    # randomly generate a new stacking vector
    def rng(self, s_max):
        '''
        s_max (float)
        '''
        
        new_sx = np.random.Generator.uniform(self, 0, s_max, 1)
        new_sy = np.random.Generator.uniform(self, 0, s_max, 1)
        
        self.s_x = new_sx
        self.s_y = new_sy
        
    def show(self):
        print(self)
        
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' creates a grid of stacking vectors based on a given step size '''

class v_grid(list):
    def __init__(self, s_max, step):
        '''
        s_max (float)
        step (float)
        '''
        
        v_pairs = []
        
        num = (s_max / step)
        
        for n in num:
            s = []
            sx = np.linspace(0, s_max, retstep=step, dtype=float)
            sy = np.linspace(0, s_max, retstep=step, dtype=float)
            s.append(sx, sy)
            v_pairs.append(s)
            
        return

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
''' creates a RNG grid of stacking vectors '''
    
class rng_grid(list):
    def __init__(self, s_max, num_v):
        '''
        s_max (float)
        num_v (float)
        '''
        
        v_pairs = []
        
        for n in num_v:
            s = []
            sx = np.random.Generator.uniform(self, 0, s_max, 1)
            sy = np.random.Generator.uniform(self, 0, s_max, 1)
            s.append(sx, sy)
            v_pairs.append(s)
            
        return