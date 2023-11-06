##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np

#---------------------------------------------------------------------------------
# calculates difference curve ----------------------------------------------------
#---------------------------------------------------------------------------------
def diffCurve(q1, q2, ints1, ints2): 
    '''
    Parameters
    ----------
    q1
        nparray : dataset 1 Q values (A^-1)
    q2
        nparray : dataset 2 Q values (A^-1)
    ints1
        nparray : dataset 1 intensity values
    ints2
        nparray : dataset 2 intensity values

    Returns
    -------
    diff_q : nparray
    diff_ints : nparray
    '''
    diff_q_list = []
    diff_ints_list = []
    for i in range(0, len(q1)):
        q1_val = float('%.3f'%(q1[i]))
        for j in range(0, len(q2)):
            q2_val = float('%.3f'%(q2[j]))
            if q1_val == q2_val:
                diff_q_list.append(q1_val)
                diff_ints_list.append(ints1[i]-ints2[j])
    diff_q = np.array(diff_q_list)
    diff_ints = np.array(diff_ints_list)
    return diff_q, diff_ints

#---------------------------------------------------------------------------------
# saves difference curve as text file --------------------------------------------
#---------------------------------------------------------------------------------
def saveDiffCurve(q, ints, path, fn):
    '''
    Parameters
    ----------
    q
        nparray : Q values (A^-1)
    ints
        nparray : intensity difference values
    path
        str : file save directory
    fn
        str : File name
    '''
    with open(path + fn + '.txt', 'w') as f:
        for (q, ints) in zip(q, ints):
            f.write('{0} {1}\n'.format(q, ints))
    f.close()