##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np

# get parameters for normalization ----------
def getNormVals(q, ints):
    
    intsMax = 0
    maxIndex = 0
    
    for i in range(len(ints)):
        if ints[i] > intsMax:
            intsMax = ints[i]
            maxIndex = i
            
    qAtIntsMax = q[maxIndex]
    
    intsMin = np.min(ints)
    
    return intsMax, qAtIntsMax, maxIndex, intsMin


# normalize PXRD to experimental data ----------
def normalizeToExpt(exptQ, exptInts, q, ints):
    
    intsMax, qAtIntsMax, maxIndex, intsMin = getNormVals(exptQ, exptInts)
    
    qRange = [qAtIntsMax-0.1, qAtIntsMax+0.1]
    
    normMax = 0
    for i in range(len(q)):
        if q[i] >= qRange[0] and q[i] <= qRange[1]:
            if ints[i] > normMax:
                normMax = ints[i]
    
    normInts = ints / normMax
            
    return normInts


# calculates a difference curve between two sets of PXRD data ----------
def diffCurve(q1, q2, ints1, ints2):
    '''
    Parameters
    ----------
    q1 (array_like) : dataset 1 Q values (Å-1)
    q2 (array_like) : dataset 2 Q values (Å-1)
    ints1 (array_like) : dataset 1 intensity values
    ints2 (array_like) : dataset 2 intensity values

    Returns
    -------
    diff_q (array_like) : Q values of difference curve
    diff_ints (array_like) : intensity values of difference curve
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

# calculates R^2 value between two sets of PXRD data ----------
def r2val(q1, q2, ints1, ints2):
    '''
    Parameters
    ----------
    q1 (array_like) : dataset 1 Q values (Å-1)
    q2 (array_like) : dataset 2 Q values (Å-1)
    ints1 (array_like) : dataset 1 intensity values
    ints2 (array_like) : dataset 2 intensity values

    Returns
    -------
    r2 (float) -- R^2 value
    '''
    
    import sklearn.metrics as skl
    
    q_list = []
    ints1_list = []
    ints2_list = []
    for i in range(len(q1)):
        q1_val = float('%.3f'%(q1[i]))
        for j in range(len(q2)):
            q2_val = float('%.3f'%(q2[j]))
            if q1_val == q2_val:
                q_list.append(q1_val)
                ints1_list.append(ints1[i])
                ints2_list.append(ints2[j])
    ints1_arr = np.array(ints1_list)
    ints2_arr = np.array(ints2_list)

    r2 = skl.r2_score(ints1_arr, ints2_arr)

    return r2

# calculates a difference curve and and R^2 value between two sets of PXRD data ----------
def diff_r2(q1, q2, ints1, ints2):
    import sklearn.metrics as skl
    
    q_list = []
    ints1_list = []
    ints2_list = []
    for i in range(len(q1)):
        q1_val = float('%.3f'%(q1[i]))
        for j in range(len(q2)):
            q2_val = float('%.3f'%(q2[j]))
            if q1_val == q2_val:
                q_list.append(q1_val)
                ints1_list.append(ints1[i])
                ints2_list.append(ints2[j])
    ints1_arr = np.array(ints1_list)
    ints2_arr = np.array(ints2_list)
    diff_ints = np.subtract(ints1_arr, ints2_arr)

    r2 = skl.r2_score(ints1_arr, ints2_arr)

    return r2, q_list, diff_ints

# calculates difference between two difference curves ----------
def fitDiff(diff_ints1, diff_ints2):
    '''
    Parameters
    ----------
    diff_ints1 (array_like) : dataset 1 difference in intensity values
    diff_ints2 (array_like) : dataset 2 difference in intensity values

    Returns
    -------
    fitDiff (array_like) : difference of differences intensities
    '''
    
    fitDiff = np.subtract(diff_ints1, diff_ints2)
    
    return fitDiff

# calculates R^2 values for each PXRD simulation in a directory against experimental data and generates a text file ----------
def simR2vals(simPath, exptPath, exptFN, exptWL, maxTT):
    '''
    Parameters
    ----------
    simPath (str) : file path of simulations directory
    exptPath (str) : file path of experimental data
    exptFN (str) : experimental data file name
    exptWL (str) : instrument wavelength
    maxTT (str) : maximum two theta

    Returns
    -------
    r2vals (array_like) : list of R2 values
    '''
    
    import pyfaults as pf
    import glob
    
    r2vals = []
    
    expt_q, expt_ints = pf.importExpt(exptPath, exptFN, exptWL, maxTT)
    
    sims = glob.glob('./simulations/*.txt')
    
    for f in sims:
        removeExt = f.split('.')
        fn = removeExt[0].split('\\')
        
        q, ints = pf.importFile('./simulations/', fn[-1])
        r2 = pf.analyze.r2val(expt_q, q, expt_ints, ints)
        
        r2vals.append([fn[-1], r2])
        
    with open('./r2vals.txt', 'w') as x:
        for i in range(len(r2vals)):
            for (fn, r2) in zip(r2vals[i][0], r2vals[i][1]):
                x.write('{0} {1}\n'.format(fn, r2))
        x.close()

    return r2vals
