##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np

#---------------------------------------------------------------------------------
# analyzes simulated PXRD data from PyFaults input file -------------------------
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
    save
        bool (optional) : saves text file with difference curve if True;
        default is False
    savePath
        str : file directory to save text file
    saveName
        str : file name to save text file

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

def r2val(q1, q2, ints1, ints2):
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

def fitDiff(diff_ints1, diff_ints2):
    fitDiff = np.subtract(diff_ints1, diff_ints2)
    
    return fitDiff

def simR2vals(simPath, exptPath, exptFN, exptWL, maxTT):
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