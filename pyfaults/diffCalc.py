##################################################################################
# Author: Sinclair R. Combs
##################################################################################

import numpy as np
import copy as cp
import os

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
    
#---------------------------------------------------------------------------------
# generates experimental vs simulated XRD difference curves ----------------------
#---------------------------------------------------------------------------------
def calcDiffs(path, simDF, expt):
    '''
    Parameters
    ----------
    path
        str : file save directory
    simDF
        DataFrame : tabulated data returned from calcSims()
    expt
        list (nparray) : experimental data as [Q, ints]

    Returns
    -------
    exptDiffDF
        DataFrame : tabulated data, includes the following
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
            'Simulated Q' -- calculated Q values (A^-1)
            'Simulated Intensity' -- calculated intensity values
            'Expt vs. Model Q' -- difference curve Q Values (A^-1)
            'Expt vs. Model Difference' -- difference curve intensity values
    '''
    # create 'diffCurves' folder in file directory
    if os.path.exists(path + 'diffCurves/') == False:
        os.mkdir(path + 'diffCurves/')
    
    # create copy of simDF dataframe
    exptDiffDF = cp.deepcopy(simDF)
    exptQ = expt[0]
    exptInts = expt[1]
    
    exptDiffQ = []
    exptDiffInts = []
    # calculate experimental vs model difference for each supercell
    for row in simDF.index:
        name = simDF['Model'][row]
        simQ = simDF['Simulated Q'][row]
        simInts = simDF['Simulated Intensity'][row]
        diffQ, diffInts = diffCurve(exptQ, simQ, exptInts, simInts)
        exptDiffQ.append(diffQ)
        exptDiffInts.append(diffInts)
        # save difference data
        with open(path + 'diffCurves/' + name + '_exptDiff.txt', 'w') as f:
            for (q, ints) in zip(diffQ, diffInts):
                f.write('{0} {1}\n'.format(q, ints))
        f.close() 
    
    # append each entry with difference data
    exptDiffDF['Expt vs. Model Q'] = exptDiffQ
    exptDiffDF['Expt vs. Model Difference'] = exptDiffInts
    return exptDiffDF

#---------------------------------------------------------------------------------
# calculates fit difference curves -----------------------------------------------
# (diff, expt vs unfaulted) - (diff, expt vs faulted) ----------------------------
#---------------------------------------------------------------------------------
def calcFitDiffs(path, exptDiffDF):
    '''
    Parameters
    ----------
    path
        str : file save directory
    exptDiffDF
        DataFrame : tabulated data returned from calcDiffs()

    Returns
    -------
    fitDiffDF
        DataFrame : tabulated data, includes the following
            'Model' -- unique identifier
            'Stacking Vector' -- stacking vector in string format
            'Stacking Probability' -- stacking fault probability in str format
            'S_x' -- x-component of stacking vector
            'S_y' -- y-component of stacking vector
            'S_z' -- z-component of stacking vector
            'P' -- fault probability (0 to 1)
            'Simulated Q' -- calculated Q values (A^-1)
            'Simulated Intensity' -- calculated intensity values
            'Expt vs. Model Q' -- difference curve Q Values (A^-1)
            'Expt vs. Model Difference' -- difference curve intensity values
            'UF vs. FLT Model' -- intensity values of fit difference curve
    '''
    # create 'fitDiffCurves' folder in file directory
    if os.path.exists(path + 'fitDiffCurves/') == False:
        os.mkdir(path + 'fitDiffCurves/')
    
    # create copy of exptDiffDF dataframe
    fitDiffDF = cp.deepcopy(exptDiffDF)
    
    # set unfaulted model
    UFdiff = exptDiffDF['Expt vs. Model Difference'][0]
    
    fitDiffs = []
    # calculate fit differences
    for i in exptDiffDF.index:
        if i == 0:
            fitDiffs.append(np.nan)
        elif i > 0:
            name = exptDiffDF['Model'][i]
            modelDiff = exptDiffDF['Expt vs. Model Difference'][i]
            diff = np.subtract(UFdiff, modelDiff)
            fitDiffs.append(diff)
            # save fit difference data
            with open(path + 'fitDiffCurves/' + name + '_fitDiff.txt', 'w') as f:
                for diff in range(len(fitDiffs)):
                    f.write('{0}\n'.format(diff))
            f.close() 
    
    # append each entry with fit difference data
    fitDiffDF['UF vs. FLT Model'] = fitDiffs
    return fitDiffDF