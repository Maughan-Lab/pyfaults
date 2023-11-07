from pyfaults.diffCurve import diffCurve, saveDiffCurve
from pyfaults import importExpt, importFile

wl = 0.459744

path = 'C:/Users/sincl/Documents/MaughanLab/pyfaults/pyfaults/test/testoutput/'
save = 'C:/Users/sincl/Documents/MaughanLab/pyfaults/pyfaults/test/testdata/'

expt_q, expt_ints = importExpt(save, 'expt_XRD', wl, 30)
UF_q, UF_ints = importFile(save, 'Unfaulted_sim')
S1_P10_q, S1_P10_ints = importFile(save, 'S1_P10_sim')
S2_P10_q, S2_P10_ints = importFile(save, 'S2_P10_sim')
S1_P20_q, S1_P20_ints = importFile(save, 'S1_P20_sim')
S2_P20_q, S2_P20_ints = importFile(save, 'S2_P20_sim')

UF_dq, UF_dints = diffCurve(expt_q, UF_q, expt_ints, UF_ints)
S1_P10_dq, S1_P10_dints = diffCurve(expt_q, S1_P10_q, expt_ints, S1_P10_ints)
S2_P10_dq, S2_P10_dints = diffCurve(expt_q, S2_P10_q, expt_ints, S2_P10_ints)
S1_P20_dq, S1_P20_dints = diffCurve(expt_q, S1_P20_q, expt_ints, S1_P20_ints)
S2_P20_dq, S2_P20_dints = diffCurve(expt_q, S2_P20_q, expt_ints, S2_P20_ints)

saveDiffCurve(UF_dq, UF_dints, save, 'UF_diff')
saveDiffCurve(S1_P10_dq, S1_P10_dints, save, 'S1_P10_diff')
saveDiffCurve(S2_P10_dq, S2_P10_dints, save, 'S2_P10_diff')
saveDiffCurve(S1_P20_dq, S1_P20_dints, save, 'S1_P20_diff')
saveDiffCurve(S2_P20_dq, S2_P20_dints, save, 'S2_P20_diff')