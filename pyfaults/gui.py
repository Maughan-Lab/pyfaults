import PySimpleGUI as sg

sg.theme('DarkGrey13')

layout = [ [sg.Text('Supercell Generator')],
           [sg.Text('Number of Stacks (N)'), sg.InputText()],
           [sg.Text('Faulted Layer'), sg.InputText()],
           [sg.Text('Fault Probabilities'), sg.InputText()],
           [sg.Text('Stacking Vectors'), sg.InputText()]
           [sg.OK(), sg.Cancel()]]

window = sg.Window('Window Title', layout)

while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, 'Cancel'):
        break
    
window.close()


#supercells = pf.genSupercells.genSupercells(unitcell, nStacks, fltLayer, 
#                                            probList, sVecList, savePath)