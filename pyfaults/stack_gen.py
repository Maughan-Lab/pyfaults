from diffpy.structure.structure import Atom
import pandas as pd
import os

def import_csv(directory, filename):
    path = os.path.join(directory, filename + ".csv")
    
    df = pd.read_csv(path)
    
    return df

def get_layer(df, layer_name, latt):
    import classes
    
    l = df[df["Layer"] == layer_name]
    
    atom_list = []
    for j in range(0, len(l.index)):
        atom = l.iloc[j]["Atom"]
        elem = l.iloc[j]["Element"]
        x = l.iloc[j]["x"]
        y = l.iloc[j]["y"]
        z = l.iloc[j]["z"]
        occ = l.iloc[j]["Occupancy"]
        
        pos = [x,y,z]
        
        new_atom = Atom(atype=elem, xyz=pos, label=atom, occupancy=occ, lattice=latt)
        
        atom_list.append(new_atom)
        
    new_layer = classes.Layer(atom_list, latt, layer_name)
    
    return new_layer

def gen_child(name, parent_layer, t_vect):
    import classes 
    child_layer = classes.Child_Layer(layer_name=name, parent=parent_layer, shift=t_vect)
    
    return child_layer

def gen_stack(*layers):
    import classes 
    unitcell = classes.Unitcell(layers)
    return unitcell