import numpy as np
import matplotlib.pyplot as plt

def cou_file_to_EDA_grid(filename):
    data=np.genfromtxt(filename, skip_header=1,delimiter=None,usecols=(1,2,3))
    rescount = int(data.max(axis=0)[1])
    cou_grid = np.zeros(shape=(rescount,rescount),dtype=float)
    for line in data:
        x = int(line[0]-1)
        y = int(line[1]-1)
        z = float(line[2])
        cou_grid[x][y] = z
    return cou_grid
    
def vdw_file_to_EDA_grid(filename):
    data=np.genfromtxt(filename,delimiter=None,usecols=(1,2,3))
    rescount = int(data.max(axis=0)[1])
    vdw_grid = np.zeros(shape=(rescount,rescount),dtype=float)
    for line in data:
        x = int(line[0]-1)
        y = int(line[1]-1)
        z = float(line[2])
        vdw_grid[x][y] = z
    return vdw_grid

def single_residue_EDA_plot(coulomb_file,vdw_file,residue_number,output_image_file):
    vdw = vdw_file_to_EDA_grid(vdw_file)
    coulomb = cou_file_to_EDA_grid(coulomb_file)
    resnum = residue_number-1
