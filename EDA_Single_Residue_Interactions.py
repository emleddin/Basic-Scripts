#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import sys

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

def single_residue_EDA_row(vdw,cou,residue_number,buffer_space=0):
    resnum = residue_number-1
    front = resnum + buffer_space
    back = resnum - buffer_space
    courow = np.zeros(rescount,dtype=float)
    courow[:back] = cou[0,:back] + cou[:back,-1]
    courow[front:] = cou[0,front:] + cou[front:,-1]
    vdwrow = np.zeros(rescount,dtype=float)
    vdwrow[:back] = vdw[0,:back] + vdw[:back,-1]
    vdwrow[front:] = vdw[0,front:] + vdw[front:,-1]
    return vdwrow,courow

def plot_EDA(courow,vdwrow,output_file):
    totalrow = courow + vdwrow
    fig = plt.figure(figsize=(10,8),dpi=300)
    x = np.arange(1,len(vdwrow)+1,1)
    ax = fig.add_subplot(1,3,1)
    ax.set_title("Van der Waals Energy - Total: "+str(round(sum(vdwrow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(vdwrow)+1)
    ax.bar(x,vdwrow,align="center",color="green")
    ax = fig.add_subplot(1,3,2)
    ax.set_title("Coulomb Energy - Total: "+str(round(sum(courow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(courow)+1)
    ax.bar(x,courow,align="center",color="green")
    ax = fig.add_subplot(1,3,3)
    ax.set_title("Combined Nonbonded Energy - Total: "+str(round(sum(totalrow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(totalrow)+1)
    ax.bar(x,totalrow,align="center",color="green")
    plt.savefig(output_file,dpi=300)
    return None

if __name__ == "__main__":
    if "help" in sys.argv:
        pass
    elif len(sys.argv) < 5:
        print("Syntax Error:  Expected at least 4 arguments.\n\nplot_single_residue_EDA <coulomb filename> <vdw filename> <output filename> <target residue> [buffer]\n\n"
    elif len(sys.argv) == 5:
        cou = cou_file_to_EDA_grid(sys.argv[1])
        vdw = vdw_file_to_EDA_grid(sys.argv[2])
        output = sys.argv[3]
        resnum = sys.argv[4]
        vdwrow,courow = single_residue_EDA_row(vdw,cou,resnum)
        plot_EDA(courow,vdwrow,output)
    elif len(sys.argv) == 6:
        cou = cou_file_to_EDA_grid(sys.argv[1])
        vdw = vdw_file_to_EDA_grid(sys.argv[2])
        output = sys.argv[3]
        resnum = sys.argv[4]
        buffer_space = sys.argv[5]
        vdwrow,courow = single_residue_EDA_row(vdw,cou,resnum,buffer_space)
        plot_EDA(courow,vdwrow,output)
