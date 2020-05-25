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

def single_residue_EDA_row(vdw,cou,residue_number,buffer_space=2):
    rescount = len(vdw)
    resnum = residue_number-1
    front = resnum + buffer_space
    back = resnum - buffer_space
    courow = np.zeros(rescount,dtype=float)
    courow[:back] = cou[resnum,:back] + cou[:back,resnum]
    courow[front:] = cou[resnum,front:] + cou[front:,resnum]
    vdwrow = np.zeros(rescount,dtype=float)
    vdwrow[:back] = vdw[resnum,:back] + vdw[:back,resnum]
    vdwrow[front:] = vdw[resnum,front:] + vdw[front:,resnum]
    return vdwrow,courow

def plot_EDA(courow,vdwrow,output_file):
    totalrow = courow + vdwrow
    fig = plt.figure(figsize=(10,8),dpi=300)
    x = np.arange(1,len(vdwrow)+1,1)
    ax = fig.add_subplot(3,1,1)
    ax.set_title("Van der Waals Energy - Total: "+str(round(sum(vdwrow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(vdwrow)+1)
    ax.bar(x,vdwrow,align="center",color="green")
    ax = fig.add_subplot(3,1,2)
    ax.set_title("Coulomb Energy - Total: "+str(round(sum(courow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(courow)+1)
    ax.bar(x,courow,align="center",color="blue")
    ax = fig.add_subplot(3,1,3)
    ax.set_title("Combined Nonbonded Energy - Total: "+str(round(sum(totalrow),2))+" kcal/mol")
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_xlim(0,len(totalrow)+1)
    ax.bar(x,totalrow,align="center",color="red")
    fig.subplots_adjust(hspace=0.5)
    plt.savefig(output_file,dpi=300)
    return None

def res_range_parser(res_list_string):
    res_list_commas = res_list_string.split(',')
    full_res_list = []
    for i in res_list_commas:
        if '-' in i:
            start = int(i.split('-')[0])
            end = int(i.split('-')[1])
            for j in np.arange(start,end+1,1):
                full_res_list.append(j)
        else:
            full_res_list.append(int(i))
    print(full_res_list)
    return full_res_list

if __name__ == "__main__":
    if "help" in sys.argv:
        print("Help Response Coming Soon...")
    elif len(sys.argv) < 5:
        print("Syntax Error:  Expected at least 4 arguments.\n\nplot_single_residue_EDA <coulomb filename> <vdw filename> <output filename> <target residues> [buffer]\n---Target Residues should be a list separated by commas, or ranges such as 13-25, without spaces.\n\n")
    elif len(sys.argv) == 5:
        cou = cou_file_to_EDA_grid(sys.argv[1])
        vdw = vdw_file_to_EDA_grid(sys.argv[2])
        output = sys.argv[3]
        reslist = res_range_parser(sys.argv[4])
        full_vdw, full_cou = single_residue_EDA_row(vdw,cou,reslist[0])
        for i in reslist[1:]:
            vdwrow,courow = single_residue_EDA_row(vdw,cou,i)
            full_vdw = full_vdw + vdwrow
            full_cou = full_cou + courow
        plot_EDA(courow,vdwrow,output)
    elif len(sys.argv) == 6:
        cou = cou_file_to_EDA_grid(sys.argv[1])
        vdw = vdw_file_to_EDA_grid(sys.argv[2])
        output = sys.argv[3]
        buffer_space = int(sys.argv[5])
        reslist = res_range_parser(sys.argv[4])
        full_vdw, full_cou = single_residue_EDA_row(vdw,cou,reslist[0],buffer_space)
        for i in reslist[1:]:
            vdwrow,courow = single_residue_EDA_row(vdw,cou,i,buffer_space)
            full_vdw = full_vdw + vdwrow
            full_cou = full_cou + courow
        plot_EDA(courow,vdwrow,output)

