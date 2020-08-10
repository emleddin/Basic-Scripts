#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys

def plot_PCA_from_NMD(filename,outputfile,title,plot_range=None):
    data = np.genfromtxt(filename,delimiter=None,skip_header=9)
    area = np.pi * (2)**2
    x = data[0,2] * data[0,3:]
    y = data[1,2] * data[1,3:]
    z = -(-x**2 - y**2)
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(1,1,1)
    ax.set_title(title,fontsize=18)
    ax.set_xlabel("Mode 1",fontsize=16)
    ax.set_ylabel("Mode 2",fontsize=16)
    if plot_range != None:
        ax.set_xlim(-plot_range,plot_range)
        ax.set_ylim(-plot_range,plot_range)
    ax.scatter(x,y,marker='o', s=area, zorder=10, alpha=0.4, c=z, edgecolors = 'black', cmap='viridis')
    plt.xticks(fontsize=14,rotation=90)
    plt.yticks(fontsize=14)
    plt.tight_layout()
    fig.savefig(outputfile,dpi=300)
    return


if __name__ == "__main__":
    if len(sys.argv) == 4:
        plot_PCA_from_NMD(sys.argv[1],sys.argv[2],sys.argv[3])
    elif len(sys.argv) == 5:
        plot_PCA_from_NMD(sys.argv[1],sys.argv[2],sys.argv[3],float(sys.argv[4]))
    else:
        print("Syntax Error:  Expected 3 or 4 arguments.\n\nplot_PCA <input file> <output file> <Image Title> [<Plot Range>]\n\n")
