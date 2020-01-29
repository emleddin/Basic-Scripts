##################
# Library Import #
##################
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy.linalg as la
import statsmodels.api as sm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm
from mpl_toolkits.axes_grid1 import AxesGrid
import sys

'''
This script is now callable with three arguments:  Reference file, Subject file, 
and output image filename.
It will process the correlation matrices given and produce a figure with 
four subplots showing the differences between systems.
'''




################################################################
# The shiftedColorMap function was taken without change from : #
# https://gist.github.com/phobson/7916777                      #
# Thanks to Paul Hobson!                                       #
################################################################

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def Correl_Compare(subject,reference,outputfile):
    WT_Filename=subject
    MU_Filename=reference
    temp=open(WT_Filename,"r")
    wildtype=temp.readlines()
    temp.close()
    temp=open(MU_Filename,"r")
    mutant=temp.readlines()
    temp.close()
    Res_WT=len(wildtype)+1
    Res_MU=len(mutant)+1

    #########################
    # Matrix Initialization #
    #########################
    Corr_WT=[]
    Corr_MU=[]
    Corr_Diff=np.array(np.zeros(Res_WT*Res_WT,float)).reshape(Res_WT,Res_WT)
    Anti_Diff=np.array(np.zeros(Res_WT*Res_WT,float)).reshape(Res_WT,Res_WT)
    Corr_to_Anti=np.array(np.zeros(Res_WT*Res_WT,float)).reshape(Res_WT,Res_WT)
    Anti_to_Corr=np.array(np.zeros(Res_WT*Res_WT,float)).reshape(Res_WT,Res_WT)
    for i in range(Res_WT-1):
        Corr_WT.append(wildtype[i].split())
    for i in range(Res_MU-1):
        Corr_MU.append(mutant[i].split())

    #####################
    # Matrix Processing #
    #####################
    for i in range(Res_WT-1):
        for j in range(Res_WT-1):
            if float(Corr_WT[i][j])>=0:
                if float(Corr_MU[i][j])>=0:
                    Corr_Diff[i][j]=float(Corr_MU[i][j])-float(Corr_WT[i][j])
                elif float(Corr_MU[i][j])<0:
                    Corr_to_Anti[i][j]=-float(Corr_MU[i][j])+float(Corr_WT[i][j])
            elif float(Corr_WT[i][j])<0:
                if float(Corr_MU[i][j])>=0:
                    Anti_to_Corr[i][j]=float(Corr_MU[i][j])-float(Corr_WT[i][j])
                elif float(Corr_MU[i][j])<0:
                    Anti_Diff[i][j]=-float(Corr_MU[i][j])+float(Corr_WT[i][j])

    Corr_Diff[Res_WT-1][Res_WT-1]=1.0
    Anti_Diff[Res_WT-1][Res_WT-1]=1.0
    Corr_to_Anti[Res_WT-1][Res_WT-1]=2.0
    Anti_to_Corr[Res_WT-1][Res_WT-1]=2.0

    Corr_Diff[Res_WT-2][Res_WT-1]=-1.0
    Anti_Diff[Res_WT-2][Res_WT-1]=-1.0
    # Corr_to_Anti[Res_WT-2][Res_WT-1]=-1.0
    # Anti_to_Corr[Res_WT-2][Res_WT-1]=-1.0

    # Corr_Diff=Corr_Diff/Corr_Diff.max()
    # Anti_Diff=Anti_Diff/Anti_Diff.max()
    # Corr_to_Anti=Corr_to_Anti/Corr_to_Anti.max()
    # Anti_to_Corr=Anti_to_Corr/Anti_to_Corr.max()

    la.norm(Anti_Diff)
    ###################
    # Figure Plotting #
    ###################
    plotted=plt.figure(figsize=(10,8),dpi=300)
    corrplot=plt.subplot(2,2,1)
    antiplot=plt.subplot(2,2,2)
    corranti=plt.subplot(2,2,3)
    anticorr=plt.subplot(2,2,4)

    corr_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Corr_Diff.max() / (Corr_Diff.max() + abs(Corr_Diff.min())),1.0,'shiftedcorr')
    corrplotax=corrplot.matshow(Corr_Diff,cmap=corr_shifted_cmap,origin="lower")
    corrplot.set_title('Correlated Movement')
    plotted.colorbar(corrplotax, ticks=[Corr_Diff.min(),0, Corr_Diff.max()],ax=corrplot)
    corrplot.xaxis.tick_bottom()

    anti_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Anti_Diff.max() / (Anti_Diff.max() + abs(Anti_Diff.min())),1.0,'shiftedcorr')
    antiplotax=antiplot.matshow(Anti_Diff,cmap=anti_shifted_cmap,origin="lower")
    antiplot.set_title('Anticorrelated Movement')
    plotted.colorbar(antiplotax, ticks=[Anti_Diff.min(),0, Anti_Diff.max()],ax=antiplot)
    antiplot.xaxis.tick_bottom()

    corrantiax=corranti.matshow(Corr_to_Anti,cmap=cm.Greys,origin="lower")
    corranti.set_title('Correlated to Anticorrelated')
    plotted.colorbar(corrantiax, ticks=[Corr_to_Anti.min(),0, Corr_to_Anti.max()],ax=corranti)
    corranti.xaxis.tick_bottom()

    anticorrax=anticorr.matshow(Anti_to_Corr,cmap=cm.Greys,origin="lower")
    anticorr.set_title('Anticorrelated to Correlated')
    plotted.colorbar(anticorrax, ticks=[0, Anti_to_Corr.max()],ax=anticorr)
    anticorr.xaxis.tick_bottom()

    plotted.tight_layout()
    plotted.savefig(outputfile, dpi=300,layout="tight")

    return

################
# File Loading #
################
if __name__ == "__main__":
    if len(sys.argv) == 4:
        Ref_File = sys.argv[1]
        Subj_File = sys.argv[2]
        Out_Image = sys.argv[3]
        Correl_Compare(Subj_File,Ref_File,Out_Image)
    else:
        print("Syntax Error:  Expected 3 arguments\n\npython <scriptname.py> <Reference File> <Subject File> <Output Image Filename>\n")
