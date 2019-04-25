'''
This script is designed to take two correlation matrix data files produced 
by cpptraj (part of the AmberTools suite) and provide multiple comparisons
between them.  Given the nature of how cpptraj represents correlated and
anti-correlated movement, it is necessary to consider the original values
as well as the resulting difference.

The code below produces four subplots in a single image file.  The top two 
subplots show the degree of change in correlated movement (left) or anti-
correlated movement (right).  The bottom plots show cases where a residue
pair has switched between correlated and anti-correlated movement, and the
magnitude of this change.

In its current form, the script requires modification by the user in three
fields: two input files to compare and an output filename to produce
'''


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
################
# File Loading #
################
# WT_Filename=input("Please enter the filename for the baseline (usually wildtype) including the full path:  \n>")
# MU_Filename=input("Please enter the filename for the comparison (usually mutant) including the full path:  \n>")
# outputfile=input("Please enter an output filename (ending in .png):  \n>")
WT_Filename="~/path/to/file/WT_correlation_mtrx.dat"
MU_Filename="~/path/to/file/MU_correlation_mtrx.dat"
outputfile='~/path/to/file/WT_MU_Correlation_Comparisons.png'
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

# Setting minimum and maximum ranges to ensure multiple data sets processed via this script are visually comparable.
Corr_Diff[Res_WT-1][Res_WT-1]=1.0
Anti_Diff[Res_WT-1][Res_WT-1]=1.0
Corr_to_Anti[Res_WT-1][Res_WT-1]=2.0
Anti_to_Corr[Res_WT-1][Res_WT-1]=2.0

Corr_Diff[Res_WT-2][Res_WT-1]=-1.0
Anti_Diff[Res_WT-2][Res_WT-1]=-1.0

###################
# Figure Plotting #
###################
plotted=plt.figure(1,figsize=(10,8),dpi=100)
corrplot=plt.subplot(2,2,1)
antiplot=plt.subplot(2,2,2)
corranti=plt.subplot(2,2,3)
anticorr=plt.subplot(2,2,4)

corr_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Corr_Diff.max() / (Corr_Diff.max() + abs(Corr_Diff.min())),1.0,'shiftedcorr')
corrplotax=corrplot.matshow(Corr_Diff,cmap=corr_shifted_cmap)
corrplot.set_title('Correlated Movement')
plotted.colorbar(corrplotax, ticks=[Corr_Diff.min(),0, Corr_Diff.max()],ax=corrplot)

anti_shifted_cmap=shiftedColorMap(cm.bwr_r,0.0,1 - Anti_Diff.max() / (Anti_Diff.max() + abs(Anti_Diff.min())),1.0,'shiftedcorr')
antiplotax=antiplot.matshow(Anti_Diff,cmap=anti_shifted_cmap)
antiplot.set_title('Anticorrelated Movement')
plotted.colorbar(antiplotax, ticks=[Anti_Diff.min(),0, Anti_Diff.max()],ax=antiplot)

corrantiax=corranti.matshow(Corr_to_Anti,cmap=cm.Greys)
corranti.set_title('Correlated to Anticorrelated')
plotted.colorbar(corrantiax, ticks=[Corr_to_Anti.min(),0, Corr_to_Anti.max()],ax=corranti)

anticorrax=anticorr.matshow(Anti_to_Corr,cmap=cm.Greys)
anticorr.set_title('Anticorrelated to Correlated')
plotted.colorbar(anticorrax, ticks=[0, Anti_to_Corr.max()],ax=anticorr)

plotted.tight_layout()
plotted.savefig(outputfile, dpi=1000,layout="tight")
