##################
# Library Import #
##################
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import cm

##########################################################
# File Loading - Filename assumes path beginning at $PWD #
##########################################################
inputfile=input("Please enter the path and filename for input: ")
outputfile=input("Please enter the path and filename for output image: ")
temp=open(inputfile,"r")
filedata=temp.readlines()
temp.close()

###################################################
# Adjust Residue Count, establish pairwise matrix #
###################################################
Res_Count=182
matrix = [[1.]*(Res_Count+1) for _ in range(Res_Count+1)]
highest=0

##################################################
# Process File, summing each pairwise            #
# interaction (regardless of atom in residue)    #
##################################################
for i in range(1,len(filedata)):
    acceptor=filedata[i].split()[0].split('_')[1].split('@')[0]
    donor=filedata[i].split()[1].split('_')[1].split('@')[0]
    framecount=int(filedata[i].split()[3])
    matrix[int(donor)][int(acceptor)]+=framecount
    if (matrix[int(donor)][int(acceptor)]>highest):
        highest=matrix[int(donor)][int(acceptor)]
        
######################################################
# Adjust to logarithmic scale and normalize from 0-1 #
######################################################
normalization=1/np.log(highest)
for i in range(Res_Count+1):
    for j in range(Res_Count+1):
        matrix[i][j]=normalization*float(np.log(matrix[i][j]))

############################
# Plot the pairwise matrix #
############################
fig = plt.figure(1,figsize=(10,8),dpi=100)
ax = fig.add_subplot(111)
cax=ax.matshow(matrix,cmap=cm.Blues_r)
plt.ylabel("Donor Residue", fontsize=12)
plt.xlabel("Acceptor Residue", fontsize=12)
ax.xaxis.set_label_position('top')
fig.colorbar(cax, ticks=[0, 1])
plt.savefig(outputfile, dpi=1000,layout="tight")
