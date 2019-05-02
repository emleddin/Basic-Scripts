import numpy as np
import matplotlib.pyplot as plt

# Energy Values as a list into this variable.  
# The script will assume each value corresponds to a reaction coordinate in sequence

energies=[-52.3,-38.5,-40.1,-28.6,-43.4,-39.2,-58.4]
outputfile="output.png"

############################ Nothing below this line needs changing.

fig=plt.figure(figsize=(12,8),dpi=300)
ax=fig.add_subplot(111)
plt.plot(flat_linex,flat_liney,"k")
plt.plot(slope_x,slope_y,"k--")
for i in range(0,len(energies)-1):
    flat_linex=np.linspace(i+.75,i+1.25,10)
    flat_liney=np.full(10,energies[i])
    slope_x=np.linspace(i+1.25,i+1.75,10)
    slope_y=np.linspace(energies[i],energies[i+1],10)
    plt.plot(flat_linex,flat_liney,"k")
    plt.plot(slope_x,slope_y,"k--")
    ax.annotate(str(energies[i]),xy=(i+0.875,energies[i]-1.0))
    
flat_linex=np.linspace(len(energies)-.25,len(energies)+.25,10)
flat_liney=np.full(10,energies[len(energies)-1])
plt.plot(flat_linex,flat_liney,"k")
plt.plot(slope_x,slope_y,"k--")
ax.annotate(str(energies[len(energies)-1]),xy=(len(energies)-0.125,energies[len(energies)-1]-1.0))
plt.savefig(outputfile,layout="tight",dpi=300)
