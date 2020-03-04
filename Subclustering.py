#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
import sys

def Cluster_Separation(prmtopfile,trajectory,cluster_number_v_time):
    """
    Description:
        Takes existing clustering output from cpptraj and processes a single trajectory 
        into multiple subtrajectories that contain only frames corresponding to their respective
        clusters.  This may be useful as a means of narrowing down trajectories to a specific
        phase space required for a reaction to occur.
    Inputs:
        topology file (.prmtop)
        trajectory file (.mdcrd, .dcd, etc.)
        cluster data (.dat, should be frame and cluster columns from cpptraj)
    Outputs:
        cpptraj command textfile "subclustering.txt" for use with AmberTools cpptraj program.
        Job Submission script for use in PBS schedulers
    Returns:
        None
    """
    temp = pd.read_csv(cluster_number_v_time,delim_whitespace=True)
    if os.path.isfile("subclustering.txt"):
        os.system("rm subclustering.txt")
    f = open("subclustering.txt","w+")
    f.write("trajin "+trajectory+"\n")
    # f.close()
    for cluster in list(set(temp["c1"])):
        sub_df = temp[temp["c1"]==cluster]
    #     print(sub_df)
        subclust_framearray = np.asarray(sub_df["#Frame"])
    #     print(subclust_framearray)
        frame_list = ",".join(str(x) for x in subclust_framearray)
    #     print(frame_list)
        f.write("trajout Sub_Cluster_" + str(cluster) + ".mdcrd mdcrd onlyframes "+ frame_list + "\n")
    f.write("run\n")
    f.close()
    f=open("Cluster_Separation.sh","w+")
    f.write(f"""
#!/bin/bash
#PBS -q gac.cpu
#PBS -l nodes=1:ppn=20,mem=120GB
#PBS -j oe
#PBS -r n
#PBS -o Subclustering.logfile
#PBS -N Subclustering

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > PBS_NODEFILE
module load amber/19-serial
$AMBERHOME/bin/cpptraj -p {prmtopfile} < subclustering.txt > subclustering.log
""")
    f.close()
    return

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Expected 3 arguments:\n subcluster <prmtop> <trajectory> <clustnumvtime.dat>\n")
    else:
        Cluster_Separation(sys.argv[1],sys.argv[2],sys.argv[3])
