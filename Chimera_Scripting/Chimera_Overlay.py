import numpy as np
import argparse
import os
from matplotlib import cm

def attribute_to_residues(array,filename,attribute,description):
    ### Takes an array of values (per residue) and created a chimera-compatible attribute file.
    ### This file is used by chimera to assign color values to each residue.
    f = open(filename,"w+")
    f.write("#"+str(description)+"\n")
    f.write("attribute: "+str(attribute)+"\n")
    f.write("match mode: 1-to-1\n")
    f.write("recipient: residues\n")
    for i in range(len(array)):
        f.write("\t:"+str(i+1)+"\t"+str(array[i])+" \n")
    f.close()
    return

def fragment_map_to_hex(cmap,ticks): 
    ### converts a matplotlib colormap into an array of hexcodes for use in Chimera
    colormap = cm.get_cmap(cmap)
    color_range=np.linspace(0,1,ticks)
    hex_map=[]
    for i in range(ticks):
        r,g,b,a = colormap(color_range[i])
        r = int(r*255)
        g = int(g*255)
        b = int(b*255)
        a = int(a*255)
        hex_map.append('#%02x%02x%02x' % (r,g,b))
    return hex_map

def Chimera_attribute_color_string(hex_map,value_array):
    """
    Chimera colorbars are defined from top to bottom, hence the need for a reversed colorstring.
    Chimera is weird, I guess?
    """
    attr_color_array=[]
    colorbar_array=[]
    max_index=len(value_array)-1
    for i in range(len(value_array)):
        attr_color_array.append(str(value_array[i]))
        attr_color_array.append(hex_map[i])
        if i == max_index or i == 0 or value_array[i] == 0:
            colorbar_array.append(str(value_array[max_index-i]))
        else:
            colorbar_array.append("-")
        colorbar_array.append(hex_map[max_index-i])
    attr_color_string = " ".join(iter(attr_color_array))
    colorbar_string = " ".join(iter(colorbar_array))
    return attr_color_string,colorbar_string

def Chimera_Overlay_Image(pdbfile,array,image_base,number_of_rotations=6,make_movie=False,
                          colormap="viridis",hide_atom_mask=":1-10000",keep_files=True,
                          value_max=1,value_min=-1,quality=0):
    ### Adjust parameters based on quality variable, 0 is default for testing, 
    ### 1-3 for larger images, 2 and 3 use raytracing 
    if quality>3:
        quality=3
    windowsize = ["640 480","800 600","1024 768","2590 1920"]
    supersample = quality+1
    raytrace_true=["","","raytrace true","raytrace true"]
    image_raytrace=["","","raytrace rtwait rtclean","raytrace rtwait rtclean"]
    
    
    ### Generate the Chimera-compatible attribute file
    attribute_to_residues(array,f"{image_base}.dat","temp_att","temporary attribute")
    
    ### Get colormaps and convert to strings for use in Chimera command file.
    hex_mapping = fragment_map_to_hex(colormap,7)
    value_array = np.linspace(value_min,value_max,7)
    attr_color_string,colorbar_string = Chimera_attribute_color_string(hex_mapping,value_array)
    
    ### Build the Chimera command file
    chi_com_file=open(f"{image_base}.com","w+")
    chi_com_file.write(f"open {pdbfile}\n")                                               #open the pdb file.
    chi_com_file.write(f"windowsize {windowsize[quality]}\n")                               #sets windowsize.
    chi_com_file.write(f"background solid white\n")                                     #background to white.
    chi_com_file.write(f"~nucleotides\n")                                      # turn off nucleotide objects.
    chi_com_file.write(f"~display {hide_atom_mask}\n")                               #hide atoms in the mask.
    chi_com_file.write(f"defattr {image_base}.dat\n")         #define the attribute from the given data file.
    chi_com_file.write(f"rangecolor temp_att {attr_color_string}\n") #apply the color scale to the attribute.
    chi_com_file.write(f"colorkey 0.96,0.05 0.99,0.25 labelSide ")    #puts the colorbar at the bottom right.
    chi_com_file.write(f"left/top labelColor black justification right {colorbar_string}\n")
    
    ### Take rotational snapshots OR make a rotational movie.
    if make_movie == False: ### Do the rotation snapshots only if the user doesn't want a rotation movie.
        for i in range(number_of_rotations):
            chi_com_file.write(f"copy file {image_base}{i+1}.png png supersample ")
            chi_com_file.write(f"{supersample} {image_raytrace[quality]}\nturn y {360/number_of_rotations}\nwait\n")

    elif make_movie == True: ### Make a movie of the rotation.
        chi_com_file.write(f"movie record {raytrace_true[quality]}\n")
        chi_com_file.write("turn y 1 360\nwait\n")
        chi_com_file.write(f"movie encode framerate 30 output {image_base}.mp4\n")

    chi_com_file.write("stop\n")
    chi_com_file.close()
    
    ### Run chimera externally using the generated command file.
    os.system(f"chimera --gui {image_base}.com")
    if keep_files == False:
        os.system(f"rm {image_base}.dat {image_base}.com")

if __name__ == "__main__":
    print("This script is currently meant for use with other modules to generate chimera overlays.")
    print("It will be updated in the future to a callable program that can work with the different \noutputs from cpptraj.")
