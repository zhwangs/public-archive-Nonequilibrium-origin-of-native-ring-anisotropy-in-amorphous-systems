import re
import numpy as np
# Read the file
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import os, sys

from scipy.optimize import minimize

def calculate_eigenvalues_and_eigenvectors(ring_arry):
    ring_center=np.mean(np.array(ring_arry),axis=0)
    diff_val=np.array(ring_arry)-ring_center
    # Calculate the covariance matrix
    covariance_matrix = np.cov(diff_val.T)
    
    # Perform eigenvalue decomposition
    eigenvalues, eigenvectors = np.linalg.eig(covariance_matrix)
    
    # Take the square root of eigenvalues to get the lengths for the eigenvectors
    eigenvalues = np.sqrt(eigenvalues)


    
    return ring_center, eigenvalues, eigenvectors
 



file_path_input='input_file.dat'
input_file=np.genfromtxt(file_path_input )[:,1]

n=int(input_file[0])

distance_=(input_file[1]) #si_o_bond_length
distance_tol=(input_file[2]) # bond tol 
resampling_num=int(input_file[3]) # sample n-member connected chains

unit_size=int(input_file[4])
reference_number=int(input_file[5])

classical_coord_path='sio.data'


unit_size=int(input_file[4])

tol_up=input_file[7]
tol_down=input_file[6]

bin_size=int(input_file[8])
length_end=int(input_file[9])


if len(sys.argv) > 1:
    # The first command-line argument (sys.argv[1]) will be the input value
    
    parent_dir='./data/'+ sys.argv[1]+'/'

else:
    parent_dir='./data/'+str(int(reference_number))+'/'


structure_coord_full= np.loadtxt(classical_coord_path)

structure_coord_name_arry=structure_coord_full[:,0]

structure_coord=structure_coord_full[:,1:4]
print('SS')
print(structure_coord)

num_structures=int(len(structure_coord)/unit_size)
print(num_structures)



#read file

 
print(structure_coord[int(1)])


roundness_tot=[]
roughness_tot=[]
for frame in range(0,num_structures):
    current_coord=structure_coord[int(frame*unit_size):int((frame+1)*unit_size) ]
    current_structure_coord_name_arry=structure_coord_name_arry[int(frame*unit_size):int((frame+1)*unit_size) ]
    sample_atom_loop= np.loadtxt(parent_dir+'/info/'+str(frame)+'_ring_index.dat') 
    current_len=len(sample_atom_loop)

    eigen_ring_info=np.zeros(3)
    roundness=[]
    roughness=[]

    for s in range(0,current_len):
        try:
            current_ring_row=sample_atom_loop[s]
            len(current_ring_row)
            
        except:
            print('single element only')
            current_ring_row=sample_atom_loop
            
        ring_arry=[]
        for i in range(0,len(current_ring_row)):
            current_index=int(current_ring_row[i])
            
            cell_arry=current_coord[int(current_index)]
            ring_arry.append(cell_arry)
        print(s)
        #print((current_ring_row))
        ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(ring_arry)
        sort_eigenvalues=np.sort(eigenvalues)
        roundness.append(sort_eigenvalues[-1]/sort_eigenvalues[-2])
        roughness.append(sort_eigenvalues[0]/np.sqrt(sort_eigenvalues[-1]*sort_eigenvalues[-2]))
        
        roundness_tot.append(sort_eigenvalues[-1]/sort_eigenvalues[-2])
        roughness_tot.append(sort_eigenvalues[0]/np.sqrt(sort_eigenvalues[-1]*sort_eigenvalues[-2]))

        eigen_ring_info=np.vstack([eigen_ring_info,ring_center,eigenvalues,eigenvectors.T])
    coord_file_path=parent_dir+'/info/'+str(frame)+'_eigen_ring_info.dat'
    np.savetxt(coord_file_path, eigen_ring_info, fmt='%s', delimiter='\t')

    roundness=np.array(roundness)
    roughness=np.array(roughness)
    roughness_roundness=np.vstack([roundness,roughness])
    # roundness first col 
    coord_file_path=parent_dir+'/info/'+str(frame)+'_round_rough.dat'
    np.savetxt(coord_file_path, roughness_roundness.T, fmt='%s', delimiter='\t')

roundness_tot=np.array(roundness_tot)
roughness_tot=np.array(roughness_tot)
roughness_roundness=np.vstack([roundness_tot,roughness_tot])
# roundness first col 
coord_file_path=parent_dir+'/info/'+'_total_round_rough.dat'
np.savetxt(coord_file_path, roughness_roundness.T, fmt='%s', delimiter='\t')
