import numpy as np
from scipy.spatial import distance_matrix
import random
import os, sys
import shutil

random.seed(42)

def path_selection(positions,n,resampling_num=10,distance_=1.62,distance_tol=0.3,closed=False):
    # Compute the distance matrix
    dist_matrix = distance_matrix(positions, positions)
    #print(dist_matrix)
    # It porduces the atomic number index [a,b] 
    connected_pairs = np.argwhere((dist_matrix > distance_-distance_tol) & (dist_matrix < distance_+distance_tol))
    if closed==False:
        loop_arry=np.zeros((1,n))
    else:
        loop_arry=np.zeros((1,n-1))
        
    for i in range(0,len(connected_pairs)):
        for s in range(0,resampling_num):
            current_loop_atom=[connected_pairs[i][0],connected_pairs[i][1]]
            counter=0
 
            initial_size=2
            
            while counter<n-initial_size:
                current_end_loc=current_loop_atom[-1]
                current_connected_pairs=connected_pairs[connected_pairs[:,0]==current_end_loc]
                keep_unique=~np.isin(current_connected_pairs[:,1], np.array(current_loop_atom))
                current_connected_pairs_unique=current_connected_pairs[keep_unique]
                len_remaining=len(current_connected_pairs)
                len_remaining_unique=len(current_connected_pairs_unique)
                last_equal_first = current_loop_atom[0] in current_connected_pairs[:,1]
                
                if len_remaining_unique>0 and counter<n-initial_size-1:
                    current_element=random.choice(current_connected_pairs_unique[:,1].tolist())
                    current_loop_atom.append(current_element)
                    counter+=1
                
                   

                elif last_equal_first and counter==n-initial_size-1 and closed==True: 
                    current_element=current_loop_atom[0]#random.choice(current_connected_pairs[:,1].tolist())
                    
                    #current_loop_atom.append(current_element)
                    counter+=1
                
                  
                elif counter==n-initial_size-1 and closed==False and len_remaining_unique>0:
                    current_element=random.choice(current_connected_pairs_unique[:,1].tolist())
                    current_loop_atom.append(current_element)
                    counter+=1
                   
                
                else:
                    #print('chain breaks. restart the sampling')
                   
                    break
            if counter==n-initial_size:
                 loop_arry=np.vstack([loop_arry,np.array(current_loop_atom)])
             
 
    return loop_arry
        

def position_process(classical_coord_path,unit_size,n,resampling_num=1,distance_=1.62,distance_tol=0.3,closed=False,parent_dir=0):
    
    # os.makedirs('./data/'+str(reference_number)+'/', exist_ok=True)
    # os.makedirs('./data/'+str(reference_number)+'/img', exist_ok=True)
    # shutil.copy(file_path_input, './data/')
    structure_coord_full= np.loadtxt(classical_coord_path)
    structure_coord=structure_coord_full[:,1:4]

    num_structures=int(len(structure_coord)/unit_size)
 
    for i in range(0,num_structures):
        print(str(i)+'/'+str(num_structures))
        current_coord=structure_coord[int(i*unit_size):int((i+1)*unit_size) ]
        loop_arry=path_selection(current_coord,n=n,resampling_num=resampling_num,distance_=distance_,distance_tol=distance_tol,closed=closed)
        #coord_file_path='data/'+str(reference_number)+'/'+str(i)+'.dat'
        coord_file_path=parent_dir+'/info/'+str(i)+'.dat'
        unique_loop_arry = np.unique(loop_arry, axis=0)
        np.savetxt(coord_file_path, unique_loop_arry, fmt='%s', delimiter='\t')


# Check if at least one command-line argument is provided
file_path_input='input_file.dat'
input_file=np.genfromtxt(file_path_input )[:,1]

n=int(input_file[0])

distance_=(input_file[1]) #si_o_bond_length
distance_tol=(input_file[2]) # bond tol 
resampling_num=int(input_file[3]) # sample n-member connected chains

unit_size=int(input_file[4])
reference_number=int(input_file[5])

classical_coord_path='sio.data'

if len(sys.argv) > 1:
    # The first command-line argument (sys.argv[1]) will be the input value
    
    parent_dir='./data/'+ sys.argv[1]+'/'

else:
    parent_dir='./data/'+str(int(reference_number))+'/'


# Create the directory
os.makedirs(parent_dir, exist_ok=True)
os.makedirs(parent_dir+'/img/', exist_ok=True)
os.makedirs(parent_dir+'/out/', exist_ok=True)
os.makedirs(parent_dir+'/info/', exist_ok=True)
shutil.copy(file_path_input, parent_dir)





 

closed=False

position_process(classical_coord_path,unit_size,n,resampling_num=resampling_num,distance_=distance_,distance_tol=distance_tol,closed=closed,parent_dir=parent_dir)
 