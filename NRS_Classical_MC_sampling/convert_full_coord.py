import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import random

from matplotlib.animation import FuncAnimation
import os, sys
from scipy.spatial.transform import Rotation as R

def recenter(coord):
    center_coord=np.mean(coord,axis=0)
    coord[:,0]=coord[:,0]-center_coord[0]
    coord[:,1]=coord[:,1]-center_coord[1]
    coord[:,2]=coord[:,2]-center_coord[2]

    return coord



def generate_ring(n, si_o_bond_length=1.62):
    # Calculate the radius of the ring
    radius = si_o_bond_length / (2 *np.sin(np.pi / n))
    
    # Initialize lists to store coordinates
    init_ring_coord=np.zeros((n,3))

    # Generate coordinates for each atom
    for i in range(n):
        angle = 2 * np.pi * i / n
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        init_ring_coord[i][0]=x
        init_ring_coord[i][1]=y

    return init_ring_coord, radius

def optimal_rotation(source_points,target_points):
 
    center_zero_coord=recenter(np.array(source_points))
    length_distance=np.linalg.norm(center_zero_coord,axis=1)
    print('ssss')
    print(center_zero_coord)
    print(length_distance)
    rot, rssd, sens = R.align_vectors(target_points,source_points, return_sensitivity=True)

    transformed_points=rot.apply(source_points)
    center_zero_coord=recenter(transformed_points)
    length_distance=np.linalg.norm(center_zero_coord,axis=1)
    print('ssss2')
    print(transformed_points)
    print(length_distance)

    return transformed_points


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


replace=True
file_round_rough=parent_dir+'out/round_rough_real.dat'
file_coord=parent_dir+'out/coord_real.dat'
file_coord_index=parent_dir+'out/index_coord_real.dat'

if replace:

    with open(file_round_rough, 'w') as file:
        pass
    with open(file_coord, 'w') as file:
        pass
    with open(file_coord_index, 'w') as file:
        pass



init_ring_coord, radius=generate_ring(n, si_o_bond_length=1.62 )
 



structure_coord_full= np.loadtxt(classical_coord_path)

structure_coord_name_arry=structure_coord_full[:,0]

structure_coord=structure_coord_full[:,1:4]

# print('SS')
# print(structure_coord)

num_structures=int(len(structure_coord)/unit_size)
print(num_structures)

 

def update(frame):

    current_coord=structure_coord[int(frame*unit_size):int((frame+1)*unit_size) ]
    current_structure_coord_name_arry=structure_coord_name_arry[int(frame*unit_size):int((frame+1)*unit_size) ]
    sample_atom_loop= np.loadtxt(parent_dir+'/info/'+str(frame)+'_ring_index.dat')
    cov_eigen_info=np.loadtxt(parent_dir+'/info/'+str(frame)+'_eigen_ring_info.dat')

    try:
        
        record_arry=[0]
        for s in range(0,10000):
            #s=random.randint(0,len(sample_atom_loop))
            current_ring_row=sample_atom_loop[s]
            # it is shifted by 0 0 0 in the file 
            current_eigen_center_coord=cov_eigen_info[5*s+1,:]
            current_eigen_eigenvalue=cov_eigen_info[5*s+2,:]
            current_eigen_vectors=cov_eigen_info[5*s+3:5*s+6,:]
            # print('value')
            # print(current_eigen_eigenvalue)

            sort_eigenvalues=np.sort(current_eigen_eigenvalue)


            current_ring_arry=[]
            current_ring_arry_index=[]
            current_rough_round=[]
            for i in range(0,len(current_ring_row)):
                current_index=int(current_ring_row[i])
                cell_arry=current_coord[int(current_index)]
                current_ring_arry_index.append(current_structure_coord_name_arry[int(current_index)])
                current_ring_arry.append(cell_arry)

     
            current_length=np.sum(np.abs(current_ring_arry))
   
            if len(record_arry)==1:
                record_arry.append(current_length)
                check_record=1000
            else:
                check_record= np.min(np.abs(np.array(record_arry)-current_length))

 
            if check_record>0.1:
 
                record_arry.append(current_length)
                
                # this rotation will keep the index order of the initial vectors. 
                optimize_coord=optimal_rotation(source_points=current_ring_arry,target_points=init_ring_coord )
                
                optimize_coord=recenter(optimize_coord) 
                
                #optimize_coord=recenter(np.array(current_ring_arry) )

                initial_diff_sum_square=np.sqrt(np.sum(np.sum((init_ring_coord-optimize_coord)**2,axis=1)))

                # roll back array if the first atom is 1. So all rings start with index 2. 
                

                if current_ring_arry_index[0]==1:
                    current_ring_arry_index=np.roll(current_ring_arry_index,-1,axis=0)

                    print(optimize_coord)
                    optimize_coord=np.roll(optimize_coord,-1,axis=0)
                    print(optimize_coord)
                else:
                    pass
                for j in range(0,len(optimize_coord)):
                    with open(file_coord, 'a') as file:
                        # Iterate through each row of the NumPy array
                        # Convert each element to string and join with a space
                        row_str = ' '.join(map(str, optimize_coord[j]))
                        # Write the row to the file
                        file.write(row_str + '\n')


                with open(file_coord_index, 'a') as file:
                    # Iterate through each row of the NumPy array
                    # Convert each element to string and join with a space
                    # Convert each element to string and join with a space
                    row_str = ' '.join(map(str,  np.array(current_ring_arry_index)))
                    # Write the row to the file
                    file.write(row_str + '\n')
                    


                current_rough_round.append(-1)
                current_rough_round.append(sort_eigenvalues[-1]/sort_eigenvalues[-2])
                current_rough_round.append(sort_eigenvalues[0]/np.sqrt(sort_eigenvalues[-1]*sort_eigenvalues[-2]))
                current_rough_round.append(initial_diff_sum_square)
                with open(file_round_rough, 'a') as file:
                    # Iterate through each row of the NumPy array
                    # Convert each element to string and join with a space
                    row_str = ' '.join(map(str, current_rough_round))
                    # Write the row to the file
                    file.write(row_str + '\n')
                
            else:
                pass
    except:
        print('passsss')
        pass
 
num_structures=2
starting_index=0
for i in range(starting_index,num_structures):
    update(i)

# if want to see the images
 
