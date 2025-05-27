import numpy as np
# Read the file
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize
import sys
import os
from scipy.spatial import distance_matrix

import shutil 

#This is for generating ideal structures with different foldings. Usually run 1 with the file input_file.dat
def optimal_rotation(source_points,target_points):
    rot, rssd, sens = R.align_vectors(target_points,source_points, return_sensitivity=True)

    transformed_points=rot.apply(source_points)

    return transformed_points


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

def alternate_repeat(string1, string2, times):
    result = []
    for i in range(times):
        if i % 2 == 0:
            result.append(string1)
        else:
            result.append(string2)
    return result


def pair_distance(x_1,x_2):
    return np.linalg.norm(x_1-x_2) 

def distance_(coord1,coord2):
    coord1_len=len(coord1)
    coord2_len=len(coord2)
    pair_full_arry=[]
    for i in range(0,coord1_len):
        current_coord1=coord1[i]
        current_coord2_right=coord2[i]
        current_coord2_left=coord2[i-1]

        pair_dis_arry=np.linalg.norm(current_coord2_right-current_coord1)
        pair_full_arry.append(pair_dis_arry)
        pair_dis_arry=np.linalg.norm(current_coord2_left-current_coord1)
        pair_full_arry.append(pair_dis_arry)

    return np.array(pair_full_arry) 


def recenter(coord):
    center_coord=np.mean(coord,axis=0)
    coord[:,0]=coord[:,0]-center_coord[0]
    coord[:,1]=coord[:,1]-center_coord[1]
    coord[:,2]=coord[:,2]-center_coord[2]

    return coord


def nearby_pair_distance(coord):
    # this is i+1 - i elements 
    coord=np.vstack([coord[-1],coord])
    
    diff = np.diff(coord,axis=0)
    
    return diff ,  np.linalg.norm(diff,axis=1) 
 


def get_element(coord,index=1):
    return coord[index::2]


def nearby_pair_distance(coord):
    # this is i+1 - i elements 
    coord=np.vstack([coord[-1],coord])
    
    diff = np.diff(coord,axis=0)
    
    return diff ,  np.linalg.norm(diff,axis=1) 

def pair_distance(x_1,x_2):
    return np.linalg.norm(x_1-x_2) 

def distance_(coord1,coord2):
    coord1_len=len(coord1)
    coord2_len=len(coord2)
    pair_full_arry=[]
    for i in range(0,coord1_len):
        current_coord1=coord1[i]
        current_coord2_right=coord2[i]
        current_coord2_left=coord2[i-1]

        pair_dis_arry=np.linalg.norm(current_coord2_right-current_coord1)
        pair_full_arry.append(pair_dis_arry)
        pair_dis_arry=np.linalg.norm(current_coord2_left-current_coord1)
        pair_full_arry.append(pair_dis_arry)

    return np.array(pair_full_arry) 

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

 
 
def random_rotation_n_element(vec,random_shift=3,l=0,mu=0, sigma=1,max_angle=np.pi/2,angle_arry=-1):
    #  axis of rotation with respect to the previous atom
    randomly_transformed_coord=np.zeros((n,3))
 
    vec=np.roll(vec, random_shift,axis=0)
 
    
    i=0
 
    check_number=0
    while i <n and check_number==0:
        if i==0:
            check_over_wind=n-l-1
 
        
        if i >check_over_wind:
            past_i_axis=i-l-1
            i=check_over_wind
            l=i-past_i_axis-1
            check_number=1

        current_vec_arry=[]
        vec_before=vec[i-l-1]
        vec_next=vec[i]
        #if l==1:
        random_axis_change = np.random.normal(mu, sigma, 3)

        # Random angle of rotation
        if angle_arry==-1:
            angle = np.random.uniform(-max_angle, max_angle)
        else:
            angle=angle_arry
        #angle=0.2*np.pi
        randomly_transformed_coord[i-l-1]=vec_before
        randomly_transformed_coord[i]=vec_next
 

        for s in range(i-l,i):

            vec_current=vec[s]
 
            centered_coord=vec_current-vec_before
            centered_coord=centered_coord#+ np.random.normal(mu, sigma, 3)
            axis= (vec_next-vec_before)#+random_axis_change
            axis=axis/np.linalg.norm(axis) 
                
            # Rodrigues' rotation formula
            cos_theta = np.cos(angle) 
            sin_theta = np.sin(angle)
            cross_mat = np.array([[0, -axis[2], axis[1]],
                                [axis[2], 0, -axis[0]],
                                [-axis[1], axis[0], 0]])
            rotation_matrix = cos_theta * np.eye(3) + sin_theta * cross_mat + (1 - cos_theta) * np.outer(axis, axis)


            transformed_coord = np.dot(centered_coord, rotation_matrix.T)+vec_before


            randomly_transformed_coord[s]=transformed_coord

            current_vec_arry.append(vec[s])
        #vec_current_list=vec[i-l:i:1]
        
        # before
        # print('before')
        # print(vec_before)
        # print('next')
        # print(vec_next)
        # print('list')
        # print(randomly_transformed_coord)

        

        i=i+l+1

        # print(i)
    randomly_transformed_coord=np.roll(randomly_transformed_coord, -random_shift,axis=0)
    return randomly_transformed_coord
#     
#     




def get_element(coord,index=1):
    return coord[index::2]

 

#     # transform the above centered coord

#    # transformed_coord=rotation_matrix@centered_coord.T#+vec_before
#     


#     return transformed_coord
 

def o_si_second_nn(coord):
    diff = np.diff(coord[::3],axis=0)
    return diff ,  np.linalg.norm(diff,axis=1) 



# input ------------------------------------------------------------------------------------------

global n 
global si_o_bond_length
global  sampling_number, max_angle,reference_number
global distance_tol, min_angle,sample_repeat,l,target_max

global seed_num_arry_si, seed_num_arry_o

 
file_path_input='input_file.dat'

input_file=np.genfromtxt(file_path_input )[:,1]



n=int(input_file[0])
si_o_bond_length=input_file[1]
string1 = 'Si'
string2 = 'O'

alternating_string = alternate_repeat(string1, string2, times=n)

sampling_number=int(input_file[2])
max_angle=np.pi*input_file[3]
min_angle=np.pi*input_file[4]

distance_tol=input_file[5]

sampling_repeat=int(input_file[6])

l=int(input_file[7])

random_shift=int(input_file[8])

reference_num=input_file[9]

max_l=n-2 # this is set to be maximum allowance 
#max_l=input_file[12] # usually set in input.
min_l=1#input_file[13]


target_min=input_file[22]
target_max=input_file[35]

only_odd=input_file[23]

sio_strength=int(input_file[41])
oo_strength=int(input_file[42])
sisi_strength=int(input_file[43])
si_o_second_nn_1_strength=int(input_file[44])
si_o_second_nn_2_strength=int(input_file[45])

angle_list=-1*np.ones(sampling_number)#np.linspace(min_angle,max_angle,sampling_number)
 


# Check if at least one command-line argument is provided
if len(sys.argv) > 1:
    # The first command-line argument (sys.argv[1]) will be the input value
    
    parent_dir='./data/'+ sys.argv[1]+'/'

else:
    parent_dir='./data/'+str(int(reference_num))+'/'


# Create the directory
os.makedirs(parent_dir, exist_ok=True)
os.makedirs(parent_dir+'/img/', exist_ok=True)
shutil.copy(file_path_input, parent_dir)

# end ------------------------------------------------------------------------------------------
 
def objective(coord):

    n=int(len(coord)/3)
    coord=coord.reshape((n,3))

    

    # controlling the strength between ring shape and line shape
    
    coord_si=get_element(coord,index=0)
    coord_O=get_element(coord,index=1)
    # Coordinates of the two atoms

    
    _, bond_length_si=nearby_pair_distance(coord_si)
    _, bond_length_o=nearby_pair_distance(coord_O)
    
    bond_length_o=bond_length_o[seed_num_arry_o]#np.random.choice(bond_length_o[0::], 5, replace=False)
    bond_length_si=bond_length_si[seed_num_arry_si]#np.random.choice(bond_length_si[0::], 6, replace=False)


    penalty_o_o= np.sum((((bond_length_o - 2.7) ) ) ** 2)   
    penalty_si_si=  np.sum((((bond_length_si - 3.14) ) ) ** 2)  
    

    SiO_distance=distance_(coord_si,coord_O)
 
    
    # Penalize deviation from the target bond length
    penalty_sio =np.sum((((SiO_distance - 1.64) ) ) ** 2) 

    _,si_o_second_nn=o_si_second_nn(coord)
 
    penalty_o_si_1= np.sum((((si_o_second_nn - 3.8) ) ) ** 2)   
    penalty_o_si_2= np.sum((((si_o_second_nn - 4.1) ) ) ** 2)   


    # Penalize deviation from the target bond length
    penalty = si_o_second_nn_1_strength*penalty_o_si_1+si_o_second_nn_2_strength*penalty_o_si_2+  sisi_strength*penalty_si_si+oo_strength*penalty_o_o+ sio_strength*penalty_sio 
 
    return penalty 

 
 

# ring structure
init_ring_coord, radius=generate_ring(n, si_o_bond_length=si_o_bond_length )
 
sample_atom_coord_arry_closed=init_ring_coord.copy()


# sampling
# calculate the pair difference for the ring
init_ring_coord_diff,init_ring_coord_diff_norm=nearby_pair_distance(init_ring_coord)
#init_ring_coord_diff_second,init_ring_coord_diff_norm_second=second_nearby_pair_distance(init_ring_coord)

roundness=[]
roughness=[]
init_diff=[]

current_sample=0
for i in range(0,sampling_number):

    # calculate the pair difference for the ring after random transformation 
    #randomly_transformed_coord=random_ring_transform(diff_vec=init_ring_coord_diff,vec=init_ring_coord,mu=mu, sigma=sigma,max_angle=max_angle)
    
    k=0
    while k==0:   
        randomly_transformed_coord=init_ring_coord 
        for s in range(0,sampling_repeat):
            # if want to specify a single l. 

            sample_shift=np.random.randint(0,n)
            
            
            # random_sd=np.random.randint(0,2*n)
            # wrapped_x = (random_sd ) % (max_l- min_l) 
            
            # sample_l=np.random.randint(min_l,min_l+(wrapped_x)+1)


            sample_l=np.random.randint(min_l,max_l+1)

            if only_odd==1:
                if sample_l % 2 == 0:
                    sample_l += 1
            else:
                pass
 
        
            current_angle=angle_list[i] 
            randomly_transformed_coord=random_rotation_n_element(vec=randomly_transformed_coord,random_shift=sample_shift,l=sample_l,mu=0, sigma=0,max_angle=max_angle,angle_arry=-1)

            # to check all distances are larger or equal to the bond length. it ensures that nothing is unphysical. 
        


        flatten_randomly_transformed_coord= randomly_transformed_coord.flatten()

    
        for w in range(0,1):
            # norm_tol=0.00001#np.random.normal(0,1, 1)
            # result = minimize(objective_sio,flatten_randomly_transformed_coord, method='CG',tol=norm_tol)
            # randomly_transformed_coord_flat=result.x
            

            rand_int=1#np.random.randint(1,2,1)
            seed_num_arry_si=np.random.choice(int(n/2),rand_int, replace=False)
            rand_int=1#np.random.randint(1,2,1)
            seed_num_arry_o=np.random.choice(int(n/2), rand_int, replace=False)

            sigma=0.3
            mu=0.3
            norm_tol=0.005#np.random.normal(mu,sigma , 1)
            if (sio_strength+oo_strength+sisi_strength+si_o_second_nn_1_strength+si_o_second_nn_2_strength)==0:
                print('no minimization')

                randomly_transformed_coord_flat=flatten_randomly_transformed_coord
            else:
                result = minimize(objective,flatten_randomly_transformed_coord, method='CG',tol=norm_tol)
                randomly_transformed_coord_flat=result.x






        randomly_transformed_coord=randomly_transformed_coord_flat.reshape((n,3))

 
        randomly_transformed_coord=optimal_rotation(source_points=randomly_transformed_coord,target_points=init_ring_coord )

        randomly_transformed_coord=recenter(randomly_transformed_coord) 


        coord_si=get_element(randomly_transformed_coord,index=0)
        coord_O=get_element(randomly_transformed_coord,index=1)
        SiO_distance=distance_(coord_si,coord_O)
   
        all_distance = distance_matrix(randomly_transformed_coord, randomly_transformed_coord)
      
        all_distance=all_distance.flatten()
   
        all_distance=all_distance[np.round(all_distance,5)>0]
        
        
        #distance_(randomly_transformed_coord,randomly_transformed_coord)
        if (np.min(all_distance)<target_min) or (np.min(all_distance)>target_max) :
            print('min distance is not within the range')
            print(np.min(all_distance))
            print('target range: ['+str(target_min)+' , '+str( target_max)+']')
            print('wdith: ')
            print( str(target_max-target_min))
            print('---------')
            print('# samples: '+str(current_sample))
            print(parent_dir)
            pass
        else:
            # this is for selection of a parituclar range of
            #lower_bound = 2.66
            #upper_bound = 2.73

            # Check if any value falls within the range
            #within_range = np.any((all_distance >= lower_bound) & (all_distance <= upper_bound))
            #values_in_range = all_distance[(all_distance >= lower_bound) & (all_distance <= upper_bound)]
            values_in_range=[1]
            if len(values_in_range)>0:
                
                k=1
                current_sample=current_sample+1
                print('---------')
                print('# samples: '+str(current_sample))
                print('within range distances')
                print(values_in_range)
            else:
                print('no distance within the provided range')
 
    
 
    optimize_coord=optimal_rotation(source_points=randomly_transformed_coord,target_points=init_ring_coord )

    optimize_coord=recenter(optimize_coord) 
    initial_diff_sum_square=np.sqrt(np.sum(np.sum((init_ring_coord-optimize_coord)**2,axis=1)))
     
    ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(optimize_coord)
    sort_eigenvalues=np.sort(eigenvalues)
    roundness.append(sort_eigenvalues[-1]/sort_eigenvalues[-2])
    roughness.append(sort_eigenvalues[0]/np.sqrt(sort_eigenvalues[-1]*sort_eigenvalues[-2]))
    init_diff.append(initial_diff_sum_square)
    sample_atom_coord_arry_closed=np.concatenate( (sample_atom_coord_arry_closed,optimize_coord),axis=0)
   
    new_row=[-1, sort_eigenvalues[-1]/sort_eigenvalues[-2],sort_eigenvalues[0]/np.sqrt(sort_eigenvalues[-1]*sort_eigenvalues[-2]),initial_diff_sum_square]
    coord_file_path=parent_dir+'round_rough.dat'
    with open(coord_file_path, 'a') as f:
        f.write(' '.join(map(str, new_row)) + '\n')
    # if i>2:
    #     roundness=np.array(roundness)
    #     roughness=np.array(roughness)
    #     init_diff=np.array(init_diff)
    #     angle_list=np.array(angle_list)
    #     roughness_roundness=np.vstack([angle_list,roundness,roughness,init_diff])

    #     coord_file_path=parent_dir+'round_rough.dat'
    #     with open(coord_file_path, 'a') as f:
    #         np.savetxt(coord_file_path, roughness_roundness.T, fmt='%s', delimiter='\t')

    #     roundness=roundness.tolist()
    #     roughness=roughness.tolist()
    #     init_diff=init_diff.tolist()
    #     angle_list=angle_list.tolist()
    
roundness=np.array(roundness)
roughness=np.array(roughness)
init_diff=np.array(init_diff)
angle_list=np.array(angle_list)
roughness_roundness=np.vstack([angle_list,roundness,roughness,init_diff])

# roundness first col 
print(parent_dir)
coord_file_path=parent_dir+'round_rough.dat'
np.savetxt(coord_file_path, roughness_roundness.T, fmt='%s', delimiter='\t')

coord_file_path=parent_dir+'coord.dat'
np.savetxt(coord_file_path, sample_atom_coord_arry_closed[n::], fmt='%s', delimiter='\t')
