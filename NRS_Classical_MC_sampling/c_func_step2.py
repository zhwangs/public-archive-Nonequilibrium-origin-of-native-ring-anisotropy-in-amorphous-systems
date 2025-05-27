import re
import numpy as np
# Read the file
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.spatial import distance_matrix
import sys
from scipy.optimize import minimize

def hist_info(data,bins=100):
    hist, bin_edges = np.histogram(data, bins=bins)
    return hist, bin_edges



def c_function(atom_coord_arry):
    data_size=len(atom_coord_arry)
    
    c_arry=np.array([])
    for current_index in range(1,data_size):
        current_val = np.tile(atom_coord_arry[0], (int(data_size-current_index), 1))
        atom_coord_arry.pop(0)
        diff_arry=np.array(atom_coord_arry)-current_val
        average_disp_vec=np.sqrt((diff_arry[:,0] )**2+(diff_arry[:,1] )**2+(diff_arry[:,2])**2)
        c_arry=np.concatenate((c_arry, average_disp_vec))

    return c_arry
    
# input ------------------------------------------------------------------------------------------




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

total_element=int(len(np.loadtxt(classical_coord_path))/int(input_file[4]))


bin_edges = np.linspace(0,length_end , int(bin_size+1))
# print('Delta length')
# print(bin_edges[1]-bin_edges[0])
bin_edges_true = (bin_edges[1:] + bin_edges[:-1]) / 2

if len(sys.argv) > 1:
    # The first command-line argument (sys.argv[1]) will be the input value
    
    parent_dir='./data/'+ sys.argv[1]+'/'

else:
    parent_dir='./data/'+str(int(reference_number))+'/'


# input end ------------------------------------------------------------------------------------------
hist_arry=bin_edges_true


structure_coord_full= np.loadtxt(classical_coord_path)
structure_coord_name_arry=structure_coord_full[:,0]
 
structure_coord=structure_coord_full[:,1:4]
 
num_structures=int(len(structure_coord)/unit_size)
# print(num_structures)

total_num_ct=[]
c_arry_sisi=0
c_arry_sio=0
c_arry_oo=0
for i in range(0,num_structures):
    current_coord=structure_coord[int(i*unit_size):int((i+1)*unit_size) ]
    current_name=structure_coord_name_arry[int(i*unit_size):int((i+1)*unit_size) ]
    #read file
    sample_atom_loop= np.loadtxt(parent_dir+'/info/'+str(i)+'.dat') 
    total_loop_num=len(sample_atom_loop)
    current_ring_row_num=len(sample_atom_loop[0])
    # print('SSS')
    # print(current_ring_row_num)
    unique_num_c_arry=int((current_ring_row_num**2-current_ring_row_num)/2)
    # print(unique_num_c_arry)
    c_arry=np.zeros((1,unique_num_c_arry))


    # print(total_loop_num)
     
     # save coords for those rings
    ring_index_tot=[] 
    dis_arry_distribution=[]
    current_num_count=0
    for j in range(1,total_loop_num):
        
        current_ring_row=sample_atom_loop[j]

        init_index=int(current_ring_row[0])
        end_index=int(current_ring_row[-1])

        dis_arry=np.linalg.norm(current_coord[init_index]-current_coord[end_index]) 
        dis_arry_distribution.append(dis_arry)
        check_element=(dis_arry < tol_up) & (dis_arry > tol_down)
        if check_element:
            # print(current_ring_row)
            # print('S')
            # print(dis_arry)
            ring_arry=[]
            ring_name_arry=[]
            for k in range(0,len(current_ring_row)):
                current_index=int(current_ring_row[k])
                #print(current_coord[int(current_index)])
                ring_arry.append(current_coord[int(current_index)])
                ring_name_arry.append(current_name[int(current_index)])
            ring_arry=np.array(ring_arry)
            ring_name_arry=np.array(ring_name_arry)
            # sio numbers
            si_coord=ring_arry[ring_name_arry==1]
            o_coord=ring_arry[ring_name_arry==2]


            c_matrix_sisi_flat=c_function((si_coord).tolist())
            hist_count_sisi,_=hist_info(np.array(c_matrix_sisi_flat),bins=bin_edges)
 

            c_arry_sisi=c_arry_sisi+(hist_count_sisi)
 

            c_matrix_oo_flat=c_function((o_coord).tolist())
            hist_count_oo,_=hist_info(np.array(c_matrix_oo_flat),bins=bin_edges)
 
            c_arry_oo=c_arry_oo+(hist_count_oo)
 
            c_matrix_sio=distance_matrix(si_coord,o_coord)
            #print(c_matrix_sio)
            c_matrix_sio_flat=c_matrix_sio.flatten()
            hist_count_sio,_=hist_info(np.array(c_matrix_sio_flat),bins=bin_edges)
            c_arry_sio=c_arry_sio+(hist_count_sio)
 

            c_matrix_flat=c_function((ring_arry).tolist())
            #print(c_matrix_flat)
            #c_matrix_flat=np.unique(c_matrix_.flatten())
            ring_index_tot.append(current_ring_row)

            c_arry=np.vstack([c_arry,c_matrix_flat])
            current_num_count=current_num_count+1
        else:
            
            pass
        if len(c_arry)<2:
            print('tol too tight')
    total_num_ct.append(current_num_count)
    ring_index_tot=np.array(ring_index_tot)
    #print(ring_index_tot.shape)
    # save the rings 
    coord_file_path=parent_dir+'/info/'+str(i)+'_ring_index.dat'
    np.savetxt(coord_file_path, ring_index_tot, fmt='%s', delimiter='\t')


    coord_file_path=parent_dir+'/info/'+str(i)+'_carry.dat'
    np.savetxt(coord_file_path, np.array(c_arry)[1::,0::], fmt='%s', delimiter='\t')

    hist_count,_=hist_info(np.array(dis_arry_distribution),bins=bin_edges)
    hist_arry=np.vstack([hist_arry,hist_count])
    
coord_file_path=parent_dir+'/info/'+'_structure_index_vs_ring_counts.dat'
np.savetxt(coord_file_path, np.array(total_num_ct), fmt='%s', delimiter='\t')


coord_file_path=parent_dir+'/info/'+'_dis_arry_distribution.dat'
np.savetxt(coord_file_path, hist_arry, fmt='%s', delimiter='\t')



coord_file_path=parent_dir+'/info/'+'hist_sio.dat'
np.savetxt(coord_file_path, np.array(c_arry_sio), fmt='%s', delimiter='\t')
coord_file_path=parent_dir+'/info/'+'hist_oo.dat'
np.savetxt(coord_file_path, np.array(c_arry_oo), fmt='%s', delimiter='\t')
coord_file_path=parent_dir+'/info/'+'hist_sisi.dat'
np.savetxt(coord_file_path, np.array(c_arry_sisi), fmt='%s', delimiter='\t')
