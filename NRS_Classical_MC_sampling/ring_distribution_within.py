import numpy as np
from scipy.spatial import distance_matrix
import networkx as nx
import random
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import os, sys,shutil
from scipy.stats import skewnorm

def hist_info(data,bins=100):
    hist, bin_edges = np.histogram(data, bins=bins)
    return hist, bin_edges
def c_function_2(atom_coord_arry_1,atom_coord_arry_2,smear=0.15):
    c_matrix_=distance_matrix(atom_coord_arry_1,atom_coord_arry_2).flatten()
    #upper_triangle = c_matrix_[np.triu_indices(len(atom_coord_arry), k=1)]
 
    c_matrix_flat=np.unique(c_matrix_+2*smear*(np.random.rand( len(c_matrix_) )))
    
    return c_matrix_flat




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


def path_selection(positions_1,positions_2, bin_size,smear=0.15):
    # Compute the distance matrix
    #c_function_val=c_function(positions.tolist())
    c_function_val=c_function_2(positions_1,positions_2,smear=smear)
    hist, bin_edges=hist_info(c_function_val,bin_size)
    print(len(c_function_val))
    #hist=hist/len(c_function_val)
 
    return hist, bin_edges
        

def position_process(classical_coord_path,unit_size,length_end=1,bin_size=100,smear=0.15 ):
    bin_edges = np.linspace(0.1,length_end , int(bin_size+1))
    print('Delta length')
    print(bin_edges[1]-bin_edges[0])
    bin_edges_true = (bin_edges[1:] + bin_edges[:-1]) / 2

    structure_coord_full= np.loadtxt(classical_coord_path)
    print(structure_coord_full)
    structure_coord=structure_coord_full[:,1:4]
    structure_coord_name=structure_coord_full[:,0]

    num_structures=int(len(structure_coord)/unit_size)
    hist_arry=np.zeros( int(bin_size))
    hist_arry_OO=np.zeros( int(bin_size))
    hist_arry_OSi=np.zeros( int(bin_size))
    hist_arry_SiSi=np.zeros( int(bin_size))

    for i in range(0,num_structures):
        current_coord=structure_coord[int(i*unit_size):int((i+1)*unit_size) ]
        current_coord_name=structure_coord_name[int(i*unit_size):int((i+1)*unit_size) ]

        current_coord_O=current_coord[current_coord_name==2]
        current_coord_Si=current_coord[current_coord_name==1]

        hist, _= path_selection(current_coord,current_coord, bin_size=bin_edges,smear=smear)
        hist_arry=hist_arry+hist

        
        hist_OO, _= path_selection(current_coord_O,current_coord_O, bin_size=bin_edges,smear=smear)
        hist_arry_OO=hist_arry_OO+hist_OO
        
        hist_OSi, _= path_selection(current_coord_O,current_coord_Si, bin_size=bin_edges,smear=smear)
        hist_arry_OSi=hist_arry_OSi+hist_OSi
        
        hist_SiSi, _= path_selection(current_coord_Si,current_coord_Si, bin_size=bin_edges,smear=smear)
        hist_arry_SiSi=hist_arry_SiSi+hist_SiSi
        
        print(i)

    hist_arry[bin_edges_true<distance_-tol_down] = 0
    hist_arry_OO[bin_edges_true<distance_-tol_down] = 0
    hist_arry_OSi[bin_edges_true<distance_-tol_down] = 0
    hist_arry_SiSi[bin_edges_true<distance_-tol_down] = 0
    
    fig = plt.figure()
    ax1 = fig.add_subplot(311)

    hist_arry=(hist_arry )/(num_structures*bin_size)
    ax1.plot(bin_edges_true ,hist_arry, color='black',alpha=0.8)
    c_count_open_net_peaks, _ = find_peaks(hist_arry)
    for peak in c_count_open_net_peaks:
        ax1.text(bin_edges_true[peak], hist_arry[peak], str(np.round(bin_edges_true[peak],1)), ha='center', va='bottom')

    ax1.grid(True)
    #ax1.set_ylim(0,1)


    ax2 = fig.add_subplot(312)
    #ax2.set_ylim(0,1)
    hist_arry_OO=(hist_arry_OO )/(num_structures*bin_size)

    ax2.plot(bin_edges_true , hist_arry_OO , color='tab:blue',alpha=0.8,label='O-O')
    c_count_open_net_peaks, _ = find_peaks(hist_arry_OO)
    # for peak in c_count_open_net_peaks:
    #     ax2.text(bin_edges_true[peak], hist_arry_OO[peak], str(np.round(bin_edges_true[peak],1)), ha='center', va='bottom')
    hist_arry_OSi=(hist_arry_OSi )/(num_structures*bin_size)

    ax2.plot(bin_edges_true , hist_arry_OSi , color='tab:red',alpha=0.8,label='Si-O')
    c_count_open_net_peaks, _ = find_peaks(hist_arry_OSi)
    # for peak in c_count_open_net_peaks:
    #     ax2.text(bin_edges_true[peak], hist_arry_OSi[peak], str(np.round(bin_edges_true[peak],1)), ha='center', va='bottom')
    hist_arry_SiSi=(hist_arry_SiSi )/(num_structures*bin_size)

    ax2.plot(bin_edges_true , hist_arry_SiSi, color='tab:green',alpha=0.8,label='Si-Si')
    c_count_open_net_peaks, _ = find_peaks(hist_arry_SiSi)
    # for peak in c_count_open_net_peaks:
    #     ax2.text(bin_edges_true[peak], hist_arry_SiSi[peak], str(np.round(bin_edges_true[peak],1)), ha='center', va='bottom')

    ax2.grid(True)
    ax2.legend()


    ax3 = fig.add_subplot(313)

    ax3.plot(bin_edges_true , hist_arry_SiSi, color='tab:green',alpha=0.8,label='Si-Si')
    c_count_open_net_peaks, _ = find_peaks(hist_arry_SiSi)
    # for peak in c_count_open_net_peaks:
    #     ax2.text(bin_edges_true[peak], hist_arry_SiSi[peak], str(np.round(bin_edges_true[peak],1)), ha='center', va='bottom')

    ax3.grid(True)
    ax3.legend()

    # plot the closed loop
    fig_file_path=parent_dir+'/out/'+'tot_hist.png'
    plt.savefig(fig_file_path,dpi=500)
    plt.close()

    coord_file_path=parent_dir+'/out/'+'hist_tot_full_data.dat'
    np.savetxt(coord_file_path, hist_arry, fmt='%s', delimiter='\t')
    coord_file_path=parent_dir+'/out/'+'hist_O-O_full_data.dat'
    np.savetxt(coord_file_path, hist_arry_OO, fmt='%s', delimiter='\t')
    coord_file_path=parent_dir+'/out/'+'hist_O-Si_full_data.dat'
    np.savetxt(coord_file_path, hist_arry_OSi, fmt='%s', delimiter='\t')
    coord_file_path=parent_dir+'/out/'+'hist_Si-Si_full_data.dat'
    np.savetxt(coord_file_path, hist_arry_SiSi, fmt='%s', delimiter='\t')
    coord_file_path=parent_dir+'/out/'+'bin_edges_full_data.dat'
    np.savetxt(coord_file_path, bin_edges, fmt='%s', delimiter='\t')
    





file_path_input='input_file.dat'
print('SS')
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


smear=0

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


ring_num_distribution=[]
num_=[]
for kk in range(10,12,2):
    # try:
    current_data_dir=parent_dir+'n'+str(kk)+'/out/round_rough_real.dat'
    data=np.loadtxt(current_data_dir)
    ring_num_distribution.append(len(data))
    print('num')
    print(len(data))
    num_.append(kk)

coord_file_path=parent_dir+'/out/out/'+'ring_stats.dat'
np.savetxt(coord_file_path,np.array([ring_num_distribution,num_]).T, fmt='%s', delimiter='\t')
    # except:
    #     print('escape'+str(kk))

 