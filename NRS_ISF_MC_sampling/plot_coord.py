import numpy as np
# Read the file
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

from scipy.optimize import minimize
import os, sys
import statsmodels.api as sm
from scipy.stats import skewnorm, norm, lognorm,loggamma
from scipy.optimize import curve_fit


# plot coordinates directly. It can be run after the file gen_ideal_structures_step1
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

def recenter(coord):
    center_coord=np.mean(coord,axis=0)
    coord[:,0]=coord[:,0]-center_coord[0]
    coord[:,1]=coord[:,1]-center_coord[1]
    coord[:,2]=coord[:,2]-center_coord[2]

    return coord

def skewnorm_pdf(x, scale, a,A,loc):
        
    return A*skewnorm.pdf(x,a=a,loc=loc,scale=scale)#A/(theta*np.sqrt(2*np.pi))*np.exp(-(x-k)**2/(2*theta**2))


def skewnorm_fit(xdata,ydata,p0):
    print(p0)
    try:
       print('gta')
       param_bounds = ([-100, -100,-100,-100],[100, 0,100,100]) 
       fitted_params,_ =curve_fit(skewnorm_pdf, xdata, ydata,p0=p0,bounds=param_bounds ,method='dogbox',maxfev = 100000)
       print(fitted_params)

    except:
        print('Error')
        #p0=[0.6501934458725108,0.02 ,0.00012406868608977807, 4.166291465288458]
        #fitted_params,_ =curve_fit(skewnorm_pdf, xdata, ydata,p0=p0,method='lm')

        fitted_params=p0

    return fitted_params

def max_boz_pdf(x, scale, a,A,loc):
    a=0
    val=A*x**2 * np.exp(-scale*(x-loc)**2)
    return val

def max_boz_fit(xdata,ydata):
    p0=[0.8,0 ,0.1, 0.1]
    
    try:
       fitted_params,_ =curve_fit(max_boz_pdf, xdata, ydata,p0=p0 ,method='lm')
    except:
        fitted_params=p0

    return fitted_params

def nearby_pair_distance(coord):
    # this is i+1 - i elements 
    coord=np.vstack([coord[-1],coord])
    
    diff = np.diff(coord,axis=0)
    
    return diff ,  np.linalg.norm(diff,axis=1) 

def generate_line(n, si_o_bond_length=1.62):
    # Initialize lists to store coordinates
    init_line_coord=np.zeros((n,3))
    for i in range(0,n):
        init_line_coord[i][0]=i*si_o_bond_length
    return recenter(init_line_coord)
    

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


def generate_zig_zag_ring(n, si_o_bond_length=1.62):
    half_n=int(n/2)
    # Calculate the radius of the ring
    radius = si_o_bond_length / (2 *np.sin(np.pi / half_n))
    
    # Initialize lists to store coordinates
    init_ring_coord=np.zeros((n,3))

    # Generate coordinates for top ring
    for i in range(0,half_n):
        index_val=int(i*2)
        angle = 2 * np.pi * index_val / n
        print(angle)

        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        z=si_o_bond_length*np.sqrt(3)/4
        init_ring_coord[index_val][0]=x
        init_ring_coord[index_val][1]=y
        init_ring_coord[index_val][2]=z

    for i in range(0,half_n):
        index_val=int(i*2+1)
        angle = 2 * np.pi * index_val / n
        

        shift_angle=2 * np.pi /half_n
        x = radius * np.cos(angle+shift_angle)
        y = radius * np.sin(angle+shift_angle)
        z= -si_o_bond_length*np.sqrt(3)/4
        
        init_ring_coord[index_val][0]=x
        init_ring_coord[index_val][1]=y
        init_ring_coord[index_val][2]=z


        

    return init_ring_coord, radius


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
 
def hist_info(data,bins=100):
    hist, bin_edges = np.histogram(data, bins=bins)
    return hist, bin_edges



def ring_stat_hist(current_ring,scaled_min,scaled_max ):
    # make it log and take the ratio with respet to the initial diff.

    resulted_ring= current_ring 
    print('log')
    # scaled_min=init_bin_min
    # scaled_max=init_bin_max


    # This is added for tol hist because they are real numbers 
    bin_edges = np.linspace(scaled_min,scaled_max, int(init_bin_size+1))
    print(bin_edges)
    print('Delta length')
    print(bin_edges[1]-bin_edges[0])
    bin_edges_true = (bin_edges[1:])


    data_y,_=hist_info(resulted_ring,bins=bin_edges)
    data_x=bin_edges_true
    print(data_y)
    return data_x, data_y

# input ------------------------------------------------------------------------------------------



def init_structure(replace=False):
    global n,file_path,bin_size,init_bin_min,init_bin_max,init_bin_size
    global si_o_bond_length,fit_hist_size,sampling_repeat,distance_max,distance_min,distance_bin_size
    global  sampling_number, max_angle,reference_number,bin_min,bin_max
    global distance_tol, min_angle,sample_repeat,l,critical_tol,parent_dir
    global max_step, max_l, min_l,total_sample_,bin_min_rough,bin_max_rough,init_diff_arry,init_diff_arry_r
    # Check if at least one command-line argument is provided


    if len(sys.argv) > 1:
        # The first command-line argument (sys.argv[1]) will be the input value
        
        parent_dir='./data/'+ sys.argv[1]+'/'
    else:
        print('you must provide argv name of the directory')
        
  

    file_path_input=parent_dir+'input_file.dat'
    input_file=np.genfromtxt(file_path_input )[:,1]
    
    n=int(input_file[0])
    si_o_bond_length=input_file[1]


    sampling_number=int(input_file[2])
    max_angle=np.pi*input_file[3]
    min_angle=np.pi*input_file[4]

    distance_tol=input_file[5]

    sampling_repeat=int(input_file[6])

    l=int(input_file[7])

    random_shift=int(input_file[8])

    reference_num=input_file[9]

    critical_tol=input_file[10]

    max_step=input_file[11]

    max_l=input_file[12]
    min_l=input_file[13]

    target_min=input_file[22]

    bin_size=int(input_file[24])
    bin_min=input_file[25]
    bin_max=input_file[26]
    init_bin_min=input_file[28]
    init_bin_max=input_file[29]
    init_bin_size=int(input_file[31])
    bin_min_rough=input_file[51]
    bin_max_rough=input_file[52]

    global init_bin_interval, bin_interval,bin_interval_rough
    init_bin_interval=init_bin_max-init_bin_min
    bin_interval=bin_max-bin_min
    bin_interval_rough=bin_max_rough-bin_min_rough
   

    try:
        distance_max= input_file[53] # max of the radial distance
        distance_min =input_file[54]# min of the radial distance
        distance_bin_size =input_file[55]  # bin size of the radial distance
    except:
        distance_max= 5 # max of the radial distance
        distance_min =0# min of the radial distance
        distance_bin_size =100  # bin size of the radial distance


   
    # overwrite 
    bin_interval=8
    bin_interval_rough=8#bin_max_rough-bin_min_rough
    init_bin_size=120
    bin_size=120
    distance_bin_size =100 




    


    # this is for the size of round and rough
    fit_hist_size=40


    file_path=parent_dir +'ring_statistics.dat'
    if replace:
        with open(file_path, 'w') as file:
            pass
    else:
        pass

    total_sample_=input_file[14]
    global full_coord_arry, angle_list_sample, round_arry,rough_arry,init_diff_arry

    coord_file_path=parent_dir+'round_rough.dat'
    rough_round_arry=np.loadtxt(coord_file_path)
    rough_round_arry = rough_round_arry[np.all(np.isfinite(rough_round_arry), axis=1)] 
    threshold = 0#10e-1
    # Keep rows where all values are greater than the threshold
    #rough_round_arry = rough_round_arry[np.all(rough_round_arry > threshold, axis=1)]

    angle_list_sample=rough_round_arry[:,0]
    round_arry=rough_round_arry[:,1]-1+threshold # add -1 so that it is zero when no anisotropy
    rough_arry=rough_round_arry[:,2]+threshold
    init_diff_arry=rough_round_arry[:,3]+threshold

    mask = np.isfinite(np.log(round_arry)) 
    round_arry= round_arry[mask]
    round_arry=round_arry[np.log(round_arry)>-30]

    mask = np.isfinite(np.log(rough_arry)) 
    rough_arry= rough_arry[mask]
    rough_arry=rough_arry[np.log(rough_arry)>-20]

    mask = np.isfinite(np.log(init_diff_arry)) 
    init_diff_arry= init_diff_arry[mask]
    init_diff_arry=init_diff_arry[np.log(init_diff_arry)>-20]



    coord_file_path=parent_dir+'coord.dat'
    full_coord_arry=np.loadtxt(coord_file_path)

    global full_coord_arry_r, angle_list_sample_r, round_arry_r,rough_arry_r

    coord_file_path_r=parent_dir+'round_rough_real.dat'
    rough_round_arry_r=np.loadtxt(coord_file_path_r)
    angle_list_sample_r=rough_round_arry_r[:,0]
    round_arry_r=rough_round_arry_r[:,1]-1 # add -1 so that it is zero when no anisotropy
    rough_arry_r=rough_round_arry_r[:,2]
    init_diff_arry_r=rough_round_arry_r[:,3]


    coord_file_path_r=parent_dir+'coord_real.dat'
    full_coord_arry_r=np.loadtxt(coord_file_path_r)




def radial_distance(full_coord):
    distance=np.sqrt(full_coord[:,0]**2+full_coord[:,1]**2+full_coord[:,2]**2)
    bin_edges = np.linspace(distance_min,distance_max, int(distance_bin_size+1))
    print(bin_edges)
    print('Delta length')
    print(bin_edges[1]-bin_edges[0])
    bin_edges_true = (bin_edges[1:])

    
    data_y,_=hist_info(distance,bins=bin_edges)

    return bin_edges_true, data_y

# end ------------------------------------------------------------------------------------------
 
#read file 
 
init_structure(replace=False)


string1 = 'Si'
string2 = 'O'
alternating_string = alternate_repeat(string1, string2, times=n)

 
_, radius=generate_ring(n, si_o_bond_length=si_o_bond_length )


# ploting 
unique_labels = np.unique(alternating_string)
cmap = plt.get_cmap('Set1')

label_color_map ={'Si': 'blue', 'O': 'red', 'Yb': 'green'}
 


fig = plt.figure(figsize=(12, 14))

ax1 = fig.add_subplot(4,3,1, projection='3d')
ax2 = fig.add_subplot(4,3,2, projection='3d')
ax3 = fig.add_subplot(4,3,3 , projection='3d')
ax4 = fig.add_subplot(4,3,4 )

ax5 = fig.add_subplot(4,3,5 )
ax6 = fig.add_subplot(4,3,6 )

ax7 = fig.add_subplot(4,3,7 )
ax8 = fig.add_subplot(4,3,8 )
ax9 = fig.add_subplot(4,3,9 )

ax10= fig.add_subplot(4,3,10 )
ax11= fig.add_subplot(4,3,11 )


def update(frame):
    global n,file_path,bin_size,init_bin_min,init_bin_max,init_bin_size
    global si_o_bond_length,fit_hist_size,sampling_repeat,distance_max,distance_min,distance_bin_size
    global  sampling_number, max_angle,reference_number,bin_min,bin_max
    global distance_tol, min_angle,sample_repeat,l,critical_tol,parent_dir
    global max_step, max_l, min_l,total_sample_,bin_min_rough,bin_max_rough,init_diff_arry,init_diff_arry_r
    global full_coord_arry_r, angle_list_sample_r, round_arry_r,rough_arry_r
    global full_coord_arry, angle_list_sample, round_arry,rough_arry,init_diff_arry

    global init_bin_interval, bin_interval,bin_interval_rough
    frame=frame
    plt.cla()
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax4.clear()

    # plot zig-zag ring 
    init_ring_coord, radius=generate_zig_zag_ring(n=n, si_o_bond_length=si_o_bond_length )
    ax1.set_title('number of ring element: '+str(n))
    cell_array=init_ring_coord
    for i in range(0,len(unique_labels)):
    
        label_bool=np.array(alternating_string)==unique_labels[i]
 
        
        ax1.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=15)

    
    ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(init_ring_coord)
 
    mean_coords=ring_center
    eigenvectors=eigenvectors
    eigenvalues=eigenvalues
    
    for i in range(3):
        # Plot eigenvector originating from the mean
        ax1.quiver(mean_coords[0], mean_coords[1], mean_coords[2],
                eigenvectors[0, i], eigenvectors[1, i], eigenvectors[2, i],
                color='r', length=eigenvalues[i], normalize=True ,linewidths=0.5)
    ax1.scatter(mean_coords[0],mean_coords[1],mean_coords[2] , marker='o',color='cyan',s=10)


    # plot ring 
    init_ring_coord, radius=generate_ring(n=n, si_o_bond_length=si_o_bond_length )
    ax1.set_title('number of ring element: '+str(n))
    cell_array=init_ring_coord
    for i in range(0,len(unique_labels)):
    
        label_bool=np.array(alternating_string)==unique_labels[i]
 
        
        ax1.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=15)

    
    ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(init_ring_coord)
 
    mean_coords=ring_center
    eigenvectors=eigenvectors
    eigenvalues=eigenvalues
    
    for i in range(3):
        # Plot eigenvector originating from the mean
        ax1.quiver(mean_coords[0], mean_coords[1], mean_coords[2],
                eigenvectors[0, i], eigenvectors[1, i], eigenvectors[2, i],
                color='r', length=eigenvalues[i], normalize=True ,linewidths=0.5)
    ax1.scatter(mean_coords[0],mean_coords[1],mean_coords[2] , marker='o',color='cyan',s=10)

    # Change view angle
    # ax1.view_init(elev=10, azim=80)  # Adjust azimuth to view from the side


    bin_edges_true, radial_dis=radial_distance(full_coord=full_coord_arry)
    ax10.scatter(bin_edges_true,radial_dis/len(radial_dis),s=1,c='tab:blue')

    # plot closed loop
    
    cell_array=full_coord_arry[ int((frame+1)*n):int((frame+2)*n) ]

    for i in range(0,len(unique_labels)):
    
        label_bool=np.array(alternating_string)==unique_labels[i]

        
        ax2.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=15)

    ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(cell_array)
 
    mean_coords=ring_center
    eigenvectors=eigenvectors
    eigenvalues=eigenvalues

    ax2.set_title('Simulated folding' )

    for i in range(3):
        # Plot eigenvector originating from the mean
        ax2.quiver(mean_coords[0], mean_coords[1], mean_coords[2],
                eigenvectors[0, i], eigenvectors[1, i], eigenvectors[2, i],
                color='r', length=eigenvalues[i], normalize=True ,linewidths=0.5)
    ax2.scatter(mean_coords[0],mean_coords[1],mean_coords[2] , marker='o',color='cyan',s=10)

    cell_array=init_ring_coord
    int_x=init_ring_coord[:,0].tolist()
    int_y=init_ring_coord[:,1].tolist()#.append(init_ring_coord[:,1][0])
    int_z=init_ring_coord[:,2].tolist()#.append(init_ring_coord[:,2][0])


    #ax1.plot(init_ring_coord[:,0] ,init_ring_coord[:,1] ,init_ring_coord[:,2]  ,color='black')
    ax1.plot(int_x,int_y ,int_z ,color='black')
    ax2.plot(int_x,int_y ,int_z ,color='black')



    bin_edges_true, radial_dis=radial_distance(full_coord=full_coord_arry_r)
    ax11.scatter(bin_edges_true,radial_dis/len(radial_dis),s=1,c='tab:red')
    # plot closed loop for r 
    cell_array=full_coord_arry_r[ int((frame+1)*n):int((frame+2)*n) ]
    # property 
    skewnorm_fit_arry=[n,max_angle,sampling_repeat]
    direct_fit_arry=[n,max_angle,sampling_repeat]
    percentile_fit_arry=[n,max_angle,sampling_repeat]
    bare_fit_arry=[n,max_angle,sampling_repeat]


    for i in range(0,len(unique_labels)):
    
        label_bool=np.array(alternating_string)==unique_labels[i]

        
        ax3.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=15)

    ring_center, eigenvalues, eigenvectors=calculate_eigenvalues_and_eigenvectors(cell_array)
 
    mean_coords=ring_center
    eigenvectors=eigenvectors
    eigenvalues=eigenvalues
    ax3.set_title('Ring Defect' )

    for i in range(3):
        # Plot eigenvector originating from the mean
        ax3.quiver(mean_coords[0], mean_coords[1], mean_coords[2],
                eigenvectors[0, i], eigenvectors[1, i], eigenvectors[2, i],
                color='r', length=eigenvalues[i], normalize=True ,linewidths=0.5)
    ax3.scatter(mean_coords[0],mean_coords[1],mean_coords[2] , marker='o',color='cyan',s=10)

    cell_array=init_ring_coord
    int_x=init_ring_coord[:,0].tolist()
    int_y=init_ring_coord[:,1].tolist()#.append(init_ring_coord[:,1][0])
    int_z=init_ring_coord[:,2].tolist()#.append(init_ring_coord[:,2][0])


    #ax1.plot(init_ring_coord[:,0] ,init_ring_coord[:,1] ,init_ring_coord[:,2]  ,color='black')
    ax1.plot(int_x,int_y ,int_z ,color='black')
    ax2.plot(int_x,int_y ,int_z ,color='black')
    ax3.plot(int_x,int_y ,int_z ,color='black')




    # real data is independent from the seeting in ring_stat_analysis

    # ax3.scatter(rough_arry,round_arry,c=angle_list_sample,s=20)

    # ax3.scatter(rough_arry[frame],round_arry[frame],c='black',s=50)
 
    range_sigma=bin_interval
    print('KSNFJKSD')

    print(np.log(round_arry))
    bin_center=np.mean(np.log(round_arry))
    bin_min=bin_center-range_sigma#*np.std(np.log(round_arry))
    bin_max=bin_center+range_sigma#*np.std(np.log(round_arry))
    print('SBA')
    print(bin_min)
    print(bin_max)

    bin_edges = np.linspace(bin_min,bin_max, int(bin_size+1))
    print(bin_edges)
    print('Delta length')
    print(bin_edges[1]-bin_edges[0])
    bin_edges_true = (bin_edges[1:])

    
    data_y,_=hist_info(np.log(round_arry),bins=bin_edges)
    data_x=bin_edges_true
    ax4.scatter(data_x,data_y/len(round_arry),s=1,c='tab:blue')

    mean_val=np.sum(data_x*(data_y/len(round_arry))/np.sum(data_y/len(round_arry)))
    p0=[2,-1.11,np.mean(data_y/len(round_arry) ), mean_val]

    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/len(round_arry),p0=p0 )
    
    #skewnorm_pdf(x, scale, a,A,loc):

    x_smooth=np.linspace(bin_min,bin_max,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    ax4.plot(x_smooth,y_smooth,c='tab:blue')

    # simulated round 1 0 
    a=fitted_params_skewnorm[-3]

    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)

    ax4.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:green", linestyle="--", linewidth=2, label="mean")
    ax4.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:green", linewidth=2, label="std")
    
    ax4.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    skewnorm_fit_arry.extend([a,loc,scale])

    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])


    # replaced with mean and variance 
    loc=np.mean(np.log(round_arry))
    scale=np.std(np.log(round_arry))
    direct_fit_arry.extend([a,loc,scale])


    data_y,_=hist_info(np.log(round_arry_r),bins=bin_edges)
    data_x=bin_edges_true
   # ax4.scatter(data_x,data_y/len(round_arry_r),s=1,c='tab:red')
    
    
    mean_val=np.sum(data_x*(data_y/len(round_arry_r) )/np.sum(data_y/len(round_arry_r) ))
    p0=[2,-1.11,np.mean(data_y/len(round_arry_r) ), mean_val]

    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/len(round_arry_r),p0=p0 )
    # real round 1 0 
    # simulated round 1 0 


    x_smooth=np.linspace(bin_min,bin_max,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    #ax4.plot(x_smooth,y_smooth,c='tab:red')
 

    a=fitted_params_skewnorm[-3]
    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)


    skewnorm_fit_arry.extend([a,loc,scale])
    # ax4.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:purple", linestyle="--", linewidth=2, label="mean")
    # ax4.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:purple", linewidth=2, label="std")
    # ax4.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])

    # replaced with mean and variance 
    loc=np.mean(np.log(round_arry_r))
    scale=np.std(np.log(round_arry_r))
    direct_fit_arry.extend([a,loc,scale])



    #ax4.text(1.1*fitted_params_skewnorm[-4], 0.8*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    #ax4.set_title('simulated: sigma: '+str(np.round(fitted_params_skewnorm[-4],2))+', a: '+str(np.round(fitted_params_skewnorm[-3],2))+', a/sigma: '+str(np.round(fitted_params_skewnorm[-3]/fitted_params_skewnorm[-4],3)) )

    # for rough

    range_sigma=bin_interval_rough
    bin_center=np.mean(np.log(rough_arry))
    bin_min_rough=bin_center-range_sigma#*np.std(np.log(rough_arry))
    bin_max_rough=bin_center+range_sigma#*np.std(np.log(rough_arry))
    print('SBA_rough')
    print(bin_min_rough)
    print(bin_max_rough)

    bin_edges = np.linspace(bin_min_rough,bin_max_rough, int(init_bin_size+1))
    print(bin_edges)
    print('Delta length')
    print(bin_edges[1]-bin_edges[0])
    bin_edges_true = (bin_edges[1:])


    data_y,_=hist_info(np.log(rough_arry),bins=bin_edges)
    data_x=bin_edges_true

    ax6.scatter(data_x,data_y/len(rough_arry),s=1,c='tab:blue')
    mean_val=np.sum(data_x*(data_y/len(rough_arry) )/np.sum(data_y /len(rough_arry)))
    p0=[1,-1.11,np.mean(data_y/len(rough_arry) ), mean_val]

    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/len(rough_arry),p0=p0 )
    
    x_smooth=np.linspace(bin_min_rough,bin_max_rough,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    ax6.plot(x_smooth,y_smooth,c='tab:blue')
 

    # simulated rough 1 1 
    a=fitted_params_skewnorm[-3]
    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)

    ax6.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:green", linestyle="--", linewidth=2, label="mean")
    ax6.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:green", linewidth=2, label="std")
    ax6.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    skewnorm_fit_arry.extend([a,loc,scale])

    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])
    # replaced with mean and variance 
    loc=np.mean(np.log(rough_arry))
    scale=np.std(np.log(rough_arry))
    direct_fit_arry.extend([a,loc,scale])

    
    

    data_y,_=hist_info(np.log(rough_arry_r),bins=bin_edges)
    data_x=bin_edges_true
    #ax6.scatter(data_x,data_y/len(rough_arry_r),s=1,c='tab:red')
    mean_val=np.sum(data_x*(data_y/len(rough_arry_r) )/np.sum(data_y/len(rough_arry_r) ))
    p0=[2,-1.11,np.mean(data_y/len(rough_arry_r) ), mean_val]

    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/len(rough_arry_r),p0=p0 )


    x_smooth=np.linspace(bin_min_rough,bin_max_rough,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    #ax6.plot(x_smooth,y_smooth,c='tab:red')
 

    a=fitted_params_skewnorm[-3]
    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)
 
    skewnorm_fit_arry.extend([a,loc,scale])
    # ax6.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:purple", linestyle="--", linewidth=2, label="mean")
    # ax6.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:purple", linewidth=2, label="std")
    # ax6.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])

    # replaced with mean and variance 
    loc=np.mean(np.log(rough_arry_r)) 
    scale=np.std(np.log(rough_arry_r))
    direct_fit_arry.extend([a,loc,scale])



   # ax6.set_title('sigma: '+str(np.round(fitted_params_skewnorm[-4],2))+', a: '+str(np.round(fitted_params_skewnorm[-3],2))+', a/sigma: '+str(np.round(fitted_params_skewnorm[-3]/fitted_params_skewnorm[-4],3)) )

    ax1.set_xlim(-radius,radius)
    ax1.set_ylim(-radius,radius)
    ax1.set_zlim(-radius,radius)

    ax2.set_xlim(-radius,radius)
    ax2.set_ylim(-radius,radius)
    ax2.set_zlim(-radius,radius)

    ax3.set_xlim(-radius,radius)
    ax3.set_ylim(-radius,radius)
    ax3.set_zlim(-radius,radius)


    ax4.set_xlabel('roundness')
    ax6.set_xlabel('roughness')

    # coord_file_path=parent_dir+'round_rough.dat'
    # rough_round_arry=np.loadtxt(coord_file_path)
    # init_diff_arry=np.log(rough_round_arry[:,3])

    init_diff_arry=np.log(init_diff_arry)
    range_sigma=init_bin_interval
    bin_center=np.mean((init_diff_arry))
    init_bin_min=bin_center-range_sigma#*np.std((init_diff_arry))
    init_bin_max=bin_center+range_sigma#*np.std((init_diff_arry))
    print('SBA')
    print(init_bin_min)
    print(init_bin_max)

    power_val=4/4#3/5
    # change into square 
    #init_diff_arry=init_diff_arry**2
    # add a power
    init_diff_arry=init_diff_arry**power_val

    data_x,data_y=ring_stat_hist(init_diff_arry,scaled_min=init_bin_min,scaled_max=init_bin_max)

    mean_val=np.sum(data_x*(data_y )/np.sum(data_y ))
    p0=[2,-1.11,np.mean(data_y/np.sum(data_y ) ), mean_val]

    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/np.sum(data_y ),p0=p0 )



    ax5.scatter(data_x,data_y/np.sum(data_y),c='tab:blue',s=3)
    
    x_smooth=np.linspace(init_bin_min,init_bin_max,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    ax5.plot(x_smooth,y_smooth,c='tab:blue')
 

    a=fitted_params_skewnorm[-3]
    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)
 
    ax5.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:green", linestyle="--", linewidth=2, label="mean")
    ax5.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:green", linewidth=2, label="std")
    ax5.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    skewnorm_fit_arry.extend([a,loc,scale])

    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])

    # replaced with mean and variance 
    loc=np.mean(init_diff_arry)
    scale=np.std(init_diff_arry)
    direct_fit_arry.extend([a,loc,scale])

    ax4.set_title('sigma: '+str(np.round(fitted_params_skewnorm[-4],2))+', a: '+str(np.round(fitted_params_skewnorm[-3],2))+', loc: '+str(np.round(fitted_params_skewnorm[-1] ,3)))
 
    ax7.set_title('Simulated folding vs skew-normal random variable')
    
    pp = sm.ProbPlot(init_diff_arry, fit=True)
    qq = pp.qqplot(marker='x', markerfacecolor='b', markeredgecolor='b', alpha=0.2, ax=ax7)
    qq = pp.qqplot(marker='x', markerfacecolor='b', markeredgecolor='b', alpha=0.2, ax=ax9)
  
    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    r_skew = skewnorm.rvs(a, size=len(init_diff_arry ))
    pp = sm.ProbPlot(r_skew, fit=True)
    qq = pp.qqplot(marker='.', markerfacecolor='k', markeredgecolor='k', alpha=0.05, ax=ax7)
    
    sm.qqline( ax=ax7, line='45', fmt='k--')


    coord_file_path=parent_dir+'round_rough_real.dat'
    rough_round_arry=np.loadtxt(coord_file_path)
    # since already log from the previous step, no need to take again. 
    init_diff_arry=np.log(rough_round_arry[:,3])
    # change into square 
    #init_diff_arry=init_diff_arry**2
    # add a power
    #init_diff_arry=init_diff_arry**power_val

    data_x,data_y=ring_stat_hist(init_diff_arry,scaled_min=init_bin_min,scaled_max=init_bin_max)

    mean_val=np.sum(data_x*(data_y )/np.sum(data_y ))
    p0=[2,-1.11,np.mean(data_y/np.sum(data_y ) ), mean_val]



    fitted_params_skewnorm=skewnorm_fit(xdata=data_x,ydata=data_y/np.sum(data_y ),p0=p0 )


    # ax5.scatter(data_x,data_y/np.sum(data_y),c='tab:red',s=3)
    
    x_smooth=np.linspace(init_bin_min,init_bin_max,1000)
    y_smooth=skewnorm_pdf(x_smooth,fitted_params_skewnorm[-4],fitted_params_skewnorm[-3],fitted_params_skewnorm[-2],fitted_params_skewnorm[-1])
    # ax5.plot(x_smooth,y_smooth,c='tab:red')
 
    # real distance 0 2
    a=fitted_params_skewnorm[-3]
    # replaced with max
    loc=np.sum(x_smooth*y_smooth/np.sum(y_smooth))
    scale=np.sqrt(np.sum((x_smooth**2)*y_smooth/np.sum(y_smooth))-loc**2)
 

    # ax5.vlines(x=loc, ymin=0, ymax=1.5*np.max(y_smooth), color="tab:purple", linestyle="--", linewidth=2, label="mean")
    # ax5.hlines(y=0.5*np.max(y_smooth), xmin=loc-scale, xmax=loc+scale, color="tab:purple", linewidth=2, label="std")
    # ax5.text(1.5*loc, 1.1*np.max(y_smooth), str(np.round(loc,2)), color="k", ha='center')  # Above the line at x=2

    skewnorm_fit_arry.extend([a,loc,scale])
    
    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    bare_fit_arry.extend([a,loc,scale])

    # replaced with mean and variance 
    loc=np.mean(init_diff_arry)
    scale=np.std(init_diff_arry)
    direct_fit_arry.extend([a,loc,scale])
    
    ax6.set_title('sigma: '+str(np.round(fitted_params_skewnorm[-4],2))+', a: '+str(np.round(fitted_params_skewnorm[-3],2))+', loc: '+str(np.round(fitted_params_skewnorm[-1] ,3)) )
    ax8.set_title('Ring Defect vs skew-normal random variable')


    pp = sm.ProbPlot(init_diff_arry, fit=True)
    qq = pp.qqplot(marker='x', markerfacecolor='r', markeredgecolor='r', alpha=0.2, ax=ax8)
    qq = pp.qqplot(marker='x', markerfacecolor='r', markeredgecolor='r', alpha=0.2, ax=ax9)
  
    a=fitted_params_skewnorm[-3]
    loc=fitted_params_skewnorm[-1]
    scale=fitted_params_skewnorm[-4]
    r_skew = skewnorm.rvs(a, size=len(init_diff_arry ))
    pp = sm.ProbPlot(r_skew, fit=True)
    qq = pp.qqplot(marker='.', markerfacecolor='k', markeredgecolor='k', alpha=0.05, ax=ax8)
    
    sm.qqline( ax=ax8, line='45', fmt='k--')
    
    ax9.set_title('Ring Defect vs simulated folding')


    # lower_cut_off=-7
    # upper_cut_off=7

    # check_range=np.array(skewnorm_fit_arry[3::])
    # condition = (check_range > upper_cut_off) | (check_range < lower_cut_off)
    # # in this case of directly calculating data mean, no need for this check. 
    # check_val=np.any(condition )

    # if not check_val:
    #     skew_out_arry=np.array(skewnorm_fit_arry).reshape(7,3)
    #     coord_file_path=parent_dir+'skew_out.dat'
    #     np.savetxt(coord_file_path, skew_out_arry.T, fmt='%s', delimiter='\t')
    # else:
    #     coord_file_path=parent_dir+'skew_out.dat'
    #     np.savetxt(coord_file_path, np.array([]))
    #     print('pass')

    skew_out_arry=np.array(skewnorm_fit_arry).reshape(7,3)
    coord_file_path=parent_dir+'skew_out.dat'
    np.savetxt(coord_file_path, skew_out_arry.T, fmt='%s', delimiter='\t')

    bare_out_arry=np.array(bare_fit_arry).reshape(7,3)
    coord_file_path=parent_dir+'bare_out.dat'
    np.savetxt(coord_file_path, bare_out_arry.T, fmt='%s', delimiter='\t')


    direct_out_arry=np.array(direct_fit_arry).reshape(7,3)
    coord_file_path=parent_dir+'direct_out.dat'
    np.savetxt(coord_file_path, direct_out_arry.T, fmt='%s', delimiter='\t')

    # percentile_out_arry=np.array(percentile_fit_arry).reshape(7,3)
    # coord_file_path=parent_dir+'percentile_out.dat'
    # np.savetxt(coord_file_path, percentile_out_arry.T, fmt='%s', delimiter='\t')

    sm.qqline( ax=ax9, line='45', fmt='k--')


    fig_file_path=parent_dir+'distance_distribution_.png'
    plt.savefig(fig_file_path,dpi=500)

    plt.close()
    
update(0)

# ani = FuncAnimation(fig, update, frames=sampling_number-1, interval=10)


 
# plt.show()
     