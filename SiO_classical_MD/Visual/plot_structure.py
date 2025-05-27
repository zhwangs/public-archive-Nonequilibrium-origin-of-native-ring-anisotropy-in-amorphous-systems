import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

x_max= 12.589462669
y_max=12.589462669
z_max=  13.222095665
def hist_info(data,bins=100):
    hist, bin_edges = np.histogram(data, bins=bins)
    return hist, bin_edges

# Load data from the text file
data = np.loadtxt('temperature.dat')

time_cut_index=0


structure_coord= np.loadtxt('time_evo_coord.txt')[1::]

structure_coord_name=np.loadtxt('time_evo_name.txt',dtype='str' )[1::]


structure_coord=structure_coord[time_cut_index::]
structure_coord_name=structure_coord_name[time_cut_index::]
# Only use in the case looking at partial structure
# unit_cell_correction=structure_coord[:,0]<0.1*x_max
# structure_coord[:,0][unit_cell_correction]=structure_coord[:,0][unit_cell_correction]+x_max

#file_path_coord_target='init_target_position.txt'


unit_step_structure_size=96


# Assuming data has two columns (x and y)
time = data[:, 0]  # First column
temp = data[:, 1]  # Second column

time=time[time_cut_index::]
temp=temp[time_cut_index::]
time_length=int(len(structure_coord)/unit_step_structure_size)
print(time_length)

unique_labels = np.unique(structure_coord_name)
cmap = plt.get_cmap('Set1')
 
#label_color_map = {unique_labels[i]: cmap(i) for i  in range(0,len(unique_labels))}
label_color_map ={'Si': 'blue', 'O': 'red', 'Yb': 'green'}
 
fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122)



size_marker=5
 
for i in range(0,len(unique_labels)):
 
    label_bool=structure_coord_name[0:unit_step_structure_size ]==unique_labels[i]
 
    cell_array=structure_coord[0:unit_step_structure_size ]
    plot=ax.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=10)
    ax.set_xlim(0,x_max)
    ax.set_ylim(0,y_max)
    ax.set_zlim(0,z_max)

ax2.plot(time,temp,c='blue')
ax2.scatter(time[0],temp[0],c='red',s=20)




imag_arry=[]#range(0,time_length)

def update(frame,structure_coord):
    frame=frame*10
    save_imag_loc='img/'
    #plt.cla()
    ax.clear()
    ax2.clear()
    current_structure_coord=structure_coord[ int(frame*unit_step_structure_size):int((frame+1)*unit_step_structure_size) ]
    for i in range(0,len(unique_labels)):
        
        label_bool=structure_coord_name[ int(frame*unit_step_structure_size):int((frame+1)*unit_step_structure_size) ]==unique_labels[i]
        cell_array=current_structure_coord
        plot=ax.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=size_marker)
        ax.set_xlim(0,x_max)
        ax.set_ylim(0,y_max)
        ax.set_zlim(0,z_max)
        ax.set_title(f'Time: {frame*time_length}')


        plt.legend()

        
        plot2=ax2.plot(time,temp,c='blue',alpha=0.3)

        scatter2=ax2.scatter(time[frame],temp[frame],c='red',s=20)
        ax2.set_xlabel('time')
        ax2.set_ylabel('Temperature')
        ax2.grid(True)



    if frame in imag_arry:
        plt.savefig(save_imag_loc+str(frame)+'.png',dpi=500)
  
    

    return [plot ,scatter2,plot2]


if time_length==len(time):
    for i in range(0,time_length):
        current_structure_coord=structure_coord[ int(i*unit_step_structure_size):int((i+1)*unit_step_structure_size) ]
        


# size_marker=5
 
# sample_loc=1000
# for i in range(0,len(unique_labels)):

#     current_structure_coord=structure_coord[ int(sample_loc*unit_step_structure_size):int((sample_loc+1)*unit_step_structure_size) ]

#     label_bool=structure_coord_name[ int(sample_loc*unit_step_structure_size):int((sample_loc+1)*unit_step_structure_size) ]==unique_labels[i]
 
#     cell_array=structure_coord[ int(sample_loc*unit_step_structure_size):int((sample_loc+1)*unit_step_structure_size) ]
#     plot=ax.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',color=label_color_map[unique_labels[i]],label=unique_labels[i],s=10)
#     ax.set_xlim(0,x_max)
#     ax.set_ylim(0,y_max)
#     ax.set_zlim(0,z_max)
 


ani = FuncAnimation(fig, update, frames=time_length, fargs=(structure_coord,), interval=1, save_count=500)

# # Set up formatting for the movie files
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)


# ani.save('lines.mp4', writer=writer)


#ani.save('T600_quenching.mp4', writer='ffmpeg', fps=10) 

plt.show()

#unit_step_structure_size


# for i in range(0,)


# # # Create a color map that maps each unique label to a unique color



 
# ax_3d.grid(True)
# plt.show()