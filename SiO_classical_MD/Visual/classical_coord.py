import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


classical_name_list=['Si','O'] # the coord number represents the type in the classical MD


def read_text_file_with_numpy(file_path ):
    try:
        # Use numpy.loadtxt to read the text file
        single_unit_cell = np.genfromtxt(file_path )
        print(single_unit_cell)
        # atom_order=single_unit_cell[:,0]
        # atom_number=len(atom_order)
        # atom_type=single_unit_cell[:,2]
        # atom_coord=single_unit_cell[:,4:7]
        # print(atom_coord)
        # combined_array=[]
        # full_unit_cell_name=[]
        # for i in range(0,atom_number):
        #     current_type=classical_name_list[int(atom_type[i]-1)]
        #     current_coord=atom_coord[i]
        #     full_unit_cell_name.append(current_type)
        #     current_combined_array = np.concatenate((np.array([current_type]), current_coord)).tolist()
        #     combined_array.append(current_combined_array)
        # print(np.array(combined_array))
        # output_file_path = 'classical_md_coord.dat'
        # np.savetxt(output_file_path, combined_array, fmt='%s', delimiter='\t')
        # full_unit_cell_name=np.array(full_unit_cell_name)
        # print(full_unit_cell_name)
        # # # # Create a color map that maps each unique label to a unique color
        # unique_labels = np.unique(full_unit_cell_name)
 
        # cmap = plt.get_cmap('tab10')
        # label_color_map = {label: cmap(i) for i, label in enumerate(unique_labels)}

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')

        # for label in unique_labels:
        #     label_bool=full_unit_cell_name==label
        #     cell_array=atom_coord
        #     ax.scatter(cell_array[:,0][label_bool],cell_array[:,1][label_bool],cell_array[:,2][label_bool] , marker='o',c=label_color_map[label],label=label)
        #     # ax.set_xlim(0,lattice_parameter[0][0]*nx)
        #     # ax.set_xlim(0,lattice_parameter[1][1]*ny)
        #     # ax.set_xlim(0,lattice_parameter[2][2]*nz)
        # plt.legend()
        # plt.show()

        # # Change this number for different lattice constant (single unit cell)
        # lat_in_A=au_2_A*4.6415377

        # lattice_parameter= np.loadtxt(file_path_lattice)
        # lattice_parameter=(lattice_parameter)*lat_in_A
        
        

        # cell_para_x = np.tile(lattice_parameter[0], (data_size, 1))
        # cell_para_y = np.tile(lattice_parameter[1], (data_size, 1))
        # cell_para_z = np.tile(lattice_parameter[2], (data_size, 1))


        # cell_val_x=np.tile(single_unit_cell[:,0], (3, 1)).T
        # cell_val_y=np.tile(single_unit_cell[:,1], (3, 1)).T
        # cell_val_z=np.tile(single_unit_cell[:,2], (3, 1)).T

        # cart_coord=(cell_para_x*cell_val_x+cell_para_y*cell_val_y+cell_para_z*cell_val_z)
        
        # print(lattice_parameter)
 
 
        # supercell_parameter=np.array([lattice_parameter[0]*nx,lattice_parameter[1]*ny,lattice_parameter[2]*nz])
 

        # cart_coord_supercell=cart_coord  

        # s=0
        # for i in range(0,nx):
        #     for j in range(0,ny):
        #         for z in range(0,nz):
        #             if int(i+j+z)==0 :
        #                 pass
        #             else:
        #                 s=s+1
        #                 added_cart_coord=((i)*cell_para_x+(j)*cell_para_y+(z)*cell_para_z)+cart_coord
        #                 cart_coord_supercell= np.vstack((cart_coord_supercell, added_cart_coord))
 


        # full_unit_cell_name=np.tile(single_unit_cell_element_name, n)
       
        # max_cart=np.sqrt(np.sum(cart_coord_supercell[:,0]**2+cart_coord_supercell[:,1]**2+cart_coord_supercell[:,2]**2))
        
        


        # combined_array = np.column_stack((full_unit_cell_name, cart_coord_supercell))

        # # Save the combined array to a text file
        # output_file_path = 'supercell.dat'
        # np.savetxt(output_file_path, combined_array, fmt='%s', delimiter='\t')

        # output_file_path = 'supercell_lattice_vectors.dat'
        # np.savetxt(output_file_path, supercell_parameter, fmt='%s', delimiter='\t')


        # # Get unique elements and their counts
        # unique_elements, counts = np.unique(full_unit_cell_name, return_counts=True)

        # # Create a dictionary to associate counts with strings
        # count_dict = dict(zip(unique_elements, counts))

        # # Save to a text file
        # output_file = "supercell_elements.txt"
        # np.savetxt(output_file, list(count_dict.items()), delimiter='\t', fmt='%s')



    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Provide the path to your text file
file_path = 'dump.lammpstrj'
# lattice_path='lattice_vectors.txt'
# nx=3
# ny=3
# nz=3
read_text_file_with_numpy(file_path )