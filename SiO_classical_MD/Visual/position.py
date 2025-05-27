import re
import numpy as np
# Read the file
file_path = 'dump.lammpstrj'

file_path_coord='cart_coordinates.txt'

classical_name_list=['Si','O'] # the coord number represents the type in the classical MD


init_name_arry= np.loadtxt(file_path_coord,dtype='str' )[:,0]
init_coord_arry=np.genfromtxt(file_path_coord )[:,1:4]
data_size=len(init_name_arry)

data_size=96
file_path_temp='temperature.dat'
temp_array=np.genfromtxt(file_path_temp )[:,1].tolist()



def read_atom_num(line):
    try: 
        parts = line.split()
        if int(parts[0])==data_size:
            return True
            
    except:
        return False
    

def read_time(line,time_array):
    if 'Timestep' in line and 'Atoms' in line:
        # Split the line by spaces and extract the numerical part
        parts = line.split()
        for part in parts:
            try:
                # Try to convert the part to a float
                time_value = float(part)

                time_array.append(time_value)
                
                break  # Stop searching after finding the numerical value
            except ValueError:
                pass  # Ignore non-numeric parts
    return time_array
# Open the file and read line by line


def atomic_position(line):
    S=line.strip() 
    return S


time_array=[]

total_atom_coord_arry=np.zeros(3) # changed to target sorted 
total_atom_name_arry = np.array(['start'])  # changed to target sorted 

 
atom_name_arry=[]
atom_coord_arry=[]

global_count=0
with open(file_path, 'r') as file:
    counter=0
    for line in file:
        bool_head=read_atom_num(line)
        #print(bool_head)
        if bool_head:
            counter=1
            #print('SSS')

        elif counter==1:
            time_array=read_time(line,time_array)
            counter=2
            atom_name_arry=[]
            atom_coord_arry=[]
            
            
        elif counter<data_size+1:
            parts = line.split()
            print(counter)
            #print(parts)
            print(parts[0])
            atom_name = classical_name_list[int(parts[0])-1]  # First part is the atom
            # print(int(parts[0]))
            coordinates = np.array([float(coord) for coord in parts[1:]]) 
            
            atom_coord_arry.append(coordinates)
            atom_name_arry.append(atom_name)
            
            counter+=1
        

        elif counter==data_size+1:
            counter=0
            #print(atom_coord_arry)
            print(np.array(atom_name_arry))
            
            total_atom_name_arry=np.hstack([total_atom_name_arry,np.array(atom_name_arry)])
            print(total_atom_name_arry)
            total_atom_coord_arry=np.vstack([total_atom_coord_arry,np.array(atom_coord_arry)])
            print(total_atom_coord_arry)
        global_count+=1


print(np.array(time_array).shape)
print(np.array(temp_array).shape)
print((total_atom_name_arry).shape)
print(total_atom_coord_arry)



coord_file_path='time_evo_coord.txt'
np.savetxt(coord_file_path, total_atom_coord_arry, fmt='%s', delimiter='\t')
name_file_path='time_evo_name.txt'
np.savetxt(name_file_path, total_atom_name_arry, fmt='%s')


