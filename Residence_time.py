'''
    Autor: Sergio Rodriguez Peña
'''



####################################################################
#                                                                  #
#   Code that analyzes the coordination times of a molecule with   # 
#   respect to the observed atoms/molecules by measuring the       #
#   continued residence time.                                      #
#                                                                  #
#   This code takes as input a process trajectory by the code      #
#   Process_Trajectory_PDB_to_TXT.py.                              #
#                                                                  #
####################################################################



''' IMPORTED PACKAGES '''

import numpy as np
import matplotlib.pyplot as plt



''' DEFINED FUNCTIONS '''

# Function that reads the first snapshot of the simulation and analyzes its composition.
# Then asks the user fot the reference atom/molecule (only 1) and the observable atoms/molecules (up to 5 different).
def Type_of_Interaction():
    molecule_id=[]
    atom_id=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])!=0:
                break
            if len(l)==7:
                if l[3] not in molecule_id:
                    molecule_id.append(l[3])
                    atom_id.append([])
                if l[1] not in atom_id[molecule_id.index(l[3])]:
                    atom_id[molecule_id.index(l[3])].append(l[1])
    print('\n')
    print('Molecule ID: Atom ID\n')
    for i in range(len(molecule_id)):
        print('%s: '%(molecule_id[i]), *atom_id[i])
    print('\n')
    reference=input('Choose Reference Molecule ID: ')
    obs=input('Choose Observable Atom/Molecule ID (up to 5 elements separated by comas): ')
    print('\n')
    observable=[]
    o=[]
    for i in range(len(obs)):
        if obs[i]==',':
            observable.append(''.join(o))
            o=[]
            continue
        o.append(obs[i])
    observable.append(''.join(o))
    return(reference,observable)


# Function that reads the trajectory and returns all the esential information of the system:
# time of the snapshots, size of the box, number of atoms, indexes of the molecules and atoms.
def System_Analysis():
    time=[]
    snap_offset=[]
    offset=0
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='t':
                time.append(float(l[2]))  
                snap_offset.append(offset)
            offset+=(len(line)+1)
    molecule_id=[]
    molecule_index=[]
    atom_id=[]
    atom_index=[]
    with open(input_file,'r') as r:
        for line in r:
            l=line.split()
            if l[0]=='Box':
                a_box=float(l[3])
            if len(l)==7:
                if l[3] not in molecule_id:
                    molecule_id.append(l[3])
                    molecule_index.append([])
                    atom_id.append([])
                    atom_index.append([])
                if l[1] not in atom_id[molecule_id.index(l[3])]:
                    atom_id[molecule_id.index(l[3])].append(l[1])
                    atom_index[molecule_id.index(l[3])].append([])
                if int(l[2]) not in molecule_index[molecule_id.index(l[3])]:
                    molecule_index[molecule_id.index(l[3])].append(int(l[2]))
                atom_index[molecule_id.index(l[3])][atom_id[molecule_id.index(l[3])].index(l[1])].append(int(l[0]))
            if l[0]=='t' and float(l[2])!=0:
                break 
    n_atoms=0
    for i in range(len(atom_index)):
        for j in range(len(atom_index[i])):
            n_atoms+=len(atom_index[i][j])
    return(time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms)


# Function that creates the hypercube. It gets the position of the reference atom/molecule without applying PBC
# and the position of the observable atoms/molecules applying PBC.
def Create_Hyper_Cube(t_index):
    hyper_cube=[]
    pos_ref=[]
    with open(input_file,'r') as r:
        r.seek(snap_offset[t_index])
        correct_snap=False
        for line in r:
            l=line.split()
            if l[0]=='t' and float(l[2])==time[t_index]:
                correct_snap=True
            if l[0]=='t' and float(l[2])!=time[t_index]:
                correct_snap=False
            if not correct_snap:
                break
            if correct_snap:
                if l[1] in observable or l[3] in observable:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    hyper_cube.append([x,y,z,l[0],l[1],l[2],l[3]])
                    pbc=Check_PBC([x,y,z])
                    if pbc!=[0,0,0]:
                        if pbc[0]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z,l[0],l[1],l[2],l[3]])
                        if pbc[1]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z,l[0],l[1],l[2],l[3]])
                        if pbc[2]!=0:
                            hyper_cube.append([x,y,z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z,l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y,z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x,y+a_box*pbc[1],z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                        if pbc[0]!=0 and pbc[1]!=0 and pbc[2]!=0:
                            hyper_cube.append([x+a_box*pbc[0],y+a_box*pbc[1],z+a_box*pbc[2],l[0],l[1],l[2],l[3]])
                if l[1]==reference or l[3]==reference:
                    x=float(l[4])
                    y=float(l[5])
                    z=float(l[6])
                    pos_ref.append([x,y,z,l[0],l[1],l[2],l[3]])
    return(hyper_cube,pos_ref)


# Function that checks if a position is close enought to the edges of the box and needs to be applied PBC.
def Check_PBC(pos):
    pbc=[]
    for i in pos:
        if i>(a_box-cutoff+0.01):
            pbc.append(-1)
        elif i<(cutoff+0.01):
            pbc.append(1)
        else:
            pbc.append(0)
    return(pbc)


# Function that creates the progress bar shown in the terminal.
def Progress_Bar(progress,total):
    percent=100*(progress/float(total))
    bar='|'+'#'*int(percent/2)+'-'*(50-(int(percent/2)))+'|'
    print(f'%s %.1f' %(bar,percent) +' %' , end='\r')


# Funcion que analiza todos los atomos del hyper_cube para comprobar si se encuentran dentro
# del cutoff del atomo referencia y guarda los datos del observable.
def Check_Observable(i,j,pos_ref,coordination,coordination_info,hyper_cube):
    for k in np.arange(len(hyper_cube)):
        pos_obs=hyper_cube[k][0:3]
        if abs(pos_obs[0]-pos_ref[0])>cutoff:
            continue
        if abs(pos_obs[1]-pos_ref[1])>cutoff:
            continue
        if abs(pos_obs[2]-pos_ref[2])>cutoff:
            continue
        else:
            dis=np.sqrt((pos_ref[0]-pos_obs[0])**2+(pos_ref[1]-pos_obs[1])**2+(pos_ref[2]-pos_obs[2])**2)
            if dis<=cutoff:
                if hyper_cube[k][4] in observable:
                    coordination[j][i].append(int(hyper_cube[k][3]))
                    coordination_info[j][i].append(hyper_cube[k][4])
                if hyper_cube[k][6] in observable:
                    coordination[j][i].append(int(hyper_cube[k][5]))
                    coordination_info[j][i].append(hyper_cube[k][6])
    return(coordination,coordination_info)



''' CODE '''

# Name of the input file.
print('\n')
entr=input('Trayectory File Name: ')
input_file=entr+'.txt'

# Ask for the reference and observables.
reference,observable=Type_of_Interaction()

# Get the information of the system.
time,snap_offset,a_box,molecule_id,atom_id,molecule_index,atom_index,n_atoms=System_Analysis()

# Ask for the cutoff radius.
print('\n')
cutoff=float(input('Cutoff Distance (A): '))

# Analysis of the trajectory.
# For each time step, create the hypercube and check the distances between reference and observables.
print('\n')
print('ANALYSING COORDINATION ENVIRONMENT:')
coordination=[]
coordination_info=[]
for i in np.arange(len(molecule_index[molecule_id.index(reference)])):
    coordination.append([])
    coordination_info.append([])
for i in np.arange(len(time)):
    for j in np.arange(len(molecule_index[molecule_id.index(reference)])):
        coordination[j].append([])
        coordination_info[j].append([])
    hyper_cube,position_ref=Create_Hyper_Cube(i)
    for j in np.arange(len(molecule_index[molecule_id.index(reference)])):
        pos_ref=position_ref[j][:3]
        coordination,coordination_info=Check_Observable(i,j,pos_ref,coordination,coordination_info,hyper_cube)
    Progress_Bar(i,len(time))  
bar='|'+'#'*50+'|'
print(f'%s 100.0' %(bar) +' %' , end='\r')

# Ordenamos los contactos en funcion de a que observable pertenecen.
contact_ordered=[]
for i in np.arange(len(coordination)):
    contact_ordered.append([])
    for j in np.arange(len(coordination[i])):
        contact_ordered[i].append([])
        for k in np.arange(len(observable)):
            contact_ordered[i][j].append([])
        for k in np.arange(len(coordination[i][j])):
            if coordination[i][j][k] not in contact_ordered[i][j][observable.index(coordination_info[i][j][k])]:
                contact_ordered[i][j][observable.index(coordination_info[i][j][k])].append(coordination[i][j][k])

# Medimos los tiempos de residencia de cada contacto.
residence_time=[]
for i in np.arange(len(observable)):
    residence_time.append([])
for i in np.arange(len(contact_ordered)):
    for j in np.arange(len(contact_ordered[i])-1):
        for k in np.arange(len(contact_ordered[i][j])):
            if len(contact_ordered[i][j][k])==0:
                continue
            for m in np.arange(len(contact_ordered[i][j][k])):
                ref=contact_ordered[i][j][k][m]
                if j!=0:
                    if ref in contact_ordered[i][j-1][k]:
                        continue
                r_t=0
                for n in np.arange(j+1,len(contact_ordered[i])-1):
                    if ref in contact_ordered[i][n][k]:
                        r_t+=1
                        if n==len(contact_ordered[i])-2:
                            residence_time[k].append(r_t)
                        else:
                            continue
                    if ref in contact_ordered[i][n+1][k]:
                        r_t+=1
                        if n==len(contact_ordered[i])-2:
                            residence_time[k].append(r_t)
                        else:
                            continue                    
                    if ref not in contact_ordered[i][n][k] and ref not in contact_ordered[i][n+1][k]:
                        if r_t>1:
                            residence_time[k].append(r_t)
                            break

# Calculamos la distribucion de los tiempos de residencia.
distribution=[]
for i in np.arange(len(observable)):
    distribution.append(np.zeros((len(time)-2)))
for i in np.arange(len(residence_time)):
    for j in np.arange(len(residence_time[i])):
        distribution[i][residence_time[i][j]-2]+=1 

# Calculamos la media.
medias=[]
error=[]
for i in np.arange(len(distribution)):
    a=0
    for j in np.arange(len(distribution[i])):
        a=a+distribution[i][j]*(j+1)*((time[1]-time[0])/1000)
    medias.append(a/np.sum(distribution[i]))
    b=0
    for j in np.arange(len(distribution[i])):
        b=b+(distribution[i][j]*(j+1)*((time[1]-time[0])/1000)-medias[i])**2    
    error.append(np.sqrt(b)/len(distribution[i]))

# Mostramos los resultados en pantalla.
print('\n')
for i in np.arange(len(observable)):
    print('Avg. Residence Time %s-%s = %.3f +- %.3f ns' %(reference,observable[i],medias[i],error[i]))

# Ploteamos la distribucion 
x_plot=[]
for i in np.arange(2,len(time)):
    x_plot.append(time[i]/1000)
for i in np.arange(len(distribution)):
    plt.figure()
    plt.scatter(x_plot,distribution[i])

plt.show(block=False)
