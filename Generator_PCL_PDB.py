'''
    Autor: Sergio Rodriguez Peña
'''



###############################################################################
#                                                                             #
#   Code that builds a PCL chain in format .pdb with the number of monomers   #
#   required by the user. The code randomly creates the atoms connecting      #
#   them with short bond lenghts.                                             #
#                                                                             #
#   PCL Chain:    CH3 - (CH2 - CH2 - CH2 - CH2 - CH2 - C(=O) - O)n - CH3      #
#                                                                             #
#   The resulting chain must be optimized with other software.                #
#                                                                             #
###############################################################################



''' IMPORTED PACKAGES '''

import numpy as np



''' DEFINED FUNCTIONS '''

# Function that updates the position of the next created atom by randomly
# displaced it a distance of 1.2 A from the previous one.
def pos_act(pos):
	x=float(np.random.rand(1))
	r=np.random.rand(1)
	if r<0.5:
		x=-x
	y_corr=False
	while not y_corr:
		y=float(np.random.rand(1))
		s=np.random.rand(1)
		if s<0.5:
			y=-y
		if np.sqrt(x**2+y**2)<1:
			y_corr=True
	z=np.sqrt(1-x**2-y**2)
	t=np.random.rand(1)
	if t<0.5:
		z=-z
	dist_enlace=1.2
	f=np.sqrt(dist_enlace)
	vector=np.array([f*x,f*y,f*z])
	pos_act=[pos[0]+vector[0],pos[1]+vector[1],pos[2]+vector[2]]
	return (pos_act)


# Function that creates the line of the new atom with its information and 
# required spaces to align the columns.
def new_atom(atom_number,atom_id,pos):
	a_n=(7-len(atom_number))*' '+atom_number
	a_i=2*' '+atom_id
	s_peo=(4-len(atom_id))*' '
	if pos[0]>-1 and pos[0]<0:
		s_x=(7-len(str(int(pos[0]))))*' '+'{:.3f}'.format(pos[0])
	else:
		s_x=(8-len(str(int(pos[0]))))*' '+'{:.3f}'.format(pos[0])
	if pos[1]>-1 and pos[1]<0:
		s_y=(3-len(str(int(pos[1]))))*' '+'{:.3f}'.format(pos[1])
	else:
		s_y=(4-len(str(int(pos[1]))))*' '+'{:.3f}'.format(pos[1])
	if pos[2]>-1 and pos[2]<0:
		s_z=(3-len(str(int(pos[2]))))*' '+'{:.3f}'.format(pos[2])
	else:
		s_z=(4-len(str(int(pos[2]))))*' '+'{:.3f}'.format(pos[2])
	new_atom='ATOM'+a_n+a_i+s_peo+'PCL     1'+s_x+s_y+s_z+'  1.00  0.00           '+str(atom_id)
	return(new_atom)


# Function that creates the new CONECT line of the .pdb file with the required
# spcaes to align the columns.
def new_connect(indexes):
	ind=[]
	for i in np.arange(len(indexes)):
		if indexes[i]<10:
			ind.append('    '+str(indexes[i]))
		if indexes[i]>=10 and indexes[i]<100:
			ind.append('   '+str(indexes[i]))			
		if indexes[i]>=100:
			ind.append('  '+str(indexes[i]))		
	a=''
	for i in np.arange(len(ind)):
		a=a+ind[i]
	new_connnect='CONECT'+a
	return(new_connnect)



''' CODE '''

# Structure of the polymer.
print('\n')
print('PCL Structure: CH3 - (CH2 - CH2 - CH2 - CH2 - CH2 - C(=O) - O)n - CH3 \n')

# Create the strcuture .pdb file. and write the atomic information.
num_monomers=int(input('Number of monomers (n): '))
f=open('pcl.pdb','w')
f.write('COMPND    UNNAMED\n')
pos=np.array([0,0,0])
new=new_atom('1','C',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom('2','H',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom('3','H',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom('4','H',pos)
f.write(new+'\n')
pos=pos_act(pos)
atoms=['C','H','H','C','H','H','C','H','H','C','H','H','C','H','H','C','O','O']
for i in np.arange(num_monomers):
	for j in np.arange(len(atoms)):
		pos=pos_act(pos)
		atom=str(j+18*i+5)
		new=new_atom(atom,atoms[j],pos)
		f.write(new+'\n')
pos=pos_act(pos)
new=new_atom(str(5+num_monomers*18),'C',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom(str(6+num_monomers*18),'H',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom(str(7+num_monomers*18),'H',pos)
f.write(new+'\n')
pos=pos_act(pos)
new=new_atom(str(8+num_monomers*18),'H',pos)
f.write(new+'\n')

# Write the CONECT section.
f.write('CONECT    1    2    3    4    5 \n')
f.write('CONECT    2    1 \n')
f.write('CONECT    3    1 \n')
f.write('CONECT    4    1 \n')
for i in np.arange(num_monomers):
	atom=[]
	for j in np.arange(18):
		atom.append(j+18*i+5)
	if i==0:
		new=new_connect([atom[0],atom[0]-4,atom[0]+1,atom[0]+2,atom[0]+3])	
	else:
		new=new_connect([atom[0],atom[0]-1,atom[0]+1,atom[0]+2,atom[0]+3])	
	f.write(new+'\n')
	new=new_connect([atom[1],atom[1]-1])	
	f.write(new+'\n')
	new=new_connect([atom[2],atom[2]-2])	
	f.write(new+'\n')
	new=new_connect([atom[3],atom[3]-3,atom[3]+1,atom[3]+2,atom[3]+3])	
	f.write(new+'\n')
	new=new_connect([atom[4],atom[4]-1])	
	f.write(new+'\n')
	new=new_connect([atom[5],atom[5]-2])	
	f.write(new+'\n')
	new=new_connect([atom[6],atom[6]-3,atom[6]+1,atom[6]+2,atom[6]+3])	
	f.write(new+'\n')
	new=new_connect([atom[7],atom[7]-1])	
	f.write(new+'\n')
	new=new_connect([atom[8],atom[8]-2])	
	f.write(new+'\n')
	new=new_connect([atom[9],atom[9]-3,atom[9]+1,atom[9]+2,atom[9]+3])	
	f.write(new+'\n')
	new=new_connect([atom[10],atom[10]-1])	
	f.write(new+'\n')
	new=new_connect([atom[11],atom[11]-2])	
	f.write(new+'\n')
	new=new_connect([atom[12],atom[12]-3,atom[12]+1,atom[12]+2,atom[12]+3])	
	f.write(new+'\n')
	new=new_connect([atom[13],atom[13]-1])	
	f.write(new+'\n')
	new=new_connect([atom[14],atom[14]-2])	
	f.write(new+'\n')
	new=new_connect([atom[15],atom[15]-3,atom[15]+1,atom[15]+2])	
	f.write(new+'\n')
	new=new_connect([atom[16],atom[16]-1])	
	f.write(new+'\n')
	new=new_connect([atom[17],atom[17]-2,atom[17]+1])	
	f.write(new+'\n')
new=new_connect([num_monomers*18+4+1,num_monomers*18+4-1,num_monomers*18+4+2,num_monomers*18+4+3,num_monomers*18+4+4])	
f.write(new+'\n')
new=new_connect([num_monomers*18+4+2,num_monomers*18+4+1])	
f.write(new+'\n')
new=new_connect([num_monomers*18+4+3,num_monomers*18+4+1])	
f.write(new+'\n')
new=new_connect([num_monomers*18+4+4,num_monomers*18+4+1])	
f.write(new+'\n')

# Write the final information of the .pdb file.
f.write('MASTER        0    0    0    0    0    0    0    0  %i    0  %i    0 \n' %(8+num_monomers*18,8+num_monomers*18))
f.write('END')
f.close()

# Clarification.
print('\n')
print('This polymer has not the correct distance betwwen atoms.')
print('Quick geomtery optimization is needed. \n')