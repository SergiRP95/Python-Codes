import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import pandas as pd
from brokenaxes import brokenaxes

def check_isolated(cl):
	iso_cation=0
	for i in np.arange(len(cl)):
		if len(cl[i])==1:
			iso_cation+=1
	a=0
	for i in np.arange(len(cl)):
		a=a+len(cl[i])
	iso_anion=len(index_cation)+len(index_anion)-a
	return(iso_anion,iso_cation)

def check_limits(numero):
	m=[]
	for i in np.arange(len(numero)):
		if numero[i]==0:
			m.append(0)
		if numero[i]!=0:
			m=[]
		if m==[0,0,0,0,0]:
			minimo=i-2
			break
	m=[]
	empezar=False
	start=True
	for i in np.arange(len(numero)):
		if start:
			if numero[-i]!=0:
				empezar=True
				start=False
				continue
		if empezar:
			if numero[-i]==0:
				m.append(0)
			if numero[-i]!=0:
				m=[]
			if m==[0,0,0,0,0]:
				maximo=len(numero)-i-2
				break			
	return(minimo,maximo)


entrada=input('Cluster File Name: ')+'.txt'
print('\n')
check=input('Make checking (y/n): ')
dat=open(entrada)
datos=[]
for lin in dat:
	l=lin.split()
	datos.append(l)

clust=[]
t=[]
for i in np.arange(len(datos)):
	if len(datos[i])==0 or datos[i][0]=='Anion:' or datos[i][0]=='Cation:':
		continue
	elif datos[i][0]=='t':
		clust.append([])
		t.append(datos[i][2])
	else:
		clust[len(t)-1].append(datos[i])
index_anion=np.arange(int(datos[0][1]),int(datos[0][3])+1)
index_cation=np.arange(int(datos[1][1]),int(datos[1][3])+1)

if check=='y':
	a=0
	b=0
	for i in np.arange(len(clust)):
		a=a+check_isolated(clust[i])[0]
		b=b+check_isolated(clust[i])[0]
		b=b+len(clust[i])
		for j in np.arange(len(clust[i])):
			a=a+len(clust[i][j])	
	print('Total # Ions = %i' %(a))
	print('Total # Clusters = %i' %(b))
	print('\n')
'''                  Columnas = # Aniones 
Filas = # Litios                              '''

clust_matrix=[]
max_n=0
for i in np.arange(len(clust)):
	big=False
	for j in np.arange(len(clust[i])):
		if len(clust[i][j])>100:
			clust_matrix.append(np.zeros((len(index_cation)+1,len(index_cation)+1)))
			big=True
	if not big:
		clust_matrix.append(np.zeros((100,100)))
	iso_anion,iso_cation=check_isolated(clust[i])
	clust_matrix[i][0][1]=iso_anion
	clust_matrix[i][1][0]=iso_cation
	for j in np.arange(len(clust[i])):
		if len(clust[i][j])==1:
			continue
		if len(clust[i][j])>1:
			n_l=0
			n_a=0
			for k in np.arange(len(clust[i][j])):
				if int(clust[i][j][k]) in index_cation:
					n_l+=1
				if int(clust[i][j][k]) in index_anion:
					n_a+=1
		while n_l>len(clust_matrix[i])-1 or n_a>len(clust_matrix[i][0])-1:
			clust_matrix[i]=np.insert(clust_matrix[i],-1,np.zeros((len(clust_matrix[i][0]))),0)
			clust_matrix[i]=np.insert(clust_matrix[i],-1,np.zeros((len(clust_matrix[i]))),1)
		clust_matrix[i][n_l][n_a]+=1
		if n_l>max_n:
			max_n=n_l
		if n_a>max_n:
			max_n=n_a
if check=='y':
	c=0
	d=0
	for i in np.arange(len(clust_matrix)):
		for j in np.arange(len(clust_matrix[i])):
			for k in np.arange(len(clust_matrix[i][j])):
				c=c+clust_matrix[i][j][k]*(j+k)
		d=d+np.sum(clust_matrix[i])
	print('Total # Ions = %i' %(c))
	print('Total # Clusters = %i' %(d))
	print('\n')

total_clust_matrix=np.zeros((max_n+1,max_n+1))
if not big:	
	for i in np.arange(max_n+1):
		for j in np.arange(max_n+1):
			a=[]
			for k in np.arange(len(clust_matrix)):
				if i<len(clust_matrix[k]):
					if j<len(clust_matrix[k][0]):
						a.append(clust_matrix[k][i][j])
			total_clust_matrix[i][j]=np.sum(a)
if big:
	for i in np.arange(len(clust_matrix)):
		if np.sum(clust_matrix[i])==1:
			for j in np.arange(len(index_cation)-5,len(index_cation)+1):
				for k in np.arange(len(index_cation)-5,len(index_cation)+1):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
		if np.sum(clust_matrix[i])>1:
			for j in np.arange(int(np.sum(clust_matrix[i]))+100):
				for k in np.arange(int(np.sum(clust_matrix[i]))+100):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
			for j in np.arange(len(index_cation)-int(np.sum(clust_matrix[i]))-100,len(index_cation)+1):
				for k in np.arange(len(index_cation)-int(np.sum(clust_matrix[i]))-100,len(index_cation)+1):
					total_clust_matrix[j][k]=total_clust_matrix[j][k]+clust_matrix[i][j][k]
if check=='y':
	e=0
	f=np.sum(total_clust_matrix)
	for i in np.arange(len(total_clust_matrix)):
		for j in np.arange(len(total_clust_matrix[i])):
			e=e+total_clust_matrix[i][j]*(i+j)
	print('Total # Ions = %i' %(e))
	print('Total # Clusters = %i' %(f))
	print('\n')

'''
if not big:
	index=[]
	for i in np.arange(len(total_clust_matrix)):
		index.append('%s' %(i))	
	df=pd.DataFrame(total_clust_matrix,index=index)
	df.columns=pd.MultiIndex.from_product([['Number of Anions'],df.columns])
	df.index=pd.MultiIndex.from_product([['Number of Lithiums'],df.index])
	negative=0
	positive=0
	neutral=0
	for i in np.arange(len(total_clust_matrix)):
		for j in np.arange(len(total_clust_matrix)):
			if i>j:
				positive=positive+total_clust_matrix[i][j]
			if i<j:
				negative=negative+total_clust_matrix[i][j]
			if i==j:
				neutral=neutral+total_clust_matrix[i][j]
	print ('\n')
	print('Total Frames: %i' %len(clust))
	print ('\n')
	print('Total Clusters: %i' %np.sum(total_clust_matrix))
	print ('\n')
	print(df)
	print('\n')
	print('Total Negative Clusters: %i (Real: %i)' %(negative,negative-total_clust_matrix[0][1]))
	print('\n')
	print('Total Positive Clusters: %i (Real: %i)' %(positive,positive-total_clust_matrix[1][0]))
	print('\n')
	print('Total Neutral Clusters: %i' %(neutral))
'''
do_print=input('Print cluster map? (y/n): ')
if do_print=='y':
	X,Y=np.meshgrid(np.arange(len(total_clust_matrix)),np.arange(len(total_clust_matrix)))
	Z=np.zeros((len(total_clust_matrix),len(total_clust_matrix)))
	for i in np.arange(len(total_clust_matrix)):
		for j in np.arange(len(total_clust_matrix)):
			Z[i][j]=100*total_clust_matrix[i][j]/np.sum(total_clust_matrix)
	minimo=10
	for i in np.arange(len(Z)):
		for j in np.arange(len(Z)):
			if Z[i][j]==0:
				continue
			else:
				if Z[i][j]<minimo:
					minimo=Z[i][j]
	fig,ax=plt.subplots()
	cl=ax.pcolormesh(X,Y,Z,norm=colors.LogNorm(vmin=1e-3,vmax=1),shading='auto')
	cbar = plt.colorbar(cl,extend='max')
	for t in cbar.ax.get_yticklabels():
	     t.set_fontsize(22)
	plt.plot(np.arange(len(total_clust_matrix)),np.arange(len(total_clust_matrix)),linewidth=6,color='k')
	plt.xlabel(r'# $TFSI^{-}$',fontsize=30)
	plt.ylabel(r'# $Li^{+}$',fontsize=30)
	plt.xticks(ticks=np.arange(len(total_clust_matrix)+1),fontsize=16)
	plt.yticks(ticks=np.arange(len(total_clust_matrix)+1),fontsize=16)
	plt.grid(visible=True,color='gray')
	plt.xlim(-0.5,len(total_clust_matrix)+0.5)
	plt.ylim(-0.5,len(total_clust_matrix)+0.5)
	plt.show(block=False)

'''
numero=np.zeros((len(index_cation)+len(index_anion)+1))
for i in np.arange(len(total_clust_matrix)):
	for j in np.arange(len(total_clust_matrix[i])):
		numero[i+j]=numero[i+j]+total_clust_matrix[i][j]
prob=100*numero/np.sum(numero)
if not big:
	x=np.arange(max_n*2+2)
	y=prob[:x[-1]+1]
	plt.figure()
	plt.bar(x,y)
	plt.yscale('log')
	plt.ylabel('Probability %',fontsize=26)
	plt.xlabel('# Ions in the cluster',fontsize=26)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlim(0,max(x)+2)
if big:
	minimo,maximo=check_limits(numero)
	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob
	bax.bar(x,y)
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
	bax.set_yscale('log')


neutral=np.zeros((len(index_cation)+len(index_anion)+1))
positive=np.zeros((len(index_cation)+len(index_anion)+1))
negative=np.zeros((len(index_cation)+len(index_anion)+1))
for i in np.arange(len(total_clust_matrix)):
	for j in np.arange(len(total_clust_matrix[i])):
		if i==j:
			neutral[i+j]=neutral[i+j]+total_clust_matrix[i][j]
		if i>j:
			positive[i+j]=positive[i+j]+total_clust_matrix[i][j]
		if i<j:
			negative[i+j]=negative[i+j]+total_clust_matrix[i][j]
prob_neutral=100*neutral/np.sum(numero)
prob_positive=100*positive/np.sum(numero)
prob_negative=100*negative/np.sum(numero)
if not big:
	x=np.arange(max_n*2+2)
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	plt.figure()
	plt.bar(x,y+y_p+y_n,color='b',label='negative')
	plt.bar(x,y+y_p,color='r',label='positive')
	plt.bar(x,y,color='k',label='neutral')
#	plt.yscale('log')
	plt.ylabel('Probability %',fontsize=20)
	plt.xlabel('# Ions in the cluster',fontsize=20)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlim(0,max(x))
	plt.legend(fontsize=18)
if big:
	minimo,maximo=check_limits(numero)
	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	bax.bar(x,y+y_p+y_n,color='b',label='negative')
	bax.bar(x,y+y_p,color='r',label='positive')
	bax.bar(x,y,color='k',label='neutral')
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
#	bax.set_yscale('log')
	bax.legend(fontsize=18)

if big:
	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	bax.bar(x,y+y_p+y_n,color='b',label='negative')
	bax.bar(x,y+y_p,color='r',label='positive')
	bax.bar(x,y,color='k',label='neutral')
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
#	bax.set_yscale('log')
	bax.legend(fontsize=18)
	bax.set_ylim(0,2)

plt.show(block=False)
'''
numero=np.zeros((len(index_cation)+len(index_anion)+1))
for i in np.arange(len(total_clust_matrix)):
	for j in np.arange(len(total_clust_matrix[i])):
		numero[i+j]=numero[i+j]+total_clust_matrix[i][j]
if check=='y':
	g=0
	h=np.sum(numero)
	for i in np.arange(len(numero)):
		g=g+numero[i]*i
	print('Total # Ions = %i' %(g))
	print('Total # Clusters = %i' %(h))
	print('\n')

prob=np.zeros((len(numero)))
for i in np.arange(len(numero)):
	prob[i]=100*numero[i]*i/((len(numero)-1)*len(t))
print(np.sum(prob))
if not big:
	x=np.arange(max_n*2+2)
	y=prob[:x[-1]+1]
	plt.figure()
	plt.bar(x,y)
	plt.yscale('log')
	plt.ylabel('Probability %',fontsize=26)
	plt.xlabel('# Ions in the cluster',fontsize=26)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlim(0,max(x)+2)
if big:
	minimo,maximo=check_limits(numero)
	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob
	bax.bar(x,y)
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
	bax.set_yscale('log')

neutral=np.zeros((len(index_cation)+len(index_anion)+1))
positive=np.zeros((len(index_cation)+len(index_anion)+1))
negative=np.zeros((len(index_cation)+len(index_anion)+1))
for i in np.arange(len(total_clust_matrix)):
	for j in np.arange(len(total_clust_matrix[i])):
		if i==j:
			neutral[i+j]=neutral[i+j]+total_clust_matrix[i][j]
		if i>j:
			positive[i+j]=positive[i+j]+total_clust_matrix[i][j]
		if i<j:
			negative[i+j]=negative[i+j]+total_clust_matrix[i][j]
prob_neutral=np.zeros((len(neutral)))
prob_positive=np.zeros((len(positive)))
prob_negative=np.zeros((len(negative)))
for i in np.arange(len(neutral)):
	prob_neutral[i]=100*neutral[i]*i/((len(numero)-1)*len(t))
print(np.sum(prob_neutral))
for i in np.arange(len(positive)):
	prob_positive[i]=100*positive[i]*i/((len(numero)-1)*len(t))
print(np.sum(prob_positive))
for i in np.arange(len(negative)):
	prob_negative[i]=100*negative[i]*i/((len(numero)-1)*len(t))
print(np.sum(prob_negative))
print(np.sum(prob_neutral)+np.sum(prob_positive)+np.sum(prob_negative))
if not big:
	x=np.arange(max_n*2+2)
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	plt.figure()
	plt.bar(x,y+y_p+y_n,color='b',label='negative')
	plt.bar(x,y+y_p,color='r',label='positive')
	plt.bar(x,y,color='k',label='neutral')
	plt.ylabel('Probability %',fontsize=20)
	plt.xlabel('# Ions in the cluster',fontsize=20)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlim(0,max(x))
	plt.legend(fontsize=18)

	plt.figure()
	plt.bar(x,y+y_p+y_n,color='b',label='negative')
	plt.bar(x,y+y_p,color='r',label='positive')
	plt.bar(x,y,color='k',label='neutral')
	plt.ylabel('Probability %',fontsize=20)
	plt.xlabel('# Ions in the cluster',fontsize=20)
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	plt.xlim(0,max(x))
	plt.ylim(0,2)
	plt.legend(fontsize=18)
if big:
	minimo,maximo=check_limits(numero)
	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	bax.bar(x,y+y_p+y_n,color='b',label='negative')
	bax.bar(x,y+y_p,color='r',label='positive')
	bax.bar(x,y,color='k',label='neutral')
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
	bax.legend(fontsize=18)

	plt.figure()
	bax=brokenaxes(xlims=((0,minimo),(maximo,len(numero))),hspace=0.1,yscale='log')
	x=np.arange(len(numero))
	y=prob_neutral[:x[-1]+1]
	y_p=prob_positive[:x[-1]+1]
	y_n=prob_negative[:x[-1]+1]
	bax.bar(x,y+y_p+y_n,color='b',label='negative')
	bax.bar(x,y+y_p,color='r',label='positive')
	bax.bar(x,y,color='k',label='neutral')
	bax.set_ylabel('Probability %',fontsize=26)
	bax.set_xlabel('# Ions in the cluster',fontsize=26)
	bax.legend(fontsize=18)
	bax.set_ylim(0,2)

plt.show(block=False)