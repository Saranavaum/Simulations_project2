import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import h5py

#1)Problem Description
#2) Initial Conditinos: como se han generado
#3) Resultados: comparar los resultados con la analitica, despues de cprres ña simulacion, verr como se conserva la energia cuantificar esto

#apendices
#A) Compilation los parametros que hemos usados para la compilacion , el scri que se haya usado para hacer la gráfica
#B) Analysis code


#el colapso dura 3 segundos,
#esfera de gas, config= selfgravity, cubic_spline_kernel, Ntype=1, evalpotential
#archivo que se llama energia total integra la energia para todas las particulas
PartNum=32
totalPartNum=PartNum**3
totalMass=1
G=1
Side=2
InterEner=0.05
grid=Side/PartNum



#con esto generamos el cubo inicial
XCoord=np.array([])
YCoord=np.array([])
ZCoord=np.array([])

for i in range(PartNum):
	for j in range(PartNum):
		for k in range(PartNum):
			XCoord=np.append(XCoord, i*grid)
			YCoord=np.append(YCoord, j*grid)
			ZCoord=np.append(ZCoord, k*grid)

XCoord=XCoord-np.median(XCoord)
YCoord=YCoord-np.median(YCoord)
ZCoord=ZCoord-np.median(ZCoord)
#XCoord.shape = (totalPartNum, 1)
#YCoord.shape = (totalPartNum, 1)
#ZCoord.shape = (totalPartNum, 1)
#Coordinates = np.concatenate((XCoord,YCoord,ZCoord), axis = 1)
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.scatter(XCoord,YCoord,ZCoord, color="black", s=0.25)
#plt.show()


dist= np.sqrt((XCoord)**2+(YCoord)**2+(ZCoord)**2)

#con esto le estamos dando la condicion de quese quede solo con las particuls que estan dentro de la esfera
R=Side/2
bol=dist<=R
x1=XCoord[bol]
x2=YCoord[bol]
x3=ZCoord[bol]
dist=dist[bol]

#vamos a redefinir la esfera para que tenga un perfil de 1/r

x1=x1*np.sqrt(dist)
x2=x2*np.sqrt(dist)
x3=x3*np.sqrt(dist)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(x1,x2,x3, color="black", s=0.25)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
#plt.savefig('InitialCondition')

fig3= plt.figure()
ax3=plt.plot(x1,x3, 'k.', ms=2)
plt.xlabel('X')
plt.ylabel('Z')
plt.show()
#plt.savefig('InitialCondition2')

#ahora calculamos el perfil de densidad del nuevo perfil, para ver que está bien
dist= np.sqrt((x1)**2+(x2)**2+(x3)**2)

shellNumber = 20
shellWidth = np.mean(dist) / shellNumber
shellDensity = np.array([])
shellRadius = np.array([])
n=len(x1)
particleMass=totalMass/n

for i in range(shellNumber):
	numberOfParticles = 0
	for j in range(len(dist)):
		if ((dist[j] > (i * shellWidth)) and (dist[j] < ((i+1) * shellWidth))):
			numberOfParticles = numberOfParticles + 1
	shellRadius = np.append(shellRadius, shellWidth*i+(shellWidth)/2)
	massOfParticles = numberOfParticles * particleMass
	shellDensity = np.append(shellDensity, (massOfParticles)/((4./3)*np.pi*(((i+1)*shellWidth)**3 - (i*shellWidth)**3)))




def densi(M,R,r):
	rho=M/(2*np.pi*R**2*r)
	return(rho)
rrho=np.linspace(0,0.7,50)
rho=densi(totalMass,R,rrho)


fig2= plt.figure()
ax2=plt.plot(shellRadius,shellDensity,'r.',label='Data')
plt.plot(rrho,rho,label='Analytical solution')
plt.xlabel('Radius')
plt.ylabel('Density')
plt.legend()
#plt.savefig('Density')
plt.show()

#generamos las condiciones iniciales
#velocidades

vx = np.zeros(n)
vy = np.zeros(n)
vz = np.zeros(n)

##Energia interna
Energy=np.zeros(n)
Energy+=0.05 ### /particlemass
'''
hf = h5py.File("myInitialConditions.hdf5", 'w')

massvec=np.zeros(6)
massvec[0]=particleMass

number_part=np.zeros(6)
number_part[0]=n

ID = np.linspace(1, n, n, dtype = int)

header = hf.create_group('Header')
header.attrs['NumPart_ThisFile']    = number_part
header.attrs['MassTable']           = massvec
header.attrs['Time']                = 0
header.attrs['Redshift']            = 0
header.attrs['NumPart_Total']       = number_part
header.attrs['NumFilesPerSnapshot'] = 1
header.attrs['BoxSize']             = 1.0
header.attrs['Omega0']              = 1.0
header.attrs['OmegaLambdda']        = 0.
header.attrs['HubbleParam']         = 0.7
header.attrs['Flag_Entropy_ICs']    = 0
header.attrs['NumPart_Total_HighWord'] = np.zeros(6)



PartType0 = hf.create_group('PartType0')
PartType0.create_dataset("Coordinates", data = np.vstack([x1, x2, x3]).T)
PartType0.create_dataset("ParticleIDs", data = ID)
PartType0.create_dataset("Velocities", data = np.vstack([vx, vy, vz]).T)
ids_d=np.arange(n)
PartType0.create_dataset("InternalEnergy", data = Energy)

hf.close()

'''
