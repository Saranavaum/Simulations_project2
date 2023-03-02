##---- Initial conditions for Gadged 4 ----##

# Packages
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import h5py

#-- Initial condition parameters --#
# 3 Second Gas Sphere Collapse ,
# config = selfgravity , cubic_spline_kernel , Ntype =1, evalpotential

PartNum=32
totalPartNum=PartNum**3
totalMass=1
G=1
Side=2
InterEner=0.05
grid=Side/PartNum



# With this we generate the initial cube
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

# We calculate the distance to the center of each particle
dist= np.sqrt((XCoord)**2+(YCoord)**2+(ZCoord)**2)

# With this we are giving it the condition that it stays only with the particles that are
inside the sphere
R=Side/2
bol=dist<=R
x1=XCoord[bol]
x2=YCoord[bol]
x3=ZCoord[bol]
dist=dist[bol]

# Let’s redefine the sphere so that it has a profile of 1/r
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

# Now we calculate the density profile of the new profile , to see that it is ok
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

# Plotting the density profile
fig2= plt.figure()
ax2=plt.plot(shellRadius,shellDensity,'r.',label='Data')
plt.plot(rrho,rho,label='Analytical solution')
plt.xlabel('Radius')
plt.ylabel('Density')
plt.legend()
#plt.savefig('Density')
plt.show()

# We create the rest of the initial conditions to generate the hdf5

#-- Initial velocities --#
vx = np.zeros(n)
vy = np.zeros(n)
vz = np.zeros(n)

#-- Internal energy --#
Energy=np.zeros(n)
Energy+=0.05 

#-- Generation of HDF5 with the initial conditions established --#
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


